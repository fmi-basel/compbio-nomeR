#' Calculate posterior probabilities and predict footprints in a
#' SummarizedExperiment containing single-molecule footprinting (SMF) data
#'
#' @description
#' Computes posterior probabilities and predicts footprint configurations for
#' SMF data stored in a \code{SummarizedExperiment} object.
#'
#' @param se A \code{\link[SummarizedExperiment]{SummarizedExperiment}} object
#'   containing read-level data returned by
#'   \code{footprintR::readModBam}, including modification probabilities.
#' @param assayName Character scalar specifying the name of the assay in
#'   \code{se} that contains read-level modification probabilities.
#' @param threshUnmod,threshMod Numeric thresholds used to binarize modification
#'   probabilities into accessible (\code{0}; probability >= \code{threshMod}),
#'   protected (\code{1}; probability < \code{threshUnmod}), or unknown
#'   (\code{NA}) states.
#' @param min_frag_data_len Ignore fragments whose genomic span (from the
#'   leftmost to rightmost non-\code{NA} data point) is shorter than
#'   \code{min_frag_data_len}.
#' @param min_frag_data_dens Ignore fragments with a density of informative
#'   (non-\code{NA}) positions below \code{min_frag_data_dens}.
#' @param profile Enable time profiling.
#' @inheritParams predict_footprints
#'
#' @return A \code{SummarizedExperiment} object containing:
#'   \itemize{
#'     \item the original \code{mod_prob} assay,
#'     \item additional assays with calculated posterior start and coverage
#'       probabilities (e.g. "Nucl_coverProb_nomeR"), and
#'     \item predicted footprint configurations stored as \code{IntegerList}
#'     objects in \code{colData} (e.g. column "Nucl_nomeR").
#'   }
#'
#' @importFrom SummarizedExperiment SummarizedExperiment rowRanges colData
#'     colData<-
#' @importFrom SparseArray NaArray
#' @importFrom GenomicRanges GPos match seqnames start end strand
#' @importFrom Seqinfo seqinfo
#' @importFrom IRanges subsetByOverlaps IRanges IRangesList
#' @importFrom S4Vectors DataFrame SimpleList metadata metadata<-
#'     make_zero_col_DFrame
#' @import data.table
#' @importFrom checkmate makeAssertCollection assert_logical assert_int
#'     reportAssertions
#' @importFrom parallel detectCores
#'
#' @references
#' Fariselli, P., Martelli, P. L., & Casadio, R. (2005).
#' *A new decoding algorithm for Hidden Markov Models improves the prediction of the topology of all-beta membrane proteins.*
#' BMC Bioinformatics, 6(S4), S12. https://doi.org/10.1186/1471-2105-6-S4-S12
#'
#' @export

predict_footprints_SE <- function(se,
                                  assayName = "mod_prob",
                                  threshMod = 0.5,
                                  threshUnmod = threshMod,
                                  min_frag_data_len = 50L,
                                  min_frag_data_dens = 0.05,
                                  footprint_models,
                                  bgprotectprob,
                                  bgcoverprior,
                                  aggrByGroup = TRUE,
                                  ftpConfigMethod = c("PV","PosteriorDecoding", "Viterbi"),

                                  ncpu = 1L,
                                  verbose = FALSE,
                                  profile = FALSE) {

    prob_group <- fragpos <- posidx_ref <- fidx_glob <- sidx <- fidx_sample <- chr <-
        refpos <- pos <- mod_prob <- gpos_idx <- ftp_name <- ftp_group <- readName <-
        sname <- NULL # due to NSE notes in R CMD check

    ftpConfigMethod <- match.arg(ftpConfigMethod)
    ### validate se object and prepare data for nomeR prediction
    dataList <- validate_prepare_SE(se,
                                    assayName,
                                    threshMod = threshMod,
                                    threshUnmod = threshUnmod,
                                    min_frag_data_len,
                                    min_frag_data_dens)

    protect_data <- dataList[["bin_protect_data"]]
    fragAnno <- dataList[["fragAnno"]]

    ### validate footprint models
    ftpvalout <- validate_footprint_models(footprint_models,
                                           bgprotectprob,
                                           bgcoverprior,
                                           aggrByGroup,
                                           verbose,
                                           add = coll)
    footprint_models <- ftpvalout[["footprint_models"]]
    start_priors <- ftpvalout[["start_priors"]]

    ### validate ncpu
    assert_int(x = ncpu, lower = 0, na.ok = TRUE, add = coll)
    avail_ncpu <- parallel::detectCores()
    if (is.na(avail_ncpu)) {
        .warning_timestamp(
            "Could not detect number of available cpu. Setting ncpu to 1L.")
        ncpu <- 1L
    } else if (ncpu > avail_ncpu || ncpu == 0) {
        .warning_timestamp(c("Number of ncpu is 0 or exceeds number of ",
                             "available cpu. Setting ncpu to number of ",
                             "available cpus."))
        ncpu <- avail_ncpu
    }

    ## restrict data.table to use only ncpu threads
    setDTthreads(threads = ncpu)
    ## protect_data is a matrix returned by validate_prepare_SE
    ## columns are:
    ## sidx - index of sample in SE
    ## fidx_glob - unique index of fragment across all samples, as if they were cbinded
    ## fidx_sample - index of fragment for the current sample
    ## posidx_ref - index of rows in SE, corresponds to reference position stored in rowRanges(se)
    ## protect - binary protection data, 0 - accessible, 1 - protected
    ## refpos - genomic position within a reference
    ## fragpos - position within a frament, 1 - based

    ## the calcStartCoverProbs_cpp needs only fidx_glob, fragpos, protect
    if (verbose) {
        .message_timestamp("Footprint prediction... ")
    }

    profiling <- .is_profiling_enabled(profile)
    timings <- new.env(parent = emptyenv())

    if(profiling)
        cli::cli_h1("predict_footprints_SE time-profiling", .envir = if (profiling) parent.frame() else NULL)

    predict_res_list <- .time_block({
        calcStartCoverProbs_cpp(
            protect_data[["fidx_glob"]], ## unique fragment ID or index
            protect_data[["fragpos"]],      ## position within fragment, 1 - based
            protect_data[["protect"]],   ## binary protection data, 0 - accessible, 1 - protected
            footprint_models,
            bgprotectprob,
            start_priors["BG"],
            ftpConfigMethod,
            aggrByGroup,
            ncpu,
            verbose)
    }, "Step 1: Calling C++ for posterior calculations ", timings, profiling)



    ## construct ouput SE
    if (all(c(!is.null(predict_res_list[["START_PROB"]]),
              !is.null(predict_res_list[["COVER_PROB"]]),
              !is.null(predict_res_list[["FOOTPRINT_CONF"]])))) {

        seOut <- .time_block({

            if (verbose) {
                .message_timestamp("Constructing output SummarizedExperiment... ")
            }

            ## convert to data.table

            ## positions are defined by START_PROB, because they run from
            ## -maxPWMlen..lastDatPos
            ## COVER_PROB run from firstDatPos...lastDatPos
            ## order of columns
            fcols <- c("fidx_glob", "sidx", "fidx_sample", "fragpos",
                       "chr", "refpos", "strand")

            predict_res <- sapply(
                c("START_PROB", "COVER_PROB"),
                function(nm) {

                    ## convert to data.table
                    prob_dt <- as.data.table(predict_res_list[[nm]])
                    ## add reference positions
                    prob_dt <- prob_dt[, refpos := pos - 1 + fragAnno[match(seq, fragAnno[["fidx_glob"]])][["refStart"]]]

                    ## add chr, strand, sidx, and fidx_sample
                    prob_dt <- fragAnno[, list(fidx_glob,sidx,
                                               fidx_sample, chr,
                                               strand)][prob_dt,
                                                        on = list(fidx_glob = seq)]
                    setnames(prob_dt, "pos", "fragpos")

                    setcolorder(prob_dt, c(fcols,
                                           setdiff(colnames(prob_dt), fcols)))
                }, simplify = FALSE, USE.NAMES = TRUE)

            ## add modprob to COVER_PROB, as it runs from firstDatPos to lastDatPos
            predict_res[["COVER_PROB"]] <-
                protect_data[, list(fidx_glob, fragpos,
                                    mod_prob)][predict_res[["COVER_PROB"]],
                                               on = list(fidx_glob = fidx_glob,
                                                         fragpos = fragpos)]
            ## create rowRanges
            posuniq <- unique(
                predict_res[["START_PROB"]][, list(chr, refpos,
                                                   strand)])[, gpos_idx := 1:.N]
            ## add gposidx to START_PROB and COVER_PROB
            predict_res[["START_PROB"]] <- posuniq[predict_res[["START_PROB"]],
                                                   on = list(chr = chr,
                                                             refpos = refpos,
                                                             strand = strand)]
            predict_res[["COVER_PROB"]] <- posuniq[predict_res[["COVER_PROB"]],
                                                   on = list(chr = chr,
                                                             refpos = refpos,
                                                             strand = strand)]

            seOutRowRanges <- GenomicRanges::GPos(seqnames = posuniq[["chr"]],
                                                  pos = posuniq[["refpos"]],
                                                  strand = posuniq[["strand"]],
                                                  seqinfo = Seqinfo::seqinfo(se))

            ftpnames <- setdiff(colnames(predict_res[["COVER_PROB"]]),
                                c(fcols, "mod_prob", "gpos_idx"))
            nomeR_assayNames <- c("mod_prob", paste(rep(ftpnames, 2),
                                                    rep(c("coverProb", "startProb"),
                                                        each = length(ftpnames)),
                                                    "nomeR",
                                                    sep = "_"))
            assayAnno <- data.frame(assayName = nomeR_assayNames,
                                    ftpName = c("mod_prob", rep(ftpnames, 2)),
                                    probName = c("COVER_PROB",
                                                 rep(c("COVER_PROB", "START_PROB"),
                                                     each = length(ftpnames))))

            ## background_startProb and background_coverProb are identical.
            ## keep only coverProb
            assayAnno <- assayAnno[assayAnno$assayName !=
                                       "background_startProb_nomeR", , drop = FALSE]

            ## create list of assays
            assayList <- lapply(
                seq_len(nrow(assayAnno)),
                function(assayI) {
                    assayMat <- make_zero_col_DFrame(nrow = length(seOutRowRanges))
                    curProbName <- assayAnno$probName[assayI]
                    curFtpName <- assayAnno$ftpName[assayI]
                    for (sI in seq_len(ncol(se))) {

                        ## get required data
                        curDat <-
                            predict_res[[curProbName]][, .SD,
                                                       .SDcols = c("fidx_glob",
                                                                   "sidx", "gpos_idx",
                                                                   curFtpName)][
                                                                       sidx == sI &
                                                                           !is.na(get(curFtpName))]

                        ### select which fragments that belong to current sample
                        ### and keep only those which passed the filtering
                        curSmpFrags <- fragAnno[sidx == sI & keep][, curFragIdx := 1:.N]
                        curDat <- curDat[, curFragIdx := curSmpFrags[match(
                            curDat$fidx_glob, curSmpFrags$fidx_glob), "curFragIdx"]]

                        namat <- NaArray(dim = c(length(seOutRowRanges),
                                                 nrow(curSmpFrags)),
                                         dimnames = list(NULL, curSmpFrags$readName),
                                         type = "double")

                        ## add data
                        namat[as.matrix(curDat[, list(gpos_idx, curFragIdx)])] <-
                            curDat[[curFtpName]]
                        assayMat[[sI]] <- namat
                    }
                    colnames(assayMat) <- colnames(se)
                    return(assayMat)
                })
            names(assayList) <- assayAnno$assayName

            seOut <- SummarizedExperiment(
                assays = assayList,
                rowRanges = seOutRowRanges,
                colData = colData(se),
                metadata = metadata(se)
            )

            ## construct IRangesLists with footprint configurations and add to colData
            footprint_conf <- as.data.table(predict_res_list[["FOOTPRINT_CONF"]])

            ## ignore background
            #footprint_conf <- footprint_conf[ftp_name != "background"]
            ## add reference positions
            footprint_conf <-
                footprint_conf[, refpos := start - 1 +
                                   fragAnno[match(seq,
                                                  fragAnno[["fidx_glob"]])][["refStart"]]]

            ## add sidx, fidx_sample, readName
            footprint_conf <- fragAnno[, list(fidx_glob,sidx,
                                              fidx_sample,
                                              readName)][footprint_conf,
                                                         on = list(fidx_glob = seq)]

            ## add sample names
            coldat <- colData(se)
            footprint_conf <- footprint_conf[, sname := coldat$sample[sidx]]


            if (aggrByGroup) {
                ftpConf_ftpnames <- unique(footprint_conf[["ftp_group"]])
            } else {
                ftpConf_ftpnames <- unique(footprint_conf[["ftp_name"]])
            }
            for (ftp in ftpConf_ftpnames) {
                lIRl <- sapply(
                    coldat$sample,
                    function(snm) {
                        if (aggrByGroup) {
                            ftpLoc <- footprint_conf[ftp_group == ftp & sname == snm]
                        } else{
                            ftpLoc <- footprint_conf[ftp_name == ftp & sname == snm]
                        }

                        irL <- IRanges(start = ftpLoc[["refpos"]],
                                       width = ftpLoc[["width"]],
                                       ftp_name = ftpLoc[["ftp_name"]],
                                       ftp_group = ftpLoc[["ftp_group"]],
                                       score = ftpLoc[["score"]])
                        irL <- IRangesList(split(irL,ftpLoc[["readName"]]))
                        return(irL)
                    }, simplify = FALSE, USE.NAMES = TRUE)

                ## remove "--" for colnames and add nomeR
                ftp_colnm <- paste0(gsub("-", "_", ftp), "_nomeR")
                coldat[[ftp_colnm]] <- lIRl
            }

            colData(seOut) <- coldat

            ## change metadata
            mtdat <- metadata(seOut)
            ## add readLevelData assayNames
            mtdat$readLevelData$assayNames <- assayNames(seOut)
            mtdat$readLevelData$colDataColumns <- c(mtdat$readLevelData$colDataColumns,
                                                    paste0(ftpConf_ftpnames, "_nomeR"))

            metadata(seOut) <- mtdat
            seOut
        }, "Step 2: Constructing output SE ", timings, profiling)

        ## add timings into metadata
        if(profiling){
            mtdat <- metadata(seOut)
            mtdat$timings <- as.list(timings)
            metadata(seOut) <- mtdat
        }
        return(seOut)
    } else {
        stop("retrieved NULL results from the C++ function.")
    }

}
