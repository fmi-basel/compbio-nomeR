#' Checks SummarizedExperiment and prepares data structure for the
#' nomeR C++ prediction functions
#'
#' @param data \code{matrix} or \code{list} with SMF data
#' @param assayName Character scalar describing the name of the assay in
#'     \code{se} containing read-level data.
#' @param min_frag_data_len \code{integer} ignore fragments that have genomic
#'     lengths from most-left to most-right data points less than
#'     \code{min_frag_data_len}.
#' @param min_frag_data_dens \code{numeric} ignore fragments that have density of
#'     data-containing positions lower than \code{min_frag_data_dens}.
#'
#' @return \code{list} with slots - mod_prob_data and fragAnno
#'
#' @keywords internal
#' @noRd
#' @importFrom SummarizedExperiment SummarizedExperiment rowRanges colData
#'     assay assayNames
#' @importFrom SparseArray NaArray nnawhich
#' @importFrom GenomicRanges GPos match seqnames start end
#' @importFrom IRanges subsetByOverlaps
#' @importFrom S4Vectors DataFrame SimpleList metadata
#' @import data.table
validate_prepare_SE <- function(se,
                                assayName = "mod_prob",
                                min_frag_data_len = 50L,
                                min_frag_data_dens = 0.05) {

    #mod_prob = fidx_sample = posidx_ref = refpos = fidx_glob = ftp_group = NULL # due to NSE notes in R CMD check

    ### The code for checking the vailidity of se is copied
    ### from the footprintR package developed by Charlotte Soneson and Michael Stadler

    stopifnot(is(se, "SummarizedExperiment"))

    stopifnot(!is.null(assayNames(se)) &&
                  all(assayNames(se) != "") &&
                  !any(duplicated(assayNames(se))))

    if (nrow(se) > 0) {
        if (!is.null(rownames(se))) {
            stopifnot(anyDuplicated(rownames(se)) == 0L)
        }
    }

    # if (nrow(se) > 0) {
    #     stopifnot(!is.null(rownames(se)) &&
    #                   !any(duplicated(rownames(se))))
    # }

    stopifnot(!is.null(metadata(se)$readLevelData) &&
                  is.list(metadata(se)$readLevelData) &&
                  all(c("assayNames", "colDataColumns") %in%
                          names(metadata(se)$readLevelData)))

    stopifnot("sample" %in% colnames(colData(se)))

    stopifnot(assayName %in% assayNames(se))

    ## extract modification probabilities from assay

    mod_prob_assays <- assay(se, assayName)
    fidx_glob_offset <- cumsum(vapply(mod_prob_assays, ncol, 0L))
    fidx_glob_offset <- c(0, fidx_glob_offset)

    ## annotation of fragments in input SE
    fragAnno <- data.table::rbindlist(
        lapply(seq_len(ncol(mod_prob_assays)),
               function(sidx) {
                   curanno <- data.table(
                       sidx = sidx,
                       fidx_sample = seq_len(ncol(mod_prob_assays[[sidx]])),
                       fragID = colnames(mod_prob_assays[[sidx]]),
                       fidx_glob = fidx_glob_offset[sidx] +
                           seq_len(ncol(mod_prob_assays[[sidx]])))
               }))


    ## collect modification probabilities for all non-NA positions

    mod_prob_data <- data.table::rbindlist(
        lapply(seq_len(ncol(mod_prob_assays)),
               function(sidx) {

                   read_naar <- mod_prob_assays[[sidx]]
                   ## get M-indices of non-NAs
                   nonNA_data <- nnawhich(read_naar, arr.ind = TRUE)
                   colnames(nonNA_data) <- c("posidx_ref", "fidx_sample")
                   ## convert to data.table
                   nonNA_data <- as.data.table(nonNA_data)
                   ## first column - positions (rows), second column - reads(columns)

                   ## add modification probability
                   nonNA_data <-
                       nonNA_data[, "mod_prob" := read_naar[as.matrix(nonNA_data)]]

                   ## construct output
                   ## sidx - index of sample in SE
                   ## fidx_glob - unique index of fragment across all samples, as if they were cbinded
                   ## fidx_sample - index of fragment for the current sample
                   ## posidx_ref - index of rows in SE, corresponds to reference position stored in rowRanges(se)
                   ## mod_prob - modification probability in [0,1]; high values indicate accessible positions
                   nonNA_data <- nonNA_data[, c("sidx", "fidx_glob") := list(
                       rep(sidx, nrow(nonNA_data)),
                       fidx_glob_offset[sidx] + fidx_sample)]

                   return(nonNA_data)
               }))

    rowGpos <- rowRanges(se)
    fragSummary <- mod_prob_data[,
                                    list("dataNpoints" = .N,
                                         "minPosIdx_ref" = min(posidx_ref),
                                         "maxPosIdx_ref" = max(posidx_ref)
                                    ),
                                    by = .(fidx_glob)]


    fragSummary <- fragSummary[, c(
        "chr",
        "strand",
        "refStart",
        "refEnd") := list(as.character(seqnames(rowGpos)[minPosIdx_ref]),
                          as.character(strand(rowGpos)[minPosIdx_ref]),
                          start(rowGpos)[minPosIdx_ref],
                          end(rowGpos)[maxPosIdx_ref])]
    fragSummary <-
        fragSummary[,"data_len" := refEnd -
                        refStart + 1][, "data_dens" := dataNpoints/data_len]
    fragAnno <- fragSummary[fragAnno, on = c(fidx_glob = "fidx_glob")]

    # ## add reference position
    mod_prob_data <-
        mod_prob_data[, "refpos" := start(rowGpos)[posidx_ref]]

    ## add position within fragments
    ## NOTE: the fragpos are 1 - based positions within fragments
    mod_prob_data <- mod_prob_data[,
                                         "fragpos" := refpos - fragAnno[["refStart"]][match(fidx_glob,fragAnno[["fidx_glob"]])] + 1]

    ## the below fails on MacOS
    # mod_prob_data <-
    #     mod_prob_data[, "fragpos" := refpos - min(refpos) + 1,
    #                      by = fidx_glob]



    ## order by fidx_glob and fragpos by setting keyv
    setkeyv(mod_prob_data, cols = c("fidx_glob", "fragpos"))

    ## sidx - index of sample in SE
    ## fidx_glob - unique index of fragment across all samples, as if they were cbinded
    ## fidx_sample - index of fragment for the current sample
    ## posidx_ref - index of rows in SE, corresponds to reference position stored in rowRanges(se)
    ## mod_prob - modification probability in [0,1]; high values indicate accessible positions
    ## refpos - genomic position within a reference
    ## fragpos - position within a frament, 1 - based

    ## reorder columns
    setcolorder(mod_prob_data,
                c("sidx", "fidx_glob", "fidx_sample", "posidx_ref",
                  "refpos", "fragpos", "mod_prob"))

    ## filter fragments by min_frag_data_len and min_frag_data_dens
    fragAnno <- fragAnno[,"keep" := !is.na(data_len) &
                             (data_len >= min_frag_data_len &
                                  data_dens >= min_frag_data_dens)]
    fragIDkeep <- fragAnno[keep == TRUE][["fidx_glob"]]
    mod_prob_data <- mod_prob_data[fidx_glob %in% fragIDkeep]

    if (nrow(mod_prob_data) == 0) {
        stop("No fragments left after filtering by frag_data_len>=",
             min_frag_data_len, "; frag_data_dens>=", min_frag_data_dens)
    }
    Nremove <- fragAnno[, sum(!keep), ]
    if (Nremove > 0) {
        .warning_timestamp(paste0(
            Nremove,
            " fragments have been removed after filtering by frag_data_len>=",
            min_frag_data_len, "; frag_data_dens>=", min_frag_data_dens))
    }
    ## order by fidx_glob and fragpos by setting keyv
    setkeyv(fragAnno, cols = c("fidx_glob"))

    return(list("mod_prob_data" = mod_prob_data,
                "fragAnno" = fragAnno))
}
