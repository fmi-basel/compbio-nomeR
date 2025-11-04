#' Calculate coverage probabilities for footprints in SummarizedExperiment object containing
#' single-molecule footprinting (SMF) data
#'
#'
#' @param se A \code{\link[SummarizedExperiment]{SummarizedExperiment}} object
#'     containing read-level data returned by \code{\link[footprintR]{readModBam}} function
#'     containing modification probabilities.
#' @param assayName Character scalar describing the name of tje assay in \code{se} containing
#'     read-level data.
#' @param threshUnmod,threshMod Numeric scalars used to classify observations
#'     as modified (modification probability >= threshMod), unmodified
#'     (modification probability < threshUnmod) or unknown (otherwise).
#' @param footprint_models A list containing footprint models for proteins.
#'     Each element must have 3 slots:
#'     \describe{
#'     \item{PROTECT_PROB}{a numeric vector with footprint emission
#'     probabilities to find protected position within a footprint}
#'     \item{COVER_PRIOR}{prior coverage probability (abundance)
#'     (\code{numeric}) reflecting what fraction of reads you expect to be
#'     covered by a footprint}
#'     \item{NAME}{name (\code{character}) of a model, e.g. "Nucleosome"}
#'     \item{GROUP}{non-unique group (\code{character}) which defines how probabilities will be aggregated 
#'     if \code{aggrByGroup} is \code{TRUE}. Namely, if "Nucleosome--149", "Nucleosome--150" etc. footprint models
#'     have identical GROUP (e.g. "Nucleosome") and \code{aggrByGroup = TRUE}, probabilities will be aggregated
#'     across all footprints with identical GROUP.}
#'     }
#' @param bgprotectprob background emission probability to find a protected
#'     position within open (accessible) regions.
#' @param bgcoverprior prior probability for percentage of all fragments to be
#'     in a free (accessible, or background) state.
#' @param aggrByGroup if \code{TRUE} probabilities are aggregated by GROUP ID defined
#'     in \code{footprint_models}. If \code{FALSE} or GROUP IDs are missing in the \code{footprint_models}
#'     probabilities are reported for each individual footprints NAME defined in the \code{footprint_models}.
#' @param report_prediction_in_flanks \code{logical} whether to return
#'     calculated start probabilities in left flanking region.
#'     In order to take into account partial footprints at left edge of
#'     fragments the algorithm extends each fragment by maximum footprint
#'     length on the left side. \code{report_prediction_in_flanks} controls
#'     whether calculated start probabilities in the left flanking region will
#'     be reported in the \code{START_PROB}.
#' @param ncpu number of threads to use.
#' @param verbose verbose mode for bug fixing.
#'
#' @return A list which contains 2 data frames:
#'     \describe{
#'     \item{\code{START_PROB}}{data frame with calculated start probabilities
#'     for each SMF molecule (column \code{seq}), each position in ROI (column
#'     \code{pos}) and each footprint model, e.g. Nucleosome, background etc.
#'     These probabilities reflect how likely it is to find a start in each
#'     fragment and at each position of a certain footprint model.}
#'     \item{\code{COVER_PROB}}{data frame with calculated coverage
#'     probabilities for each SMF molecule (column \code{seq}), each position
#'     in ROI (column \code{pos}) and each footprint model, e.g. Nucleosome,
#'     background etc. These probabilities reflect how likely it is that a
#'     certain position in an amplicon and certain fragment is covered by a
#'     certain footprint model.}
#'     }
#' @importFrom SummarizedExperiment SummarizedExperiment rowRanges colData colData<-
#' @importFrom SparseArray NaArray
#' @importFrom GenomicRanges GPos match seqnames start end strand seqinfo
#' @importFrom IRanges subsetByOverlaps IRanges IRangesList
#' @importFrom S4Vectors DataFrame SimpleList metadata metadata<- make_zero_col_DFrame
#' @import data.table
#' @export
#'
#' @examples
#' set.seed(3346)
#' nc <- 50
#' nr <- 50
#' rmatr <- matrix(data = as.integer(rnorm(nc * nr) >= 0.5),
#'                 ncol = nc,nrow=nr)
#'
#' ## create dummy footprints
#' bg.pr <- 0.5
#' ft.pr <- 1-bg.pr
#' ft.len <- 15
#'
#' ## creating a list of binding models for nomeR
#' ftp.models <- list(list("PROTECT_PROB" = rep(0.99,ft.len),
#'                         "COVER_PRIOR" = ft.pr,
#'                         "NAME" = "FOOTPRINT"))
#'
#' nomeR.out <- predict_footprints(data=rmatr,
#'                                 footprint_models = ftp.models,
#'                                 bgprotectprob = 0.05,
#'                                 bgcoverprior = bg.pr)
#'
#' @importFrom checkmate makeAssertCollection assert_logical assert_int
#'     reportAssertions
#' @importFrom parallel detectCores
predict_footprints_SE <- function(se,
																	assayName = "mod_prob",
																	threshUnmod = 0.5,
																	threshMod = 0.5,
																	footprint_models,
																	bgprotectprob,
																	bgcoverprior,
																	aggrByGroup = FALSE,
																	report_prediction_in_flanks = FALSE,
																	ncpu = 1L,
																	verbose = FALSE) {
	
	prob_group = fragpos = posidx_ref = fidx_glob = sidx = fidx_sample = chr = refpos = pos = mod_prob = gpos_idx = ftp_name = ftp_group = readName = sname = NULL # due to NSE notes in R CMD check
	
	## check arguments
	coll <- makeAssertCollection()
	### validate se object and prepare data for nomeR prediction
	protect_data <- validate_prepare_SE(se,
																			assayName,
																			threshUnmod,
																			threshMod)
	
	if(nrow(protect_data) == 0)
		stop("No data satisfy thresholds for modified and unmodified bases. Check parameters threshUnmod and threshMod")
	
	### validate footprint models
	ftpvalout <- validate_footprint_models(footprint_models,
																				 bgprotectprob,
																				 bgcoverprior,
																				 aggrByGroup,
																				 verbose,
																				 add = coll)
	footprint_models <- ftpvalout[["footprint_models"]]
	start_priors <- ftpvalout[["start_priors"]]
	### validate report_prediction_in_flanks
	assert_logical(report_prediction_in_flanks,
								 any.missing = FALSE, all.missing = FALSE,
								 len = 1, add = coll)
	
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
	
	## finish argument check
	reportAssertions(coll)
	
	if (verbose) {
		.message_timestamp("Calling run_cpp_nomeR...")
	}
	
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
	
	predict_res_list <- calcStartCoverProbs_cpp(protect_data[["fidx_glob"]], ## unique fragment ID or index
																							protect_data[["fragpos"]],      ## position within fragment, 1 - based
																							protect_data[["protect"]],   ## binary protection data, 0 - accessible, 1 - protected
																							footprint_models,
																							bgprotectprob,
																							start_priors["BG"],
																							report_prediction_in_flanks,
																							ncpu,
																							verbose)
	
	
	
	## construct ouput SE
	if (all(c(!is.null(predict_res_list[["START_PROB"]]),
						!is.null(predict_res_list[["COVER_PROB"]]),
						!is.null(predict_res_list[["VITERBI_CONF"]])))) {
		if (verbose) {
			.message_timestamp("Constructing output SummarizedExperiment... ")
		}
		
		## convert to data.table and rbind
		predict_res <- rbindlist(lapply(c("START_PROB",
																			"COVER_PROB"),
																		function(nm){
																			x <- as.data.table(predict_res_list[[nm]])
																			x <- x[,prob_group := nm]
																			return(x)
																		}))
		
		
		## create annotation of reads
		rowGpos <- rowRanges(se)
		
		frag2sample_anno <- protect_data[fragpos == 1][,
																									 c("strand","chr") := list(as.character(strand(rowRanges(se))[posidx_ref]),
																									 											 as.character(seqnames(rowRanges(se))[posidx_ref]))][,
																									 											 																										list(fidx_glob,
																									 											 																											sidx,
																									 											 																											fidx_sample,
																									 											 																											posidx_ref,
																									 											 																											chr,
																									 											 																											refpos,
																									 											 																											strand)]
		## add readNames
		mod_prob_assays <- assay(se,"mod_prob")
		readNames <- data.table::rbindlist(lapply(1:ncol(mod_prob_assays),
																							function(sidx){
																								data.table(sidx = sidx,
																													 fidx_sample = 1:ncol(mod_prob_assays[[sidx]]),
																													 readName = colnames(mod_prob_assays[[sidx]]))
																							}))
		frag2sample_anno <- readNames[frag2sample_anno,on = list(sidx == sidx,fidx_sample == fidx_sample)]
		## add reference positions
		predict_res <- predict_res[,refpos := pos - 1 + frag2sample_anno[match(seq,frag2sample_anno[["fidx_glob"]])][["refpos"]]]
		
		## add chr, strand, sidx, and fidx_sample
		predict_res <- frag2sample_anno[,list(fidx_glob,sidx, fidx_sample,chr,strand)][predict_res,
																																								on = list(fidx_glob = seq)]
		
		## add modprob 
		predict_res <- protect_data[,list(fidx_glob,fragpos,mod_prob)][predict_res, on = list(fidx_glob = fidx_glob,
																																										fragpos = pos)]
		fcols <- c("prob_group","fidx_glob","sidx","fidx_sample","fragpos",
							 "chr","refpos","strand",
							 "mod_prob")
		setcolorder(predict_res,c(fcols,
															setdiff(colnames(predict_res),fcols)))
		
		## create rowRanges
		posuniq <- unique(predict_res[,list(chr,refpos,strand)])[,gpos_idx := 1:.N]
		## add gposidx
		predict_res <- posuniq[predict_res,
													 on = list(chr=chr,refpos=refpos,strand=strand)]
		seOutRowRanges <- GenomicRanges::GPos(seqnames = posuniq[["chr"]],
																					pos = posuniq[["refpos"]],
																					strand = posuniq[["strand"]],
																					seqinfo = GenomicRanges::seqinfo(rowGpos))
		
		
		ftpnames <- setdiff(colnames(predict_res),c(fcols,"gpos_idx"))
		nomeR_assayNames <- c("mod_prob",paste(rep(ftpnames,2),
																					 rep(c("coverProb","startProb"),each = length(ftpnames)),
																					 "nomeR",
																					 sep="_"))
		assayAnno <- data.frame(assayName = nomeR_assayNames,
														ftpName = c("mod_prob",rep(ftpnames,2)),
														probName = c("START_PROB",rep(c("COVER_PROB","START_PROB"),each = length(ftpnames))))
		
		## extract readNames
		
		## create list of assays
		assayList <- lapply(1:nrow(assayAnno),
												function(assayI){
													assayMat <- make_zero_col_DFrame(nrow = length(seOutRowRanges))
													for(sI in 1:ncol(se)){
														
														### select which fragments belong to current sample. 
														curDat <- predict_res[sidx == sI & prob_group == assayAnno$probName[assayI]]
														curDat <- curDat[!is.na(curDat[[assayAnno$ftpName[assayI]]])]
														maxFidx <- frag2sample_anno[sidx == sI,max(fidx_sample)]
														curFragNames <- frag2sample_anno[sidx == sI][match(1:maxFidx,fidx_sample)][["readName"]]
														## get read names
														
														namat <- NaArray(dim = c(length(seOutRowRanges), maxFidx),
																						 dimnames = list(NULL,curFragNames),
																						 type = "double")
														## add data
														
														namat[as.matrix(curDat[,list(gpos_idx,fidx_sample)])] <- curDat[[assayAnno$ftpName[assayI]]]
														assayMat[[assayAnno$assayName[assayI]]] <- namat
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
		
		
		#browser()
		## construct IRangesLists with MAP configurations and add to colData
		viterbi_conf <- as.data.table(predict_res_list[["VITERBI_CONF"]])
		## ignore background
		viterbi_conf <- viterbi_conf[ftp_name != "background"]
		## add reference positions
		viterbi_conf <- viterbi_conf[,refpos := start - 1 + frag2sample_anno[match(seq,frag2sample_anno[["fidx_glob"]])][["refpos"]]]
		## add sidx, fidx_sample, readName
		viterbi_conf <- frag2sample_anno[,list(fidx_glob,sidx, fidx_sample,readName)][viterbi_conf,
																																								on = list(fidx_glob = seq)]
		coldat <- colData(se)
		## add sample names
		viterbi_conf <- viterbi_conf[,sname := coldat$sample[sidx]]
		
		if(aggrByGroup){
			vit_ftpnames <- unique(viterbi_conf[["ftp_group"]])
		} else{
			vit_ftpnames <- unique(viterbi_conf[["ftp_name"]])
		}
		for(ftp in vit_ftpnames){
			
			
			lIRl <- sapply(coldat$sample,
										 function(snm){
										 	if(aggrByGroup){
										 		ftpLoc <- viterbi_conf[ftp_group == ftp & sname == snm]
										 	} else{
										 		ftpLoc <- viterbi_conf[ftp_name == ftp & sname == snm]
										 	}
										 	
										 	irL <- IRanges(start = ftpLoc[["refpos"]],
										 								 width = ftpLoc[["width"]],
										 								 ftp_name = ftpLoc[["ftp_name"]],
										 								 ftp_group = ftpLoc[["ftp_group"]],
										 								 start_prob = ftpLoc[["start_prob"]])
										 	irL <- IRangesList(split(irL,ftpLoc[["readName"]]))
										 	return(irL)
										 },simplify=F,USE.NAMES = T)
			
			## remove "--" for colnames and add nomeR
			ftp_colnm <- paste0(gsub("-","_",ftp),"_nomeR")
			coldat[[ftp_colnm]] <- lIRl
		}
		
		colData(seOut) <- coldat
		
		## change metadata
		mtdat <- metadata(seOut)
		## add readLevelData assayNames
		mtdat$readLevelData$assayNames <- assayNames(seOut)
		mtdat$readLevelData$colDataColumns <- c(mtdat$readLevelData$colDataColumns,
																						paste0(vit_ftpnames,"_nomeR"))
		
		metadata(seOut) <- mtdat
		return(seOut)
	} else {
		stop("retrieved NULL results from the C++ function.")
	}
}
