#' Calculate coverage probabilities for footprints in SummarizedExperiment object containing
#' single-molecule footprinting (SMF) data
#'
#'
#' @param se A \code{\link[SummarizedExperiment]{SummarizedExperiment}} object
#'     containing read-level data returned by \code{\link[readModBam]{footprintR}} function
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
#' @importFrom SummarizedExperiment SummarizedExperiment rowRanges colData
#' @importFrom SparseArray NaArray
#' @importFrom GenomicRanges GPos match seqnames start end
#' @importFrom IRanges subsetByOverlaps
#' @importFrom S4Vectors DataFrame SimpleList
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
	browser()
	## protect_data is a matrix returned by validate_prepare_SE
	## columns are:
	## sidx - index of sample in SE
	## fidx_glob - unique index of fragment across all samples, as if they were cbinded
	## fidx_sample - index of fragment for the current sample
	## posidx_ref - index of rows in SE, corresponds to reference position stored in rowRanges(se)
	## protect - binary protection data, 0 - accessible, 1 - protected
	## refpos - genomic position within a reference
	## fragpos - position within a frament, 1 - based

	## the C++ needs only fidx_glob, fragpos, protect
	out.list <- calcStartCoverProbs_cpp(protect_data[,"fidx_glob"], ## unique fragment ID or index
																			protect_data[,"fragpos"],      ## position within fragment, 1 - based
																			protect_data[,"protect"],   ## binary protection data, 0 - accessible, 1 - protected
																			footprint_models,
																			bgprotectprob,
																			start_priors["BG"],
																			report_prediction_in_flanks,
																			ncpu,
																			verbose)
	
	
	

	### TODO: the calcStartCoverProbs_cpp will return calculated probabilities, START_PROB and COVER_PROB for EACH position within a fragment
	### Convert this output to a format that could be added as assay into SE input object
	

	if (all(c(!is.null(out.list[["START_PROB"]]),
						!is.null(out.list[["COVER_PROB"]])))) {
		if (verbose) {
			.message_timestamp("convert cpp_nomeR output to data.frame...")
		}
		return(lapply(out.list,as.data.frame,
									stringsAsFactors = FALSE,
									check.names = FALSE))
	} else {
		stop("retrieved NULL results from C++ function.")
	}
}
