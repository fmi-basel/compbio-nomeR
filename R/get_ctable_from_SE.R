#' Create count table for 00, 01, 10, 11, etc occurrences in SummarizedExperiment object
#' containing modificaton probabilities for SMF data
#'
#' @param se A \code{\link[SummarizedExperiment]{SummarizedExperiment}} object
#'     containing read-level data returned by \code{\link[footprintR]{readModBam}} function
#'     containing modification probabilities.
#' @param assayName Character scalar describing the name of tje assay in \code{se} containing
#'     read-level data.
#' @param threshUnmod,threshMod Numeric scalars used to classify observations
#'     as modified (modification probability >= threshMod), unmodified
#'     (modification probability < threshUnmod) or unknown (otherwise).
#' @param max_spacing Maximum spacing between positions.
#' @param aggrSamples \code{logical} return frequencies for each sample separately or aggregate them for all samples.
#' @param ncpu number of cores to use.
#' @param verbose verbose mode for bug fixing.
#'
#' @return if \code{aggrSamples = T} \code{matrix} containing aggregated frequencies across all samples for
#'     0,0; 0,1; etc for each spacing between 1 (total number of 0s and 1s in the data) and \code{max_spacing}.
#'     If \code{aggrSamples = F} a \code{list} of matrices with frequences for each sample separately.
#'
#'
#' @export
#'
get_ctable_from_SE <- function(se,
															 assayName = "mod_prob",
															 threshUnmod = 0.5,
															 threshMod = 0.5,
															 max_spacing = 200L,
															 aggrSamples = F,
															 ncpu = 1L,
															 verbose = FALSE) {

	### validate se object and prepare data for nomeR prediction


	protect_data <- validate_prepare_SE(se,
																			assayName,
																			threshUnmod,
																			threshMod)

	## collect pair stats for each sample
	ctable_list <- lapply(1:ncol(se),
												function(csidx){
													count_spacing_freq_cpp(protect_data[sidx == csidx][["fidx_glob"]], ## unique fragment ID or index
																								 protect_data[sidx == csidx][["fragpos"]],      ## position within fragment, 1 - based
																								 protect_data[sidx == csidx][["protect"]],   ## binary protection data, 0 - accessible, 1 - protected
																								 max_spacing,
																								 ncpu,
																								 verbose)
												})
	names(ctable_list) <- colnames(se)

	if(aggrSamples){
		return(Reduce("+",ctable_list))
	}
	return(ctable_list)
}
