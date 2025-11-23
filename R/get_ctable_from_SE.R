#' Create count table for 00, 01, 10, 11, etc occurrences in SummarizedExperiment object
#' containing modificaton probabilities for SMF data
#'
#' @param se A \code{\link[SummarizedExperiment]{SummarizedExperiment}} object
#'     containing read-level data returned by the 
#'     \code{\link[footprintR]{readModBam}} function and containing 
#'     modification probabilities.
#' @param assayName Character scalar describing the name of the assay in 
#'     \code{se} containing read-level data.
#' @param threshUnmod,threshMod Numeric scalars used to classify observations
#'     as modified (modification probability >= threshMod), unmodified
#'     (modification probability < threshUnmod) or unknown (otherwise).
#' @param min_frag_data_len \code{integer} Ignore fragments that have genomic 
#'     lengths from most-left to most-right data points less than 
#'     \code{min_frag_data_len}.
#' @param min_frag_data_dens \code{numeric} Ignore fragments that have density 
#'     of data-containing positions lower than \code{min_frag_data_dens}.
#' @param max_spacing Maximum spacing between positions.
#' @param aggrSamples \code{logical} Return frequencies for each sample 
#'     separately or aggregate them for all samples.
#' @param ncpu Number of cores to use.
#' @param verbose Verbose mode for bug fixing.
#'
#' @return If \code{aggrSamples = TRUE}, a \code{matrix} containing aggregated 
#'     frequencies across all samples for 0,0; 0,1; etc for each spacing 
#'     between 1 (total number of 0s and 1s in the data) and \code{max_spacing}.
#'     If \code{aggrSamples = FALSE}, a \code{list} of matrices with frequencies 
#'     for each sample separately.
#'
#' @export
#'
get_ctable_from_SE <- function(se,
                               assayName = "mod_prob",
                               threshUnmod = 0.5,
                               threshMod = 0.5,
                               min_frag_data_len = 50L,
                               min_frag_data_dens = 0.05,
                               max_spacing = 200L,
                               aggrSamples = FALSE,
                               ncpu = 1L,
                               verbose = FALSE) {
    
    ### validate se object and prepare data for nomeR prediction
    dataList <- validate_prepare_SE(se = se,
                                    assayName = assayName,
                                    threshMod = threshUnmod,
                                    threshUnmod = threshMod,
                                    min_frag_data_len = min_frag_data_len,
                                    min_frag_data_dens = min_frag_data_dens)
    protect_data <- dataList[["bin_protect_data"]]
    
    ## collect pair stats for each sample
    ctable_list <- lapply(
        seq_len(ncol(se)),
        function(csidx) {
            count_spacing_freq_cpp(
                protect_data[sidx == csidx][["fidx_glob"]], ## unique fragment ID or index
                protect_data[sidx == csidx][["fragpos"]],   ## position within fragment, 1 - based
                protect_data[sidx == csidx][["protect"]],   ## binary protection data, 0 - accessible, 1 - protected
                max_spacing,
                ncpu,
                verbose)
        })
    names(ctable_list) <- colnames(se)
    
    if (aggrSamples) {
        return(Reduce("+", ctable_list))
    }
    return(ctable_list)
}
