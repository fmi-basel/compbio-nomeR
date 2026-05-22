#' Create expected co-occurrence count tables (N00, N01, N10, N11) from a
#' SummarizedExperiment object containing SMF modification probabilities
#'
#' @description
#' Computes expected co-occurrence counts for accessibility state pairs
#' (\code{N00}, \code{N01}, \code{N10}, \code{N11}) at varying spacings within a
#' \code{SummarizedExperiment} object containing single-molecule
#' footprinting (SMF) modification probabilities. Counts are calculated
#' for distances up to \code{max_spacing}, either per sample or aggregated
#' across all samples. Expected counts are computed as sums of products of
#' modification probabilities: e.g. \code{N00 += p_i * p_j} for each valid
#' position pair at spacing \code{s}, where high \code{p} = accessible.
#'
#' @inheritParams ftp_spectral_analysis_SE
#' @inheritParams predict_footprints_SE
#'
#' @param max_spacing Integer specifying the maximum spacing (in bases)
#'     between positions when computing expected co-occurrence counts for state
#'     pairs (\code{N00}, \code{N01}, \code{N10}, \code{N11}) across the SMF dataset.
#' @param aggrSamples Logical. If \code{TRUE}, expected co-occurrence counts
#'     are aggregated across all samples. If \code{FALSE}, a separate count
#'     matrix is returned for each sample.
#' @param ncpu Number of CPU cores to use for parallel processing.
#' @param verbose Logical. If \code{TRUE}, prints additional progress messages
#'     for debugging.
#'
#' @return
#' If \code{aggrSamples = TRUE}, a \code{matrix} containing aggregated
#' expected co-occurrence counts for \code{N00}, \code{N01}, \code{N10},
#' and \code{N11} at spacings from 1 to \code{max_spacing}.
#'
#' If \code{aggrSamples = FALSE}, a \code{list} of matrices, one per sample,
#' where each matrix contains the corresponding expected co-occurrence counts.
#'
#' @export

get_ctable_from_SE <- function(se,
                               assayName = "mod_prob",
                               min_frag_data_len = 50L,
                               min_frag_data_dens = 0.05,
                               max_spacing = 200L,
                               aggrSamples = FALSE,
                               ncpu = 1L,
                               verbose = FALSE) {

    ### validate se object and prepare data for footBayes prediction
    dataList <- validate_prepare_SE(se = se,
                                    assayName = assayName,
                                    min_frag_data_len = min_frag_data_len,
                                    min_frag_data_dens = min_frag_data_dens)
    mod_prob_data <- dataList[["mod_prob_data"]]

    ## collect pair stats for each sample
    ctable_list <- lapply(
        seq_len(ncol(se)),
        function(csidx) {

            count_spacing_freq_cpp(
                mod_prob_data[sidx == csidx][["fidx_glob"]], ## unique fragment ID or index
                mod_prob_data[sidx == csidx][["fragpos"]],   ## position within fragment, 1 - based
                mod_prob_data[sidx == csidx][["mod_prob"]], ## modification probability in [0,1]
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
