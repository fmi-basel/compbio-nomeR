#' Create count tables of `00`, `01`, `10`, and `11` co-occurrences in a
#' SummarizedExperiment object containing SMF modification probabilities
#'
#' @description
#' Computes co-occurrence frequencies of binary accessibility states
#' (`00`, `01`, `10`, `11`) at varying spacings within a
#' \code{SummarizedExperiment} object containing single-molecule
#' footprinting (SMF) modification probabilities. Frequencies are calculated
#' for distances up to \code{max_spacing}, either per sample or aggregated
#' across all samples.
#'
#' @inheritParams ftp_spectral_analysis_SE
#' @inheritParams predict_footprints_SE
#'
#' @param max_spacing Integer specifying the maximum spacing (in bases)
#'     between positions when counting co-occurrences of state pairs
#'     (`00`, `01`, `10`, `11`) across the SMF dataset.
#' @param aggrSamples \code{logical}. If \code{TRUE}, co-occurrence frequencies
#'     are aggregated across all samples. If \code{FALSE}, a separate frequency
#'     matrix is returned for each sample.
#' @param ncpu Number of CPU cores to use for parallel processing.
#' @param verbose Logical. If \code{TRUE}, prints additional progress messages
#'     for debugging.
#'
#' @return
#' If \code{aggrSamples = TRUE}, a \code{matrix} containing aggregated
#' co-occurrence frequencies for `00`, `01`, `10`, and `11` at spacings
#' from 1 to \code{max_spacing}.
#'
#' If \code{aggrSamples = FALSE}, a \code{list} of matrices, one per sample,
#' where each matrix contains the corresponding co-occurrence frequencies.
#'
#' @export

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
