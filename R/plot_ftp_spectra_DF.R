#' Plot footprint spectra from a DataFrame
#'
#' @description
#' Generates a panel of footprint spectrum plots for all samples in a
#' \code{DataFrame} obtained from \code{\link{ftp_spectral_analysis_SE}}.
#' The input \code{DataFrame} must contain a column named \code{ftp_spectrum}.
#'
#' @param DF A \code{DataFrame} containing footprint spectra for each sample.
#'   Must include a column \code{ftp_spectrum} as returned by
#'   \code{\link{ftp_spectral_analysis_SE}}.
#'
#' @return A \code{patchwork} object containing a panel of footprint spectrum
#'   plots, one for each sample in \code{DF}.
#'
#' @importFrom patchwork wrap_plots
#'
#' @export
plot_ftp_spectra_DF <- function(DF) {
    stopifnot("ftp_spectrum" %in% colnames(DF))

    ## create list of plots
    pllist <- lapply(seq_len(nrow(DF)),
                     function(i) {
                         plot_ftp_spectrum(DF$ftp_spectrum[[i]],
                                           title = DF$sample[i])
                     })
    wrap_plots(pllist,
               axis_title = "collect")
}
