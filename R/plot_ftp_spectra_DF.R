#' Plot footprint spectra provided in the DataFrame
#'
#' @param DF \code{DataFrame} obtained using
#'       [ftp_spectral_analysis_SE()] and must contain column \code{ftp_spectrum}
#'
#' @returns \code{patchwork} object with panel of plots for footprint spectra for all samples in DF.
#' @importFrom patchwork wrap_plots
#' @export
#'
plot_ftp_spectra_DF <- function(DF){
	stopifnot("ftp_spectrum" %in% colnames(DF))

	## create list of plots
	pllist <- lapply(1:nrow(DF),
									 function(i){
									 	plot_ftp_spectrum(DF$ftp_spectrum[[i]],
									 										title = DF$sample[i])
									 })
	wrap_plots(pllist,
						 guides = "auto",
						 axes = "collect",
						 axis_title = "collect")

}
