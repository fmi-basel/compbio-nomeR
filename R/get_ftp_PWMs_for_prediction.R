#' Create footprint models for \code{\link{predict_footprints}} using footprints
#' spectrum inferred by \code{\link{infer_footprints_vb}} and summarized by
#' \code{\link{get_ftp_inference_summary}}
#'
#' This utility function creates a list of parameters required by
#' \code{\link{predict_footprints}} to predict footprint positions in data.
#'
#' @param ftp_spectrum `data.frame` with inferred abundances of footprints.
#' @param ftp_len_mat a \code{matrix} with footprint lengths where 1st and 2nd columns
#' represent minimum and maximum footprint lengths. The 3rd column will be interpreted as
#' increment. Each row will result in PWMs with lengths seq(min_ftp_length,max_ftp_length,by).
#' Row names are interpreted as footprint groups.
#' @param bg_cover estimated percentage of accessible positions in SMF data
#' @param ftp_protect_prob emission probability of protected position within footprints
#'
#' @returns \code{list} with footprint models (PWM) required for the function \code{\link{predict_footprints}}.
#' @importFrom magrittr %>%
#' @importFrom dplyr filter mutate select
#' @importFrom ggplot2 ggplot aes geom_line labs theme theme_bw scale_y_continuous scale_x_continuous sec_axis
#' @importFrom checkmate assertDataFrame assertSubset assertMatrix
#'
#' @export
#'
# @examples

get_ftp_PWMs_for_prediction <- function(ftp_spectrum,
																				ftp_len_mat,
																				bg_cover,
																				ftp_protect_prob) {

	assertSubset(x = c("ftp_length","mean"),
							 choices = colnames(ftp_spectrum))
	ftp_spectrum <- ftp_spectrum %>% select(c("ftp_length","mean"))
	assertDataFrame(ftp_spectrum,
									all.missing = F)

	if(any(duplicated(ftp_spectrum$ftp_length)))
		warning("Found duplicated ftp_length. Please make sure that the input ftp_spectrum contain only one spectrum")

	assertMatrix(ftp_len_mat,any.missing = F,
							 min.cols = 2)

	## ftp_len_mat: 1 column - min_ftp_length, 2nd column - max_ftp_length,
	## if there is a third column - by, i.e increment from min to max

	stopifnot(all(ftp_len_mat[,2] >= ftp_len_mat[,1]))
	if(ncol(ftp_len_mat) == 2){
		if(is.null(colnames(ftp_len_mat)))
			colnames(ftp_len_mat) <- c("min_ftp_length","max_ftp_length")
		ftp_len_mat <- cbind(ftp_len_mat,
												 "by" = rep(1,nrow(ftp_len_mat)))
	} else if(ncol(ftp_len_mat) == 3){
		if(is.null(colnames(ftp_len_mat)))
			colnames(ftp_len_mat) <- c("min_ftp_length","max_ftp_length","by")
	}
	if(is.null(row.names(ftp_len_mat)))
		row.names(ftp_len_mat) <- paste0("ftp",ftp_len_mat[,1],"_",ftp_len_mat[,2])


	ftp_lengths <- sapply(row.names(ftp_len_mat),
												function(ridx) {
													seq(from = ftp_len_mat[ridx, 1],
															to = ftp_len_mat[ridx, 2],
															by = ftp_len_mat[ridx,3])
												}, simplify = FALSE, USE.NAMES = TRUE)
	ftp_cov <- sapply(row.names(ftp_len_mat),
										function(ridx) {
											ftp_spectrum <- ftp_spectrum %>%
												filter(ftp_length >= ftp_len_mat[ridx, 1] & ftp_length <= ftp_len_mat[ridx, 1])
											sum(ftp_spectrum$mean)
										}, simplify = T, USE.NAMES = TRUE)
	ftp_cov <- ftp_cov/sum(ftp_cov) * (1 - bg_cover)

# browser()
# 	## select required footprints and renormalize
# 	ftp_spectrum <- ftp_spectrum %>% filter(ftp_length %in% unlist(ftp_lengths)) %>%
# 		mutate(mean = mean/sum(mean) * (1 - bg_cover))

	## create models
	ftp_models <- do.call(
		c,
		lapply(
			names(ftp_lengths),
			function(nm) {
				cov_prior <- ftp_cov[nm]/length(ftp_lengths[[nm]])
				lapply(
					ftp_lengths[[nm]],
					function(flen) {


						list("PROTECT_PROB" = rep(ftp_protect_prob,flen),
								 "COVER_PRIOR" =cov_prior,
								 "NAME" = paste0(nm, "--", flen),
								 "GROUP" = nm)
					})
			}))
	return(ftp_models)


}
