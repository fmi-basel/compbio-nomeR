#' Create footprint models for \code{\link{predict_footprints}} from an inferred footprint spectrum
#'
#' @description
#' Utility function to generate a list of footprint models required by
#' \code{\link{predict_footprints}}, using a footprint spectrum inferred by
#' \code{\link{infer_footprints_vb}} or \code{\link{ftp_spectral_analysis_SE}}
#' and summarized with \code{\link{get_ftp_inference_summary}}. The models encode footprint lengths,
#' emission probabilities, and background coverage for downstream prediction.
#'
#' @param ftp_spectrum A \code{data.frame} containing inferred abundances of footprints.
#' @param ftp_len_mat A \code{matrix} specifying footprint lengths:
#'   \describe{
#'     \item{Column 1}{Minimum footprint length.}
#'     \item{Column 2}{Maximum footprint length.}
#'     \item{Column 3}{Increment for generating lengths from minimum to maximum.}
#'   }
#'   Each row corresponds to a footprint group, with row names interpreted as group labels.
#'   PWMs will be generated for all lengths in \code{seq(min, max, by)}.
#' @param bg_cover Numeric value indicating the estimated fraction of accessible positions
#'   in the SMF dataset (used as background coverage).
#' @param ftp_protect_prob Numeric value of the emission probability for protected positions
#'   within footprints.
#'
#' @return A \code{list} of footprint models (position weight matrices) suitable for the
#'   \code{footprint_models} parameter in \code{\link{predict_footprints}} and \code{\link{predict_footprints_SE}}.
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr filter mutate select
#' @importFrom ggplot2 ggplot aes geom_line labs theme theme_bw
#'     scale_y_continuous scale_x_continuous sec_axis
#' @importFrom checkmate assertDataFrame assertSubset assertMatrix
#'
#' @export

# @examples

get_ftp_PWMs_for_prediction <- function(ftp_spectrum,
                                        ftp_len_mat,
                                        bg_cover,
                                        ftp_protect_prob) {

    assertSubset(x = c("ftp_length", "mean"),
                 choices = colnames(ftp_spectrum))
    ftp_spectrum <- ftp_spectrum %>% select(c("ftp_length", "mean"))
    assertDataFrame(ftp_spectrum,
                    all.missing = FALSE)

    if (any(duplicated(ftp_spectrum$ftp_length))) {
        warning("Found duplicated ftp_length. Please make sure that the ",
                "input ftp_spectrum contain only one spectrum")
    }

    assertMatrix(ftp_len_mat, any.missing = FALSE,
                 min.cols = 2)

    ## ftp_len_mat: 1 column - min_ftp_length, 2nd column - max_ftp_length,
    ## if there is a third column - by, i.e increment from min to max

    stopifnot(all(ftp_len_mat[, 2] >= ftp_len_mat[, 1]))
    if (ncol(ftp_len_mat) == 2) {
        if (is.null(colnames(ftp_len_mat))) {
            colnames(ftp_len_mat) <- c("min_ftp_length", "max_ftp_length")
        }
        ftp_len_mat <- cbind(ftp_len_mat,
                             "by" = rep(1, nrow(ftp_len_mat)))
    } else if (ncol(ftp_len_mat) == 3) {
        if (is.null(colnames(ftp_len_mat))) {
            colnames(ftp_len_mat) <- c("min_ftp_length", "max_ftp_length", "by")
        }
    }
    if (is.null(row.names(ftp_len_mat))) {
        row.names(ftp_len_mat) <- paste0("ftp", ftp_len_mat[, 1], "_",
                                         ftp_len_mat[, 2])
    }

    ftp_lengths <- sapply(row.names(ftp_len_mat),
                          function(ridx) {
                              seq(from = ftp_len_mat[ridx, 1],
                                  to = ftp_len_mat[ridx, 2],
                                  by = ftp_len_mat[ridx, 3])
                          }, simplify = FALSE, USE.NAMES = TRUE)
    ftp_cov <- sapply(row.names(ftp_len_mat),
                      function(ridx) {
                          selspc <- ftp_spectrum %>%
                              filter(ftp_length >= ftp_len_mat[ridx, 1] &
                                         ftp_length <= ftp_len_mat[ridx, 2])
                          sum(selspc$mean)
                      }, simplify = TRUE, USE.NAMES = TRUE)
    ftp_cov <- ftp_cov/sum(ftp_cov) * (1 - bg_cover)


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
                             "COVER_PRIOR" = cov_prior,
                             "NAME" = paste0(nm, "--", flen),
                             "GROUP" = nm)
                    })
            }))
    return(ftp_models)
}
