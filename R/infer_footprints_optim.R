#' Find point estimate for footprint abundance using STAN optimization
#' algorithm
#'
#' This function uses `rstan::optimizing` to obtain point estimate for
#' footprint abundances by maximizing the joint posterior from the model.
#'
#' @inheritParams infer_footprints_vb
#' @param ... parameters passed to \code{\link[rstan]{optimizing}} function from the \code{rstan} package. Please
#'     refer to \code{\link[rstan]{optimizing}} documentation.
#'
#' @return A list with components described in \code{\link[rstan]{optimizing}}
#' function.
#' The attribute \code{attr(<list>,"ftp_lengths")} contains a vector
#' of footprint lengths for which inference was run, i.e. parameter
#' \code{ftp_lengths} provided by user.
#'
#' @export
#'
#' @examples
#'
#' ## Simple data with two footprints of lengths 5 and 10 bps.
#' ## The table below is a count table of observed occurrences of
#' ## 00, 01, 10, and 11 at spacings from 1 until 15.
#'
#' ftp_5_10_data <- data.frame("S" = 1:15,
#'                             "N00" = c(1626964, 1508381,
#'                                       1420066, 1336897,
#'                                       1258679, 1185045,
#'                                       1136784, 1090029,
#'                                       1045066, 1001670,
#'                                       959927, 983020,
#'                                       1000410, 1012326,
#'                                       1019633),
#'                             "N01" = c(0, 113856,
#'                                       197657, 276539,
#'                                       350664, 420399,
#'                                       464873, 507988,
#'                                       549466, 589487,
#'                                       627997, 601629,
#'                                       580965, 565776,
#'                                       555217),
#'                             "N10" = c(0, 113921,
#'                                       197854, 276862,
#'                                       351147, 421042,
#'                                       465716, 509005,
#'                                       550645, 590837,
#'                                       629533, 603328,
#'                                       582814, 567737,
#'                                       557220),
#'                             "N11" = c(873036, 758842,
#'                                       674423, 594702,
#'                                       519510, 448514,
#'                                       402627, 357978,
#'                                       314823, 273006,
#'                                       232543, 257023,
#'                                       275811, 289161,
#'                                       297930))
#'
#'
#' ## finding MAP estimate
#' inf_output <- infer_footprints_optim(cooc_ctable = ftp_5_10_data,
#'                                      ftp_lengths = 2:15)
#'
#' ## plot footprint spectrum
#' get_ftp_inference_summary(inf_output, plot = TRUE)
#'
#' @importFrom rstan optimizing
infer_footprints_optim <- function(
        cooc_ctable,
        ftp_lengths,
        ftp_prior_cover = NULL,
        bg_prior_cover = 0.5,
        total_cnt_prior_dirich = NULL,
        ftp_bg_model = c("informative_prior", "bg_fixed", "ftp_bg_fixed"),
        bg_model_params = list("bg_protect_prob_fixed" = 0.05,
                               "bg_protect_min" = 0.01,
                               "bg_protect_max" = 0.2,
                               "bg_protect_mean" = 0.05,
                               "bg_protect_totcount" = 100),
        ftp_model_params = list("ftp_protect_prob_fixed" = 0.95,
                                "ftp_protect_min" = 0.8,
                                "ftp_protect_max" = 0.99,
                                "ftp_protect_mean" = 0.95,
                                "ftp_protect_totcount" = 100),
        ...) {

    ftp_bg_model <- match.arg(ftp_bg_model)

    ## validate and construct input for inference
    stan_input <- .validate_construct_stan_input(cooc_ctable,
                                                 ftp_lengths,
                                                 bg_prior_cover,
                                                 ftp_prior_cover,
                                                 total_cnt_prior_dirich,
                                                 ftp_bg_model,
                                                 bg_model_params,
                                                 ftp_model_params)

    ## get initial values for fitting
    stan_initvals <- .init_param_from_prior_distr(stan_input = stan_input,
                                                  nchains = 1)
    ## Stan optimizing
    stanfit_out <- rstan::optimizing(
        object = stanmodels[[stan_input$stan_model_name]],
        data = stan_input$stan_inputdata,
        init = stan_initvals[[1]],
        ...)

    attr(stanfit_out,"ftp_lengths") <- ftp_lengths

    return(stanfit_out)
}
