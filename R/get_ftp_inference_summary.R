#' Extract and/or plot estimates from footprint inference and optionally suggest footprints
#'
#' @description
#' Utility function to extract footprint abundance estimates from a
#' \code{stanfit} or optimization result (from \code{infer_footprints_vb},
#' \code{infer_footprints_sampling}, or \code{infer_footprints_optim}) and,
#' optionally, detect potential footprints from the abundance spectrum using
#' \code{\link{suggest_footprints}}.
#'
#' @param infer_stanfit A \code{\link[rstan]{stanfit}} object returned by
#'   \code{\link{infer_footprints_vb}} or \code{\link{infer_footprints_sampling}},
#'   or a \code{list} returned by \code{\link{infer_footprints_optim}}.
#' @param ftp_abundance_name Character scalar specifying which footprint
#'   abundance value to extract or display. Currently only \code{"ftp_abundances"} is supported.
#' @param plot Logical. If \code{TRUE}, generate a plot of the footprint abundance spectrum.
#' @param show_plot Logical. If \code{TRUE}, display the plot of the footprint spectrum.
#' @param suggest_ftps Logical. If \code{TRUE}, return suggested footprints based on spectrum peaks.
#' @param plot_posterior_range Character vector of length 2 specifying which posterior credible intervals to plot.
#'   Options provided by \code{\link[rstan]{summary,stanfit-method}} include "2.5\%", "25\%", "75\%", "97.5\%".
#'   Standard deviation (SD) and standard error (SE) are not currently supported.
#' @param spline_spar Numeric in (0, 1]. Smoothing parameter for
#'   \code{\link{suggest_footprints}} controlling the smoothness of the
#'   spline (\code{\link[stats]{smooth.spline}}). Recommended to test values
#'   like 0.1, 0.3, 0.5, 0.75 to ensure reasonable footprint detection.
#' @param max_abund_log2drop Numeric. Maximum decrease in log2 abundance relative
#'   to the local maximum when extending peaks in \code{\link{suggest_footprints}}.
#' @param max_peak_width Numeric. Maximum allowed width of detected peaks in
#'   the footprint abundance spectrum (used by \code{\link{suggest_footprints}}).
#' @param ... Additional parameters passed to \code{\link{suggest_footprints}},
#'   including those for \code{\link[stats]{smooth.spline}}.
#'
#' @return A \code{list} containing:
#' \describe{
#'   \item{\code{ESTIMATES}}{A \code{list} containing:
#'     \itemize{
#'       \item \code{ftp_abundance_estimates}: \code{data.frame} of footprint abundance estimates,
#'       \item \code{ftp_protect_prob_estimate}: estimated footprint protection probability (or \code{NA} if unavailable),
#'       \item \code{bg_protect_prob_estimate}: estimated background protection probability (or \code{NA} if unavailable).
#'     }}
#'
#'   \item{\code{FTP_SUGGEST}}{\code{matrix} with coordinates of suggested footprints returned by \code{\link{suggest_footprints}}.}
#'   \item{\code{PLOT}}{\code{ggplot} object of the footprint abundance spectrum.}
#' }
#'
#' @export
#'
#' @examples
#'
#' ## simple data with two ftp of 5 and 10 bps.
#' ## The table below is a count table of observed occurrences of
#' ## 00, 01, 10 and 11 at spacings from 1 until 15.
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
#' ## variational inference for footprints in the data
#' inf <- infer_footprints_vb(cooc_ctable = ftp_5_10_data,
#'                            ftp_lengths = 2:15)
#'
#' ## get estimates and plot footprint spectrum
#' inference_summary_list <- get_ftp_inference_summary(inf, plot = TRUE)
#'
#' @importFrom ggplot2 ggplot geom_ribbon aes labs geom_line scale_y_log10
#'     scale_color_manual scale_fill_manual guides guide_legend theme_bw theme
#'     element_text annotate geom_point scale_x_continuous
#' @importFrom graphics plot
#' @importFrom stringr str_extract
#' @importFrom rlang .data
#' @importFrom rstan summary
get_ftp_inference_summary <- function(
        infer_stanfit,
        plot = FALSE,
        show_plot = plot,
        suggest_ftps = FALSE,
        plot_posterior_range = c("2.5%", "97.5%"),
        max_peak_width = 30,
        spline_spar = 0.6,
        max_abund_log2drop = 1.5,
        ftp_abundance_name = c("ftp_abundances"),
        ...) {

    ftp_abundance_name <- match.arg(ftp_abundance_name)
    ftp_abund_pattern <- paste0(ftp_abundance_name, "\\[\\d+\\]")

    ## get ftp_lengths
    if (!is.null(attr(infer_stanfit, "ftp_lengths"))) {
        ftp_lengths <- attr(infer_stanfit, "ftp_lengths")
        ## add background length
        ftp_lengths <- c(1, ftp_lengths)
    } else {
        stop("infer_stanfit must have attribute ftp_lengths")
    }

    ## get summary from stanfit
    if (inherits(infer_stanfit, "stanfit")) {
        ftp_infer_summary <- as.data.frame(
            rstan::summary(infer_stanfit)$summary)
        ftp_infer_summary$param <- gsub("\\[\\d+\\]$", "",
                                        row.names(ftp_infer_summary))

        ## get vector of ftp coverages
        infer_ftp_abund_probs <-
            ftp_infer_summary[grep(ftp_abund_pattern,
                                   row.names(ftp_infer_summary),
                                   perl = TRUE), ]

        infer_ftp_abund_probs$ftp_length <-
            ftp_lengths[as.numeric(gsub("[\\[\\]]", "",
                                        str_extract(
                                            row.names(infer_ftp_abund_probs),
                                            pattern = "\\[(\\d+)\\]"),
                                        perl = TRUE))]

        row.names(infer_ftp_abund_probs) <- infer_ftp_abund_probs$ftp_length
        # estimate for ftp_protect_prob
        infer_ftp_protect_prob <- ftp_infer_summary["ftp_protect_prob", ]

        ## if model with informative prior for bg_protect_prob get estimate,
        ## otherwise NA
        if (infer_stanfit@model_name == "ftp_inference_informative_prior") {
            infer_bg_protect_prob <- ftp_infer_summary["bg_protect_prob", ]
        } else if (infer_stanfit@model_name == "ftp_inference_bg_fixed") {
            infer_bg_protect_prob <- NA
        } else {
            warning("Model name: ", infer_stanfit@model_name,
                    ". Setting infer_bg_protect_prob, infer_ftp_protect_prob ",
                    "to NA")
            infer_bg_protect_prob <- NA
            infer_ftp_protect_prob <- NA
        }
    } else if (inherits(infer_stanfit, "list")) {
        opt_ftp_cover <- infer_stanfit$par[grep(ftp_abund_pattern,
                                                names(infer_stanfit$par),
                                                perl = TRUE)]
        infer_ftp_abund_probs <- data.frame(mean = opt_ftp_cover)
        infer_ftp_abund_probs$param <- gsub("\\[\\d+\\]$", "",
                                            names(opt_ftp_cover))

        infer_ftp_abund_probs$ftp_length <- as.numeric(
            gsub("[\\[\\]]", "",
                 str_extract(names(opt_ftp_cover),
                             pattern = "\\[(\\d+)\\]"), perl = TRUE))

        tmp_na_mat <- matrix(NA, ncol = 7, nrow = nrow(infer_ftp_abund_probs))
        colnames(tmp_na_mat) <- c("se_mean", "sd", "2.5%", "25%", "50%",
                                  "75%", "97.5%")
        infer_ftp_abund_probs <- cbind(infer_ftp_abund_probs, tmp_na_mat)

        if ("bg_protect_prob" %in% names(infer_stanfit$par)) {
            infer_bg_protect_prob <- infer_stanfit$par["bg_protect_prob"]
        } else {
            infer_bg_protect_prob <- NA
        }

        if ("ftp_protect_prob" %in% names(infer_stanfit$par)) {
            infer_ftp_protect_prob <- infer_stanfit$par["ftp_protect_prob"]
        } else {
            infer_ftp_protect_prob <- NA
        }
    } else {
        stop("infer_stanfit must be stanfit object returned by ",
             "rstan::sampling or rstan::vb, or list returned by ",
             "rstan::optimizing")
    }

    if (plot) {
        infer_plot <- plot_ftp_spectrum(
            infer_ftp_abund_probs[infer_ftp_abund_probs$ftp_length != 1, ])
    } else {
        infer_plot <- NULL
    }

    if (suggest_ftps) {
        ## get ftp suggestions
        infer_ftp_abund_probs_subset <-
            infer_ftp_abund_probs[infer_ftp_abund_probs$ftp_length != 1, ]
        tryCatch(ftp_suggestions <- suggest_footprints(
            S = infer_ftp_abund_probs_subset$ftp_length,
            y = log2(infer_ftp_abund_probs_subset$mean),
            isLog = TRUE,
            max_peak_width = max_peak_width,
            spline_spar = spline_spar,
            max_abund_log2drop = max_abund_log2drop,
            ...),
            error = function(e) {
                print(e)
                ftp_ranges <- matrix(nrow = 0, ncol = 2)
                colnames(ftp_ranges) <- c("min_ftp_length", "max_ftp_length")
                ftp_suggestions <- list("ftp_ranges" = ftp_ranges,
                                        "smoothed_signal" = NULL)
            })

        if (nrow(ftp_suggestions[["ftp_ranges"]]) > 0) {
            ftp_lengths_suggest <- unlist(apply(
                ftp_suggestions[["ftp_ranges"]], 1,
                function(ftp_rng) {
                    seq(ftp_rng[1], ftp_rng[2])
                }))

            if (plot) {
                infer_plot <- infer_plot +
                    annotate(geom = "rect",
                             xmin = ftp_suggestions$ftp_ranges[, 1],
                             xmax = ftp_suggestions$ftp_ranges[, 2],
                             ymin = 0, ymax = Inf, alpha = 0.2,
                             color = "NA", fill = "grey") +
                    geom_point(
                        data = infer_ftp_abund_probs[
                            infer_ftp_abund_probs$ftp_length %in%
                                ftp_lengths_suggest, ],
                        mapping = aes(x = .data$ftp_length,
                                      y = .data$mean,
                                      color = .data$param),
                        color = "red") +
                    annotate(geom = "text",
                             x = rowMeans(ftp_suggestions$ftp_ranges),
                             y = 0,
                             hjust = 0.5, vjust = 0,
                             label = row.names(ftp_suggestions$ftp_ranges))
                ## add x breaks for suggested footprints
                x_breaks <- pretty(range(
                    infer_ftp_abund_probs[infer_ftp_abund_probs$ftp_length != 1,
                                          "ftp_length"], na.rm = TRUE),
                    n = 5)
                x_breaks <- sort(c(x_breaks,
                                   as.vector(ftp_suggestions$ftp_ranges)))
                infer_plot <- infer_plot +
                    scale_x_continuous(breaks = x_breaks)
            }
        }
    } else {
        ftp_suggestions <- NULL
    }

    if (plot && show_plot) {
        plot(infer_plot)
    }

    return(invisible(list(
        "ESTIMATES" = list("ftp_abundance_estimates" = infer_ftp_abund_probs,
                           "ftp_protect_prob_estimate" = infer_ftp_protect_prob,
                           "bg_protect_prob_estimate" = infer_bg_protect_prob),
        "FTP_SUGGEST" = ftp_suggestions[["ftp_ranges"]],
        "PLOT" = infer_plot)
    ))
}
