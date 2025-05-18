#' Extract and/or plot estimates from footprint inference and (optionally) 
#' identify potential footprints from spectrum
#'
#' This utility function extract estimates from \code{infer_stanfit} and uses 
#' function \code{\link{suggest_footprints}} (optionally)
#' to detect peaks in spectrum and return lengths of potential footprints.
#'
#' @param infer_stanfit \code{\link[rstan]{stanfit}} object returned by
#'     \code{\link{infer_footprints_vb}}, 
#'     \code{\link{infer_footprints_sampling}} or a \code{list} returned by
#'     \code{\link{infer_footprints_optim}}.
#' @param ftp_abundance_name \code{character} which defines which footprint 
#'     abundance value to return/display. Currently only ftp_abundances is 
#'     accepted.
#' @param plot \code{logical} return plot footprint abundance spectrum.
#' @param show_plot \code{logical} show footprint abundance spectrum.
#' @param suggest_ftps \code{logical} return suggestion for footprints.
#' @param plot_posterior_range \code{vector} of \code{character} of length 2 
#'     which specifies which posterior range to plot. 
#'     \code{\link[rstan]{summary,stanfit-method}} provides "2.5\%", "25\%" 
#'     and "75\%", "97.5\%" credible intervals as well as SD and SE. Displaying 
#'     SD and SE is not implemented at the moment.
#' @param spline_spar parameter for function \code{\link{suggest_footprints}} 
#'     controlling smoothness of spline ((0,1], the higher the smoother). Please 
#'     see \code{\link{suggest_footprints}} \code{\link[stats]{smooth.spline}}.
#'     It is recommended to test different values for \code{spline_spar}, for 
#'     example 0.1, 0.3, 0.5, 0.75 to check whether suggested footprints look 
#'     as expected.
#' @param max_abund_log2drop parameter for function 
#'     \code{\link{suggest_footprints}}, namely maximum decrease in abundance 
#'     relative to value at local maxima until which peaks are extended. Please 
#'     see \code{\link{suggest_footprints}}
#' @param max_peak_width parameter for function 
#'     \code{\link{suggest_footprints}}, namely maximum width of detected peaks 
#'     in spectrum. Please see \code{\link{suggest_footprints}}.
#' @param ... parameters for \code{\link{suggest_footprints}} which are 
#'     transmitted to function \code{\link[stats]{smooth.spline}}.
#'
#' @return \code{list} containing elements:
#' \code{"ESTIMATES"} - \code{list} which contains a \code{ftp_abundance_estimates}  - \code{data.frame} for footprint abundances estimates,
#' \code{"ftp_protect_prob_estimate"} and \code{bg_protect_prob_estimate} which are estimates for \code{ftp_protect_prob} and \code{bg_protect_prob}.
#' \code{"bg_protect_prob_estimate"} and/or \code{"ftp_protect_prob_estimate"} are \code{NA} if inference for these parameters are not available.
#' \code{"FTP_SUGGEST"} - \code{matrix} with coordinates for suggested footprints as returned by \code{\link{suggest_footprints}}.
#' \code{"PLOT"} - a \code{ggplot} object with plot for footprint spectrum.
#'
#' @export
#'
#' @examples
#' \dontrun{
#'
#' ## simple data with two ftp of 5 and 10 bps.
#' ## The table below is a count table of observed occurrences of
#' ## 00, 01, 10 and 11 at spacings from 1 until 15.
#'
#' ftp_5_10_data <- data.frame("S" = 1:15,
#'                             "N00" = c(1626964,1508381,
#'                                       1420066,1336897,
#'                                       1258679,1185045,
#'                                       1136784,1090029,
#'                                       1045066,1001670,
#'                                       959927,983020,
#'                                       1000410,1012326,
#'                                       1019633),
#'                             "N01" = c(0,113856,
#'                                       197657,276539,
#'                                       350664,420399,
#'                                       464873,507988,
#'                                       549466,589487,
#'                                       627997,601629,
#'                                       580965,565776,
#'                                       555217),
#'                             "N10" = c(0, 113921,
#'                                       197854,276862,
#'                                       351147,421042,
#'                                       465716,509005,
#'                                       550645,590837,
#'                                       629533,603328,
#'                                       582814,567737,
#'                                       557220),
#'                             "N11" = c(873036,758842,
#'                                       674423,594702,
#'                                       519510,448514,
#'                                       402627,357978,
#'                                       314823,273006,
#'                                       232543,257023,
#'                                       275811,289161,
#'                                       297930))
#'
#' ## variational inference for footprints in the data
#' inf <- infer_footprints_vb(joint_freq_table = ftp_5_10_data,
#'                                       footprint_prior_diralphas = c(10,rep(1,14)))
#'
#' ## get estimates and plot footprint spectrum
#' inference_summary_list <- get_ftp_inference_summary(inf,plot=T)
#' }
#'
#' @importFrom ggplot2 ggplot geom_ribbon aes labs geom_line scale_y_log10 
#'     scale_color_manual scale_fill_manual guides guide_legend theme_bw theme 
#'     element_text annotate geom_point scale_x_continuous
#' @importFrom graphics plot
#' @importFrom stringr str_extract
#' @importFrom rlang .data
#' @importFrom rstan summary
get_ftp_inference_summary <- function(infer_stanfit,
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
        ## create plot, do not show ftp_length = 1, i.e. background
        infer_plot <- ggplot(
            infer_ftp_abund_probs[infer_ftp_abund_probs$ftp_length != 1, ])
        
        if (inherits(infer_stanfit, "stanfit")) {
            infer_plot <- infer_plot +
                geom_ribbon(mapping = aes(
                    x = .data$ftp_length,
                    y = .data$mean,
                    ymin = .data[[plot_posterior_range[1]]],
                    ymax =  .data[[plot_posterior_range[2]]],
                    fill = .data$param),
                    alpha = 0.2, color = "grey") +
                labs(x = "footprint length, bp",
                     y = paste0("mean, ", plot_posterior_range[1], " - ",
                                plot_posterior_range[2], " interval"))
        } else if (inherits(infer_stanfit, "list")) {
            infer_plot <- infer_plot +
                labs(x = "footprint length, bp",
                     y = "MAP point estimate")
        }
        infer_plot <- infer_plot +
            geom_line(mapping = aes(x = .data$ftp_length,
                                    y = .data$mean,
                                    color = .data$param),
                      size = 1.5
            ) +
            scale_y_log10() +
            scale_color_manual(values = c("ftp_abundances" = "darkgreen")) +
            scale_fill_manual(values = c("ftp_abundances" = "lightgreen")) +
            guides(color = guide_legend(title = NULL),
                   fill = guide_legend(title = NULL)) +
            theme_bw() +
            theme(
                legend.position = "right",
                axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5,
                                           size = 7),
                axis.title = element_text(face = "bold", size = 12)
            )
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
