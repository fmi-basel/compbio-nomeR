#' Plot a footprint abundance spectrum
#'
#' @description
#' Generates a \code{ggplot2} plot of inferred footprint abundances across
#' different footprint lengths. The plot visualizes both the mean coverage
#' estimates and a measure of statistical significance of estimated values.
#'
#' @param ftp_spectrum A \code{data.frame} containing inferred footprint
#'   abundances. Must include the columns \code{ftp_length} and \code{mean}.
#' @param title Optional character string to use as the plot title.
#'
#' @return A \code{ggplot2} object.
#'   - The **dark green line** (left Y-axis, log10 scale) represents the
#'     estimated coverage (mean across samples from the posterior distribution)
#'     for each footprint length on the X-axis.
#'   - The **red line** (right Y-axis, linear scale) represents the ratio of
#'     mean to standard deviation across samples for each footprint length,
#'     which can be interpreted as approximate Z-scores. These values can be
#'     used to assess the statistical significance of the inferred footprint
#'     coverage estimates.
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr filter mutate select
#' @importFrom ggplot2 ggplot aes geom_line labs theme theme_bw
#'     scale_y_continuous scale_x_continuous sec_axis
#' @importFrom checkmate assertDataFrame assertSubset assertNumeric
#'     assertIntegerish
#'
#' @export
plot_ftp_spectrum <- function(ftp_spectrum,
                              title = NULL) {

    assertSubset(x = c("ftp_length", "mean", "sd"),
                 choices = colnames(ftp_spectrum))
    ftp_spectrum <- ftp_spectrum %>% select(c("ftp_length", "mean", "sd"))
    assertNumeric(ftp_spectrum[["mean"]], lower = 0, upper = 1,
                  all.missing = FALSE)
    assertIntegerish(ftp_spectrum[["ftp_length"]], lower = 1,
                     all.missing = FALSE)

    plotZscore <- !all(is.na(ftp_spectrum[["sd"]]))

    if (any(duplicated(ftp_spectrum$ftp_length))) {
        warning("Found duplicated ftp_length. Please make sure that the ",
                "input ftp_spectrum contains only one spectrum.")
    }

    ## calculate log10mean and zscore and scale them
    ftp_spectrum <- ftp_spectrum %>%
        mutate(log10mean = log10(mean),
               zscore = mean / sd)

    ## scale them
    lgmn_range <- range(ftp_spectrum$log10mean, na.rm = TRUE)
    ftp_spectrum <- ftp_spectrum %>%
        mutate(log10mean_scaled = (log10mean - lgmn_range[1]) /
                   diff(lgmn_range))
    if (plotZscore) {
        zsc_range <- range(ftp_spectrum$zscore, na.rm = TRUE)
        ftp_spectrum <- ftp_spectrum %>%
            mutate(zscore_scaled = (zscore - zsc_range[1]) / diff(zsc_range))
    }

    ## get axis breaks
    lgmn_breaks <- pretty(lgmn_range)
    lgmn_scl_breaks <- (lgmn_breaks - lgmn_range[1]) / diff(lgmn_range)
    mean_labs <- sapply(lgmn_breaks, function(x) parse(text = paste0("10^", x)))

    if (plotZscore) {
        zsc_breaks <- pretty(zsc_range)
        zsc_scl_breaks <- (zsc_breaks - zsc_range[1]) / diff(zsc_range)
    }

    xbreaks <- sort(c(pretty(ftp_spectrum$ftp_length,
                             bounds = FALSE),
                      range(ftp_spectrum$ftp_length, na.rm = TRUE)))

    ## create plot
    ftp_spec_pl <- ggplot(ftp_spectrum,aes(x = ftp_length)) +
        geom_line(aes(y = log10mean_scaled),
                  color = "darkgreen",
                  linewidth = 1.1) +
        labs(x = "footprint length, bp", title = title) +
        scale_x_continuous(breaks = xbreaks,
                           limits = range(ftp_spectrum$ftp_length)) +
        theme_bw() +
        theme(legend.position = "right",
              axis.text.x = element_text(angle = 0, hjust = 0.5,
                                         vjust = 0.5, size = 12, face = "bold"),
              axis.title.y.left = element_text(color = "darkgreen",
                                               face = "bold", size = 18),
              axis.text.y.left = element_text(color = "darkgreen", size = 8,
                                              face = "bold"),
              axis.title.y.right = element_text(color = "red2", face = "bold",
                                                size = 18),
              axis.text.y.right = element_text(color = "red2", size = 8,
                                               face = "bold"),
              axis.title.x = element_text(face = "bold", size = 18)

        )

    if (plotZscore) {
        ftp_spec_pl <- ftp_spec_pl +
            geom_line(aes(y = zscore_scaled),
                      color = "red2",
                      linewidth = 1.1,
                      alpha = 0.75) +
            scale_y_continuous(
                name = "estimated coverage (mean)",
                breaks = lgmn_scl_breaks,
                labels = mean_labs,
                sec.axis = sec_axis(
                    transform = ~ .* diff(zsc_range) + zsc_range[1],
                    name = "mean/sd",
                    breaks = zsc_breaks))
    } else {
        ftp_spec_pl <- ftp_spec_pl +
            scale_y_continuous(name = "estimated coverage (mean)",
                               breaks = lgmn_scl_breaks,
                               labels = mean_labs)
    }
    return(ftp_spec_pl)
}
