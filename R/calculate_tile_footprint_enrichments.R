#' Enrichment of Genomic Footprints Across Genomic Tiles
#'
#' @description
#' Calculates enrichment and statistical significance for specified footprint
#' types across genomic tiles. Unlike \code{\link{calculate_tile_BG_TF_enrichments}},
#' which uses log-ratio (ilr) scores that mix footprint types, this function
#' computes an independent log-score for each footprint type:
#' \code{log(coverage + psc)}. This avoids the correlation artefact between
#' Z-statistics that arises from log-ratio scoring when footprint prior
#' abundance is low.
#'
#' @param cover_dt A \code{data.table} containing posterior coverages returned by
#'   \code{\link{predict_footprints_SE}} with \code{returnAs="data.table"} (the
#'   \code{COVER_PROB} element).
#' @param tile_width Width of the sliding windows (tiles).
#' @param tile_step Step size for sliding windows.
#' @param ftp_names Character vector of footprint names (column names in
#'   \code{cover_dt}) for which to calculate enrichments (e.g.
#'   \code{c("background", "TF")}).
#' @param psc Pseudo-count added to posterior coverages before log-transformation.
#' @param threads Number of threads used by \code{data.table}.
#'   \code{NULL} (default) re-reads environment variables. \code{0} uses all
#'   available logical CPUs. Otherwise a positive integer.
#'
#' @returns A \code{data.table} with enrichment values and statistical
#'   significance for each footprint in \code{ftp_names}, calculated using
#'   a one-sided Z-test. The null hypothesis is that the mean log-score for a
#'   given footprint within a tile equals the genome-wide mean.
#' \describe{
#'   \item{seqnames, start, end, tile_ID}{Genomic coordinates of the sliding window (tile).}
#'   \item{n_data_points}{Number of data points (position-molecule pairs with
#'     modification data) aggregated across all molecules overlapping the tile.}
#'   \item{n_inf_pos}{Number of unique informative positions overlapping the tile.
#'     Read coverage is \code{n_data_points / n_inf_pos}.}
#'   \item{\{ftp\}_score_mean}{Mean log-score \code{mean(log(coverage + psc))}
#'     for footprint \code{ftp} across all molecules and informative positions
#'     overlapping the tile.}
#'   \item{\{ftp\}_score_mean_Zstat}{Z-statistic from a one-sided Z-test
#'     comparing the tile mean log-score to the genome-wide mean.}
#'   \item{\{ftp\}_score_mean_Ztest_pval}{P-value from the Z-test.}
#'   \item{\{ftp\}_score_mean_FDR}{FDR-adjusted p-value computed using
#'     \code{\link[stats]{p.adjust}} with \code{method="fdr"}.}
#' }
#'
#' @importFrom GenomicRanges GRanges slidingWindows
#' @importFrom IRanges IRanges
#' @importFrom stats p.adjust pnorm
#' @importFrom checkmate assertDataTable assert_character
#' @import data.table
#' @export
calculate_tile_footprint_enrichments <- function(cover_dt,
                                                  tile_width = 500,
                                                  tile_step = 250,
                                                  ftp_names,
                                                  psc = 0.1,
                                                  threads = NULL) {

    mod_prob <- i.start <- tile_ID <- n_data_points <- width <- strand <- NULL

    setDTthreads(threads = threads)

    assertDataTable(x = cover_dt, min.cols = 4)
    assert_character(ftp_names, min.len = 1, any.missing = FALSE)
    required_cols <- c("seqnames", "start", "mod_prob", ftp_names)
    if (!all(required_cols %in% colnames(cover_dt))) {
        stop("cover_dt must contain columns ",
             paste0(required_cols, collapse = ", "))
    }

    ## keep only informative positions
    cover_dt <- cover_dt[!is.na(mod_prob)]

    ## compute log-score for each footprint independently
    score_cols <- paste0(ftp_names, "_score")
    for (i in seq_along(ftp_names)) {
        cover_dt[, (score_cols[i]) := log(get(ftp_names[i]) + psc)]
    }


    ## global mean and population SD per footprint score (for Z-test denominator)
    N_total <- nrow(cover_dt)
    global_stats <- lapply(score_cols, function(sc) {
        vals <- cover_dt[[sc]]
        list(mean_total = mean(vals, na.rm = TRUE),
             sd_total   = sd(vals) * sqrt((N_total - 1L) / N_total))
    })
    names(global_stats) <- score_cols

    ## create sliding windows spanning observed positions per chromosome
    span_reg <- range(GRanges(seqnames = cover_dt[["seqnames"]],
                              IRanges(start = cover_dt[["start"]], width = 1L)))
    tiles_loci <- unlist(slidingWindows(span_reg,
                                        width = tile_width,
                                        step  = tile_step))
    tiles_loci$tile_ID <- as.character(tiles_loci)

    ## overlap cover_dt positions with tiles
    cover_dt <- cover_dt[, end := start]
    setkeyv(cover_dt, cols = c("seqnames", "start", "end"))
    tiles_dt <- data.table(as.data.frame(tiles_loci),
                           key = c("seqnames", "start", "end"))
    tiles_dt[, c("width", "strand") := list(NULL, NULL)]

    smftile_ov <- foverlaps(x = cover_dt,
                            y = tiles_dt,
                            which   = FALSE,
                            nomatch = NULL)

    ## aggregate per tile: total data points, unique positions, mean log-scores
    mean_cols <- paste0(score_cols, "_mean")
    tile_aggr_stats <- smftile_ov[, {
        res <- lapply(.SD, mean, na.rm = TRUE)
        names(res) <- mean_cols
        c(list(n_data_points = .N, n_inf_pos = uniqueN(i.start)), res)
    }, by = .(seqnames, start, end, tile_ID),
       .SDcols = score_cols]

    ## one-sided Z-test per footprint: tile mean vs genome-wide mean
    ## All arithmetic is done in plain R to avoid data.table scope ambiguity with
    ## dynamic column names; results are assigned back via (col) := vector
    nd <- tile_aggr_stats[["n_data_points"]]
    den <- as.numeric(N_total - nd)
    ratios <- (as.numeric(N_total) - 1) / as.numeric(den)
    size_factor <- ifelse(den > 0L, sqrt(nd * ratios), NA_real_)

    for (i in seq_along(score_cols)) {
        gs    <- global_stats[[score_cols[i]]]
        zsc   <- paste0(score_cols[i], "_mean_Zstat")
        zpc   <- paste0(score_cols[i], "_mean_Ztest_pval")
        fdc   <- paste0(score_cols[i], "_mean_FDR")

        zstat <- if (gs$sd_total > 0) {
            (tile_aggr_stats[[mean_cols[i]]] - gs$mean_total) /
                gs$sd_total * size_factor
        } else {
            rep(NA_real_, nrow(tile_aggr_stats))
        }

        tile_aggr_stats[, (zsc) := zstat]
        tile_aggr_stats[, (zpc) := pnorm(zstat, lower.tail = FALSE)]
        tile_aggr_stats[, (fdc) := p.adjust(tile_aggr_stats[[zpc]], method = "fdr")]
    }

    return(tile_aggr_stats)
}
