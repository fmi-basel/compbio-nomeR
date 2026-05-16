#' Enrichment of Accessible Sites and TF Footprints Across Genomic Tiles
#'
#' @description
#' Calculates enrichments and statistical significance for accessible positions
#' (background) and transcription factor (TF) footprints across genomic tiles.
#'
#' @param cover_dt A \code{data.table} containing posterior coverages returned by
#'   \code{\link{predict_footprints_SE}} with \code{returnAs="data.table"}.
#' @param tile_width Width of the sliding windows (tiles).
#' @param tile_step Step size for sliding windows.
#' @param bg_colname Name of the column containing posterior coverage for background (accessible positions).
#' @param tf_colname Name of the column containing posterior coverage for TF footprints.
#' @param nucl_colname Name of the column containing posterior coverage for nucleosome footprints.
#' @param psc Pseudo-count added when calculating background and TF scores.
#' @param threads Number of threads used by \code{data.table}.
#'   \code{NULL} (default) re-reads environment variables. \code{0} uses all
#'   available logical CPUs. Otherwise a positive integer.
#'
#' @returns A \code{data.table} with enrichment values and statistical
#'   significance calculated using a Z-test on BG and TF scores.
#'   The null hypothesis is that the mean BG/TF score for each tile equals the
#'   mean across all tiles.
#' The table contains the following columns:
#' \describe{
#'   \item{seqnames, start, end, tileID}{Genomic coordinates of sliding windows (tiles).}
#'   \item{n_data_points}{Number of informative positions with modification data aggregated across all molecules.}
#'   \item{bg_score_mean, tf_score_mean}{Mean background and TF scores across all molecules and informative positions overlapping each tile.}
#'   \item{bg_pos_cnt, tf_pos_cnt}{Number of positions classified as "positive" for accessibility (background) or TF footprints,
#'         after applying \code{bg_score_thresh} and \code{tf_score_thresh}, aggregated across all molecules overlapping the tile.}
#'   \item{bg_log2enr, tf_log2enr}{Log2 ratios of observed over expected counts for accessible and TF-covered positions.}
#'   \item{bg_binom_pval, tf_binom_pval, bg_FDR, tf_FDR}{P-values (and FDR-adjusted p-values) from one-tailed binomial tests
#'         (\code{\link[stats]{binom.test}}) under the null hypothesis that the fraction of accessible or TF-covered positions
#'         in each tile equals the genome-wide expected fraction. The alternative hypothesis is that the observed fraction
#'         is greater than expected. Adjusted p-values are computed using \code{\link[stats]{p.adjust}} with \code{method="fdr"}.}
#' }

#' @importFrom GenomicRanges GRanges reduce slidingWindows
#' @importFrom stats p.adjust binom.test pnorm
#' @import data.table
#' @export
#'
calculate_tile_BG_TF_enrichments <- function(cover_dt,
                                              tile_width = 500,
                                              tile_step = 250,
                                              bg_colname = "background",
                                              tf_colname = "TF",
                                              nucl_colname = "Nucl",
                                              psc = 0.1,
                                              threads = NULL){

    mod_prob <- fragID <- seqnames <- start <- end <- width <-
        bg_score <- tf_score <- bg_score_Z <- tf_score_Z <- n_data_points <-
        bg_score_mean <- tf_score_mean <- bg_pos_cnt <- tf_pos_cnt <- tile_ID <-
        bg_prob <- bg_binom_pval <- tf_binom_pval <-
        bg_log2enr <- tf_log2enr <- bg_FDR <- tf_FDR<- NULL


    setDTthreads(threads = threads)

    assertDataTable(x = cover_dt,min.cols = 4)
    if(!all(c("seqnames","start","mod_prob",bg_colname,tf_colname,nucl_colname) %in% colnames(cover_dt))){
        stop("cover_dt must contain columns ",paste0(c("seqnames","start","mod_prob",bg_colname,tf_colname,nucl_colname),collapse = ", "))
    }


    ### keep only informative positions
    cover_dt <- cover_dt[!is.na(mod_prob)]


    ### calculate BG and TF scores (isometric log-ratio transformation)
    cover_dt <- cover_dt[,`:=`(bg_score = sqrt(2/3) * log((get(bg_colname) + psc) / sqrt((get(tf_colname) + psc) * (get(nucl_colname) + psc))),
                               tf_score = sqrt(1/2) * log((get(tf_colname) + psc)/ (get(nucl_colname) + psc))
    )]

    ### create sliding windows
    span_reg <- range(GRanges(seqnames = cover_dt[["seqnames"]],
                              IRanges(start = cover_dt[["start"]],
                                      width=1)))
    tiles_loci <- unlist(slidingWindows(span_reg,
                                        width = tile_width,
                                        step = tile_step))
    tiles_loci$tile_ID <- as.character(tiles_loci)


    ## overlap cover_dt with tiles
    cover_dt <- cover_dt[,end := start]
    setkeyv(cover_dt,cols = c("seqnames","start","end"))
    smftile_ov <- foverlaps(x = cover_dt,
                            y = data.table(as.data.frame(tiles_loci),
                                           key = c("seqnames","start","end"))[,c("width","strand") := list(NULL,NULL)],
                            which=FALSE,
                            nomatch=NULL
    )


    ## calculate mean BG/TF scores across all molecules and positions overlapping tiles

    tile_aggr_stats <- smftile_ov[,.(n_data_points = .N, ## total number of data points
                                     n_inf_pos = length(unique(i.start)),
                                     bg_score_mean = mean(bg_score,na.rm=T),
                                     tf_score_mean = mean(tf_score,na.rm=T)
    ),
    .(seqnames,start,end,tile_ID)]
    ## get mean and SD across all positions
    bg_tf_cnts <- cover_dt[,.(N_total = .N,
                              bg_score_mean_total = mean(bg_score,na.rm=T),
                              bg_score_sd_total = sd(bg_score) * sqrt((.N - 1) /.N),
                              tf_score_mean_total = mean(tf_score,na.rm=T),
                              tf_score_sd_total = sd(tf_score) * sqrt((.N - 1) /.N))]

    ## calculate enrichments and run binomial test

    tile_aggr_stats <- tile_aggr_stats[,(c(
        "bg_score_mean_Zstat",
        "bg_score_mean_Ztest_pval",
        "tf_score_mean_Zstat",
        "tf_score_mean_Ztest_pval"
    )) := {

        ## Z-test for BG/TF score means in tiles
        denom <- bg_tf_cnts$N_total - n_data_points
        size_factor <- ifelse(denom > 0,
                              sqrt(n_data_points * (bg_tf_cnts$N_total - 1) / denom),
                              NA_real_)
        bg_score_mean_zstat <- ifelse(bg_tf_cnts$bg_score_sd_total > 0,
                                      (bg_score_mean - bg_tf_cnts$bg_score_mean_total) /
                                          bg_tf_cnts$bg_score_sd_total * size_factor,
                                      NA_real_)
        bg_score_mean_zpval <- pnorm(bg_score_mean_zstat, lower.tail = FALSE)  # one-sided

        tf_score_mean_zstat <- ifelse(bg_tf_cnts$tf_score_sd_total > 0,
                                      (tf_score_mean - bg_tf_cnts$tf_score_mean_total) /
                                          bg_tf_cnts$tf_score_sd_total * size_factor,
                                      NA_real_)
        tf_score_mean_zpval <- pnorm(tf_score_mean_zstat, lower.tail = FALSE)

        list(
            bg_score_mean_zstat,
            bg_score_mean_zpval,
            tf_score_mean_zstat,
            tf_score_mean_zpval
        )
    },by = tile_ID][,`:=`(bg_score_mean_FDR = p.adjust(bg_score_mean_Ztest_pval,method = "fdr"),
                          tf_score_mean_FDR = p.adjust(tf_score_mean_Ztest_pval,method = "fdr"))]

    return(tile_aggr_stats)


}
