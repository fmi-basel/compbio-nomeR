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
#' @param bg_score_thresh Threshold for background (accessible) scores.
#' @param tf_score_thresh Threshold for TF footprint scores.
#' @param psc Pseudo-count added when calculating background and TF scores.
#' @param threads Number or threads used by \code{data.table}.
#' NULL (default) rereads environment variables. 0 means to use all logical CPUs available. Otherwise a number >= 1
#'
#' @returns A \code{data.table} with enrichment values and statistical significance calculated using a binomial test.
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
#' @importFrom stats p.adjust binom.test
#' @import data.table
#' @export
#'
calculate_tile_BG_TF_enrichments2 <- function(cover_dt,
																						 tile_width = 500,
																						 tile_step = 250,
																						 bg_colname = "background",
																						 tf_colname = "TF",
																						 nucl_colname = "Nucl",
																						 bg_score_thresh = 0.1,
																						 tf_score_thresh = bg_score_thresh,
																							bgscoreZ_thresh = 2,
																							tfscoreZ_thresh = 2,
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
	# bg_score_smpl <- sqrt(2/3) * log(bg_smpl/sqrt(tf_smpl * nucl_smpl))
	# tf_score_smpl <- sqrt(1/2) * log(tf_smpl/nucl_smpl)
	cover_dt <- cover_dt[,`:=`(bg_score = sqrt(2/3) * log((get(bg_colname) + psc) / sqrt((get(tf_colname) + psc) * (get(nucl_colname) + psc))),
														 tf_score = sqrt(1/2) * log((get(tf_colname) + psc)/ (get(nucl_colname) + psc))
	)]
	
	### calculate Z-transformed BGand TF scores
	cover_dt <- cover_dt[,':='(bg_score_Z = scale(bg_score),
														 tf_score_Z = scale(tf_score))]
	
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
	
	
	## apply thresholds, calculate aggregated statistics for each tile
	## the reason for selecting tf positives not only on tf score but also on bg score lies
	## in the distribution of points on bg_score vs tf_score scatter.
	## it looks like a triangle |> where
	## top left corner correspond to positions with TF posterior coverage near 1.
	## right corner corresponds to positions with bg posteriors near 1
	## and bottom left correspond to Nucl posteriors near 1.
	## therefore, thresholding only on tf_score may select many points from the right corner, i.e. bg positive points.
	
	tile_aggr_stats <- smftile_ov[,.(n_data_points = .N, ## total number of data points
																	 n_inf_pos = length(unique(i.start)),
																	 bg_score_mean = mean(bg_score,na.rm=T),
																	 tf_score_mean = mean(tf_score,na.rm=T),
																	 #bg_pos_cnt = sum(bg_score >= bg_score_thresh), ## number of positions with bg_score above threshold
																	 #tf_pos_cnt = sum(tf_score >= tf_score_thresh & bg_score < bg_score_thresh), ## number of positions  tf_score  above cutoff and not bg positive
																	 bg_pos_cnt = sum(bg_score_Z >= bgscoreZ_thresh),
																	 tf_pos_cnt = sum(tf_score_Z >= tfscoreZ_thresh)
	),
	.(seqnames,start,end,tile_ID)]
	#### run binomial test
	## get total numbers
	
	# bg_tf_cnts <- cover_dt[,.(bg_prob = sum(bg_score >= bg_score_thresh)/.N,
	# 													tf_prob = sum(tf_score >= tf_score_thresh & bg_score < bg_score_thresh)/.N)]
	bg_tf_cnts <- cover_dt[,.(bg_prob = sum(bg_score_Z >= bgscoreZ_thresh)/.N,
														tf_prob = sum(tf_score_Z >= tfscoreZ_thresh)/.N)]
	
	## calculate enrichments and run binomial test
	
	tile_aggr_stats <- tile_aggr_stats[,(c(
		"bg_log2enr",
		"bg_binom_pval",
		"tf_log2enr",
		"tf_binom_pval"
	)) := {
		bgbinom_res <- binom.test(
			x = bg_pos_cnt,
			n = n_data_points,
			p = bg_tf_cnts$bg_prob,
			alternative = "greater"
		)
		bg_log_enr <- unname(log2((bgbinom_res$statistic + 1)/(bgbinom_res$null.value * bgbinom_res$parameter + 1)))
		bg_pval <- bgbinom_res$p.value
		
		tfbinom_res <- binom.test(
			x = tf_pos_cnt,
			n = n_data_points,
			p = bg_tf_cnts$tf_prob,
			alternative = "greater"
		)
		tf_log_enr <- unname(log2((tfbinom_res$statistic + 1)/(tfbinom_res$null.value * tfbinom_res$parameter + 1)))
		tf_pval <- tfbinom_res$p.value
		
		list(
			bg_log_enr,
			bg_pval,
			tf_log_enr,
			tf_pval
		)
	},by = tile_ID][,`:=`(bg_FDR = p.adjust(bg_binom_pval,method = "fdr"),
												tf_FDR = p.adjust(tf_binom_pval,method = "fdr"))]
	
	return(tile_aggr_stats)
	
	
}
