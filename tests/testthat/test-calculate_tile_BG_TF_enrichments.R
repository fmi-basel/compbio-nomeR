## Tests for calculate_tile_BG_TF_enrichments

## ── shared fixture ────────────────────────────────────────────────────────────

.ctbe_dt <- local({
    set.seed(101)
    n <- 300
    data.table::data.table(
        seqnames   = "chr1",
        start      = seq_len(n),
        mod_prob   = runif(n, 0.1, 0.9),
        background = runif(n, 0.2, 0.6),
        TF         = runif(n, 0.1, 0.4),
        Nucl       = runif(n, 0.3, 0.7)
    )
})

## ── input validation ──────────────────────────────────────────────────────────

test_that("calculate_tile_BG_TF_enrichments: errors on non-data.table input", {
    expect_error(calculate_tile_BG_TF_enrichments(as.data.frame(.ctbe_dt)))
})

test_that("calculate_tile_BG_TF_enrichments: errors when required columns are missing", {
    dt_no_seqnames <- data.table::copy(.ctbe_dt)[, seqnames := NULL]
    expect_error(calculate_tile_BG_TF_enrichments(dt_no_seqnames))

    dt_no_bg <- data.table::copy(.ctbe_dt)[, background := NULL]
    expect_error(calculate_tile_BG_TF_enrichments(dt_no_bg))

    dt_no_tf <- data.table::copy(.ctbe_dt)[, TF := NULL]
    expect_error(calculate_tile_BG_TF_enrichments(dt_no_tf))
})

## ── output structure ──────────────────────────────────────────────────────────

test_that("calculate_tile_BG_TF_enrichments: returns a data.table", {
    result <- calculate_tile_BG_TF_enrichments(.ctbe_dt, tile_width = 100, tile_step = 50)
    expect_s3_class(result, "data.table")
})

test_that("calculate_tile_BG_TF_enrichments: output contains all expected columns", {
    result <- calculate_tile_BG_TF_enrichments(.ctbe_dt, tile_width = 100, tile_step = 50)

    expected_cols <- c("seqnames", "start", "end", "tile_ID",
                       "n_data_points", "n_inf_pos",
                       "bg_score_mean", "tf_score_mean",
                       "bg_score_mean_Zstat",    "bg_score_mean_Ztest_pval", "bg_score_mean_FDR",
                       "tf_score_mean_Zstat",    "tf_score_mean_Ztest_pval", "tf_score_mean_FDR")
    missing <- setdiff(expected_cols, names(result))
    expect_equal(length(missing), 0L,
                 info = paste("Missing columns:", paste(missing, collapse = ", ")))
})

test_that("calculate_tile_BG_TF_enrichments: p-values and FDR are in [0, 1]", {
    result <- calculate_tile_BG_TF_enrichments(.ctbe_dt, tile_width = 100, tile_step = 50)

    expect_true(all(result$bg_score_mean_Ztest_pval >= 0 & result$bg_score_mean_Ztest_pval <= 1, na.rm = TRUE))
    expect_true(all(result$tf_score_mean_Ztest_pval >= 0 & result$tf_score_mean_Ztest_pval <= 1, na.rm = TRUE))
    expect_true(all(result$bg_score_mean_FDR        >= 0 & result$bg_score_mean_FDR        <= 1, na.rm = TRUE))
    expect_true(all(result$tf_score_mean_FDR        >= 0 & result$tf_score_mean_FDR        <= 1, na.rm = TRUE))
})

test_that("calculate_tile_BG_TF_enrichments: FDR >= raw p-value (BH adjustment inflates)", {
    result <- calculate_tile_BG_TF_enrichments(.ctbe_dt, tile_width = 100, tile_step = 50)

    ok_bg <- !is.na(result$bg_score_mean_FDR) & !is.na(result$bg_score_mean_Ztest_pval)
    expect_true(all(result$bg_score_mean_FDR[ok_bg] >= result$bg_score_mean_Ztest_pval[ok_bg] - 1e-10))

    ok_tf <- !is.na(result$tf_score_mean_FDR) & !is.na(result$tf_score_mean_Ztest_pval)
    expect_true(all(result$tf_score_mean_FDR[ok_tf] >= result$tf_score_mean_Ztest_pval[ok_tf] - 1e-10))
})

test_that("calculate_tile_BG_TF_enrichments: n_data_points and n_inf_pos are positive integers", {
    result <- calculate_tile_BG_TF_enrichments(.ctbe_dt, tile_width = 100, tile_step = 50)

    expect_true(all(result$n_data_points > 0L))
    expect_true(all(result$n_inf_pos     > 0L))
    ## coverage (n_data_points / n_inf_pos) must be >= 1
    expect_true(all(result$n_data_points >= result$n_inf_pos))
})

## ── NA handling ───────────────────────────────────────────────────────────────

test_that("calculate_tile_BG_TF_enrichments: NA mod_prob rows are excluded from tile counts", {
    dt_na <- data.table::copy(.ctbe_dt)
    data.table::set(dt_na, i = 1:20, j = "mod_prob", value = NA_real_)

    result_na    <- calculate_tile_BG_TF_enrichments(dt_na,                       tile_width = 1000, tile_step = 1000)
    result_clean <- calculate_tile_BG_TF_enrichments(dt_na[!is.na(mod_prob)],    tile_width = 1000, tile_step = 1000)

    expect_equal(result_na$n_data_points, result_clean$n_data_points)
    expect_equal(result_na$n_inf_pos,     result_clean$n_inf_pos)
})

## ── tile geometry ─────────────────────────────────────────────────────────────

test_that("calculate_tile_BG_TF_enrichments: single tile when tile_width exceeds data range", {
    result <- calculate_tile_BG_TF_enrichments(.ctbe_dt, tile_width = 1000, tile_step = 1000)

    expect_equal(nrow(result), 1L)
    expect_equal(result$n_data_points, nrow(.ctbe_dt))
})

test_that("calculate_tile_BG_TF_enrichments: non-overlapping tiles sum to total data points", {
    ## tile_width == tile_step and tile_width == range → no overlap, one tile
    result <- calculate_tile_BG_TF_enrichments(.ctbe_dt, tile_width = 150, tile_step = 150)

    expect_equal(sum(result$n_data_points), nrow(.ctbe_dt))
})

test_that("calculate_tile_BG_TF_enrichments: overlapping tiles inflate total data-point sum", {
    result <- calculate_tile_BG_TF_enrichments(.ctbe_dt, tile_width = 150, tile_step = 75)

    ## overlapping tiles count each position more than once
    expect_gt(sum(result$n_data_points), nrow(.ctbe_dt))
})

## ── enrichment signal ─────────────────────────────────────────────────────────

test_that("calculate_tile_BG_TF_enrichments: enriched tile has the highest bg Z-stat", {
    set.seed(202)
    n <- 600
    ## First 200 bp are strongly bg-enriched; remaining 400 bp are low
    dt_signal <- data.table::data.table(
        seqnames   = "chr1",
        start      = seq_len(n),
        mod_prob   = runif(n, 0.1, 0.9),
        background = c(runif(200, 0.85, 1.0), runif(400, 0.0, 0.15)),
        TF         = runif(n, 0.05, 0.20),
        Nucl       = runif(n, 0.20, 0.40)
    )

    result <- calculate_tile_BG_TF_enrichments(dt_signal, tile_width = 200, tile_step = 200)

    ## tile containing positions 1-200 should have the highest bg Z-stat
    expect_equal(which.max(result$bg_score_mean_Zstat), 1L)
    expect_lt(result$bg_score_mean_Ztest_pval[1], 0.05)
})

test_that("calculate_tile_BG_TF_enrichments: enriched tile has the highest tf Z-stat", {
    set.seed(303)
    n <- 600
    ## First 200 bp are strongly tf-enriched; remaining 400 bp are low
    dt_signal <- data.table::data.table(
        seqnames   = "chr1",
        start      = seq_len(n),
        mod_prob   = runif(n, 0.1, 0.9),
        background = runif(n, 0.2, 0.5),
        TF         = c(runif(200, 0.80, 1.0), runif(400, 0.0, 0.10)),
        Nucl       = runif(n, 0.05, 0.15)
    )

    result <- calculate_tile_BG_TF_enrichments(dt_signal, tile_width = 200, tile_step = 200)

    expect_equal(which.max(result$tf_score_mean_Zstat), 1L)
    expect_lt(result$tf_score_mean_Ztest_pval[1], 0.05)
})

## ── custom column names ───────────────────────────────────────────────────────

test_that("calculate_tile_BG_TF_enrichments: works with custom bg/tf/nucl column names", {
    dt_custom <- data.table::copy(.ctbe_dt)
    data.table::setnames(dt_custom,
                         c("background", "TF",     "Nucl"),
                         c("bg_cov",     "tf_cov", "nucl_cov"))

    result <- calculate_tile_BG_TF_enrichments(dt_custom,
                                               tile_width   = 100,
                                               tile_step    = 50,
                                               bg_colname   = "bg_cov",
                                               tf_colname   = "tf_cov",
                                               nucl_colname = "nucl_cov")

    expect_s3_class(result, "data.table")
    expect_gt(nrow(result), 0L)
    expect_true("bg_score_mean" %in% names(result))
})
