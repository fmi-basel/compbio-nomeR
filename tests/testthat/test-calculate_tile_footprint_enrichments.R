## Tests for calculate_tile_footprint_enrichments

## ── shared fixture ────────────────────────────────────────────────────────────

.ctfe_dt <- local({
    set.seed(42)
    n <- 300
    data.table::data.table(
        seqnames   = "chr1",
        start      = seq_len(n),
        mod_prob   = runif(n, 0.1, 0.9),
        background = runif(n, 0.2, 0.6),
        TF         = runif(n, 0.05, 0.25),
        Nucl       = runif(n, 0.30, 0.60)
    )
})

## ── input validation ──────────────────────────────────────────────────────────

test_that("calculate_tile_footprint_enrichments: errors on non-data.table input", {
    expect_error(
        calculate_tile_footprint_enrichments(as.data.frame(.ctfe_dt),
                                             ftp_names = "background")
    )
})

test_that("calculate_tile_footprint_enrichments: errors when required columns are missing", {
    dt_no_seqnames <- data.table::copy(.ctfe_dt)[, seqnames := NULL]
    expect_error(
        calculate_tile_footprint_enrichments(dt_no_seqnames, ftp_names = "background")
    )

    dt_no_start <- data.table::copy(.ctfe_dt)[, start := NULL]
    expect_error(
        calculate_tile_footprint_enrichments(dt_no_start, ftp_names = "background")
    )

    dt_no_mod <- data.table::copy(.ctfe_dt)[, mod_prob := NULL]
    expect_error(
        calculate_tile_footprint_enrichments(dt_no_mod, ftp_names = "background")
    )
})

test_that("calculate_tile_footprint_enrichments: errors when ftp_names column is absent", {
    expect_error(
        calculate_tile_footprint_enrichments(.ctfe_dt, ftp_names = "nonexistent_col")
    )
})

test_that("calculate_tile_footprint_enrichments: errors on invalid ftp_names type", {
    expect_error(
        calculate_tile_footprint_enrichments(.ctfe_dt, ftp_names = 123)
    )
    expect_error(
        calculate_tile_footprint_enrichments(.ctfe_dt, ftp_names = character(0))
    )
})

## ── output structure ──────────────────────────────────────────────────────────

test_that("calculate_tile_footprint_enrichments: returns a data.table", {
    result <- calculate_tile_footprint_enrichments(.ctfe_dt,
                                                   tile_width = 100,
                                                   tile_step  = 50,
                                                   ftp_names  = "background")
    expect_s3_class(result, "data.table")
})

test_that("calculate_tile_footprint_enrichments: single ftp produces expected columns", {
    result <- calculate_tile_footprint_enrichments(.ctfe_dt,
                                                   tile_width = 100,
                                                   tile_step  = 50,
                                                   ftp_names  = "background")


    expected_cols <- c("seqnames", "start", "end", "tile_ID",
                       "n_data_points", "n_inf_pos",
                       "background_score_mean",
                       "background_score_mean_Zstat",
                       "background_score_mean_Ztest_pval",
                       "background_score_mean_FDR")
    missing <- setdiff(expected_cols, names(result))
    expect_equal(length(missing), 0L,
                 info = paste("Missing columns:", paste(missing, collapse = ", ")))
})

test_that("calculate_tile_footprint_enrichments: multiple ftps each produce a full column set", {
    result <- calculate_tile_footprint_enrichments(.ctfe_dt,
                                                   tile_width = 100,
                                                   tile_step  = 50,
                                                   ftp_names  = c("background", "TF", "Nucl"))

    for (ftp in c("background", "TF", "Nucl")) {
        expected <- paste0(ftp, "_score_mean",
                           c("", "_Zstat", "_Ztest_pval", "_FDR"))
        missing <- setdiff(expected, names(result))
        expect_equal(length(missing), 0L,
                     info = paste(ftp, "- missing:", paste(missing, collapse = ", ")))
    }
    ## no unexpected columns from other footprints
    expect_false("nonexistent_score_mean" %in% names(result))
})

test_that("calculate_tile_footprint_enrichments: p-values and FDR are in [0, 1]", {
    result <- calculate_tile_footprint_enrichments(.ctfe_dt,
                                                   tile_width = 100,
                                                   tile_step  = 50,
                                                   ftp_names  = c("background", "TF"))

    for (col in c("background_score_mean_Ztest_pval", "background_score_mean_FDR",
                  "TF_score_mean_Ztest_pval",         "TF_score_mean_FDR")) {
        vals <- result[[col]]
        expect_true(all(vals[!is.na(vals)] >= 0 & vals[!is.na(vals)] <= 1),
                    info = paste(col, "out of [0,1]"))
    }
})

test_that("calculate_tile_footprint_enrichments: FDR >= raw p-value (BH non-decreasing)", {
    result <- calculate_tile_footprint_enrichments(.ctfe_dt,
                                                   tile_width = 100,
                                                   tile_step  = 50,
                                                   ftp_names  = c("background", "TF"))

    for (ftp in c("background", "TF")) {
        pval <- result[[paste0(ftp, "_score_mean_Ztest_pval")]]
        fdr  <- result[[paste0(ftp, "_score_mean_FDR")]]
        ok   <- !is.na(pval) & !is.na(fdr)
        expect_true(all(fdr[ok] >= pval[ok] - 1e-10),
                    info = paste(ftp, "FDR < pval"))
    }
})

test_that("calculate_tile_footprint_enrichments: n_data_points >= n_inf_pos", {
    result <- calculate_tile_footprint_enrichments(.ctfe_dt,
                                                   tile_width = 100,
                                                   tile_step  = 50,
                                                   ftp_names  = "background")

    expect_true(all(result$n_data_points > 0L))
    expect_true(all(result$n_inf_pos     > 0L))
    expect_true(all(result$n_data_points >= result$n_inf_pos))
})

## ── NA handling ───────────────────────────────────────────────────────────────

test_that("calculate_tile_footprint_enrichments: NA mod_prob rows excluded", {
    dt_na <- data.table::copy(.ctfe_dt)
    data.table::set(dt_na, i = 1:30, j = "mod_prob", value = NA_real_)

    result_na    <- calculate_tile_footprint_enrichments(dt_na,
                                                         tile_width = 1000,
                                                         tile_step  = 1000,
                                                         ftp_names  = "background")
    result_clean <- calculate_tile_footprint_enrichments(dt_na[!is.na(mod_prob)],
                                                         tile_width = 1000,
                                                         tile_step  = 1000,
                                                         ftp_names  = "background")

    expect_equal(result_na$n_data_points, result_clean$n_data_points)
    expect_equal(result_na$n_inf_pos,     result_clean$n_inf_pos)
})

## ── tile geometry ─────────────────────────────────────────────────────────────

test_that("calculate_tile_footprint_enrichments: single tile when tile_width spans all data", {
    result <- calculate_tile_footprint_enrichments(.ctfe_dt,
                                                   tile_width = 1000,
                                                   tile_step  = 1000,
                                                   ftp_names  = "background")

    expect_equal(nrow(result), 1L)
    expect_equal(result$n_data_points, nrow(.ctfe_dt))
})

test_that("calculate_tile_footprint_enrichments: non-overlapping tiles sum to total data points", {
    result <- calculate_tile_footprint_enrichments(.ctfe_dt,
                                                   tile_width = 150,
                                                   tile_step  = 150,
                                                   ftp_names  = "background")

    expect_equal(sum(result$n_data_points), nrow(.ctfe_dt))
})

test_that("calculate_tile_footprint_enrichments: overlapping tiles inflate total data-point sum", {
    result <- calculate_tile_footprint_enrichments(.ctfe_dt,
                                                   tile_width = 150,
                                                   tile_step  = 75,
                                                   ftp_names  = "background")

    expect_gt(sum(result$n_data_points), nrow(.ctfe_dt))
})

## ── enrichment signal ─────────────────────────────────────────────────────────

test_that("calculate_tile_footprint_enrichments: enriched region has highest Z-stat", {
    set.seed(77)
    n <- 600
    ## First 200 bp have strongly elevated background coverage
    dt_signal <- data.table::data.table(
        seqnames   = "chr1",
        start      = seq_len(n),
        mod_prob   = runif(n, 0.1, 0.9),
        background = c(runif(200, 0.85, 1.0), runif(400, 0.00, 0.10)),
        TF         = runif(n, 0.05, 0.20)
    )

    result <- calculate_tile_footprint_enrichments(dt_signal,
                                                   tile_width = 200,
                                                   tile_step  = 200,
                                                   ftp_names  = "background")

    expect_equal(which.max(result$background_score_mean_Zstat), 1L)
    expect_lt(result$background_score_mean_Ztest_pval[1], 0.05)
})

test_that("calculate_tile_footprint_enrichments: TF enriched region detected independently of background", {
    set.seed(88)
    n <- 600
    ## First 200 bp have elevated TF; background is uniform
    dt_signal <- data.table::data.table(
        seqnames   = "chr1",
        start      = seq_len(n),
        mod_prob   = runif(n, 0.1, 0.9),
        background = runif(n, 0.3, 0.6),
        TF         = c(runif(200, 0.80, 1.0), runif(400, 0.00, 0.10))
    )

    result <- calculate_tile_footprint_enrichments(dt_signal,
                                                   tile_width = 200,
                                                   tile_step  = 200,
                                                   ftp_names  = c("background", "TF"))

    expect_equal(which.max(result$TF_score_mean_Zstat), 1L)
    expect_lt(result$TF_score_mean_Ztest_pval[1], 0.05)
})

## ── Z-stat independence ───────────────────────────────────────────────────────

test_that("calculate_tile_footprint_enrichments: TF Z-stats near zero when TF coverage is constant", {
    set.seed(99)
    n <- 600
    ## Background varies strongly; TF is near-constant everywhere
    dt_indep <- data.table::data.table(
        seqnames   = "chr1",
        start      = seq_len(n),
        mod_prob   = runif(n, 0.1, 0.9),
        background = c(runif(200, 0.85, 1.0), runif(400, 0.00, 0.10)),
        TF         = rep(0.10, n)   ## constant → sd = 0 → Zstat should be NA
    )

    result <- calculate_tile_footprint_enrichments(dt_indep,
                                                   tile_width = 200,
                                                   tile_step  = 200,
                                                   ftp_names  = c("background", "TF"))

    ## background Z-stat for tile 1 must be strongly positive
    expect_gt(result$background_score_mean_Zstat[1], 3)

    ## TF has zero variance → Zstat is NA for all tiles
    expect_true(all(is.na(result$TF_score_mean_Zstat)))
})

test_that("calculate_tile_footprint_enrichments: BG and TF Z-stats uncorrelated when signals differ", {
    n <- 800
    ## Deterministic: BG high in positions 1-200; TF high in positions 401-600; rest flat
    dt_mixed <- data.table::data.table(
        seqnames   = "chr1",
        start      = seq_len(n),
        mod_prob   = rep(0.5, n),
        background = c(rep(0.95, 200), rep(0.05, 600)),
        TF         = c(rep(0.05, 400), rep(0.95, 200), rep(0.05, 200))
    )

    result <- calculate_tile_footprint_enrichments(dt_mixed,
                                                   tile_width = 200,
                                                   tile_step  = 200,
                                                   ftp_names  = c("background", "TF"))

    ## filter tiles by coordinate
    tile1 <- result[start == 1]
    tile3 <- result[start == 401]

    ## bg signal in tile 1: should be strongly positive
    expect_gt(tile1$background_score_mean_Zstat, 3)
    ## TF signal in tile 3: should be strongly positive
    expect_gt(tile3$TF_score_mean_Zstat, 3)
    ## no TF signal in tile 1: should be negative
    expect_lt(tile1$TF_score_mean_Zstat, 0)
    ## no bg signal in tile 3: should be negative
    expect_lt(tile3$background_score_mean_Zstat, 0)
})

## ── psc parameter ─────────────────────────────────────────────────────────────

test_that("calculate_tile_footprint_enrichments: larger psc shrinks score range", {
    result_small_psc <- calculate_tile_footprint_enrichments(.ctfe_dt,
                                                              tile_width = 100,
                                                              tile_step  = 50,
                                                              ftp_names  = "background",
                                                              psc        = 1e-6)
    result_large_psc <- calculate_tile_footprint_enrichments(.ctfe_dt,
                                                              tile_width = 100,
                                                              tile_step  = 50,
                                                              ftp_names  = "background",
                                                              psc        = 10)

    range_small <- diff(range(result_small_psc$background_score_mean, na.rm = TRUE))
    range_large <- diff(range(result_large_psc$background_score_mean, na.rm = TRUE))

    expect_gt(range_small, range_large)
})
