## ---- helpers -----------------------------------------------------------------

make_correction_se <- function(n_pos = 15, n_frags = 8,
                               seqconts = NULL, seed = 42) {
    set.seed(seed)
    if (is.null(seqconts))
        seqconts <- sample(c("GCA","GCT","GCC","GCG"), n_pos, replace = TRUE)

    mod_vals <- matrix(runif(n_pos * n_frags), n_pos, n_frags)
    mod_vals[sample(n_pos * n_frags, 10)] <- NA

    arr <- SparseArray::NaArray(dim = c(n_pos, n_frags), type = "double")
    idx <- which(!is.na(mod_vals), arr.ind = TRUE)
    arr[idx] <- mod_vals[idx]

    assay_df <- S4Vectors::make_zero_col_DFrame(nrow = n_pos)
    assay_df[["sampleA"]] <- arr

    SummarizedExperiment::SummarizedExperiment(
        assays  = list(mod_prob = assay_df),
        rowData = S4Vectors::DataFrame(sequenceContext = seqconts))
}

## Build BetaCorrect control_params (combined data.table).
## All contexts get the same shape values for simplicity.
make_bc_params <- function(alpha_pos, beta_pos, alpha_neg, beta_neg,
                           ctxts = c("GCA","GCT","GCC","GCG","OTHER")) {
    data.table::data.table(
        seqcont   = ctxts,
        n_pos     = 1000L,
        alpha_pos = alpha_pos,
        beta_pos  = beta_pos,
        n_neg     = 1000L,
        alpha_neg = alpha_neg,
        beta_neg  = beta_neg)
}

## Build BetaUniform control_params.
## Requires MLRP: alpha_pos >= alpha_neg, beta_pos <= beta_neg.
make_bu_params <- function(alpha_pos = 8, beta_pos = 1, eps_pos = 0.05,
                           alpha_neg = 1, beta_neg = 8, eps_neg = 0.05,
                           ctxts = c("GCA","GCT","GCC","GCG","OTHER")) {
    data.table::data.table(
        seqcont   = ctxts,
        alpha_pos = alpha_pos, beta_pos = beta_pos, eps_pos = eps_pos,
        alpha_neg = alpha_neg, beta_neg = beta_neg, eps_neg = eps_neg,
        n_pos     = 1000L, n_neg = 1000L,
        converged = TRUE, mlrp_ok = TRUE)
}


## ---- .calc_beta_corrected_mod_prob_dbeta ------------------------------------

test_that(".calc_beta_corrected_mod_prob_dbeta: boundary inputs produce no NaN", {
    res <- nomeR:::.calc_beta_corrected_mod_prob_dbeta(
        mod_prob  = c(0, 1, 0.5),
        alpha_neg = 2,  beta_neg = 50,
        alpha_pos = 50, beta_pos = 2)
    expect_false(any(is.nan(res) | is.na(res)))
    expect_true(all(res >= 0 & res <= 1))
})

test_that(".calc_beta_corrected_mod_prob_dbeta: corrects in the right direction", {
    res <- nomeR:::.calc_beta_corrected_mod_prob_dbeta(
        mod_prob  = c(0.05, 0.95),
        alpha_neg = 2,  beta_neg = 50,
        alpha_pos = 50, beta_pos = 2)
    expect_lt(res[1], 0.5)
    expect_gt(res[2], 0.5)
    expect_gt(res[2], res[1])
})

test_that(".calc_beta_corrected_mod_prob_dbeta: monotone in mod_prob", {
    probs <- seq(0.05, 0.95, by = 0.1)
    res <- nomeR:::.calc_beta_corrected_mod_prob_dbeta(
        mod_prob  = probs,
        alpha_neg = 2,  beta_neg = 50,
        alpha_pos = 50, beta_pos = 2)
    expect_true(all(diff(res) >= 0))
    expect_gt(res[length(res)], res[1])
})

test_that(".calc_beta_corrected_mod_prob_dbeta: mod_prior shifts output", {
    args <- list(mod_prob  = 0.5,
                 alpha_neg = 5, beta_neg = 5,
                 alpha_pos = 5, beta_pos = 5)
    res_low  <- do.call(nomeR:::.calc_beta_corrected_mod_prob_dbeta,
                        c(args, mod_prior = 0.1))
    res_high <- do.call(nomeR:::.calc_beta_corrected_mod_prob_dbeta,
                        c(args, mod_prior = 0.9))
    expect_lt(res_low, res_high)
})

test_that(".calc_beta_corrected_mod_prob_dbeta: simultaneous underflow falls back to prior", {
    mod_prior <- 0.3
    res <- nomeR:::.calc_beta_corrected_mod_prob_dbeta(
        mod_prob  = 0.5,
        alpha_neg = 500, beta_neg = 0.1,
        alpha_pos = 500, beta_pos = 0.1,
        mod_prior = mod_prior)
    expected <- stats::plogis(log(mod_prior) - log(1 - mod_prior))
    expect_equal(res, expected, tolerance = 1e-9)
})


## ---- correct_modprob_SE: BetaCorrect ----------------------------------------

test_that("correct_modprob_SE BetaCorrect: returns SE with corrected assay added", {
    se     <- make_correction_se()
    params <- make_bc_params(50, 2, 2, 50)
    res    <- correct_modprob_SE(se, control_params = params, method = "BetaCorrect")
    expect_true("mod_prob_corrected" %in% SummarizedExperiment::assayNames(res))
    expect_equal(dim(SummarizedExperiment::assay(res, "mod_prob_corrected")),
                 dim(SummarizedExperiment::assay(res, "mod_prob")))
})

test_that("correct_modprob_SE BetaCorrect: respects corrected_assayName", {
    se     <- make_correction_se()
    params <- make_bc_params(50, 2, 2, 50)
    res    <- correct_modprob_SE(se, control_params = params, method = "BetaCorrect",
                                 corrected_assayName = "my_corrected")
    expect_true("my_corrected" %in% SummarizedExperiment::assayNames(res))
    expect_false("mod_prob_corrected" %in% SummarizedExperiment::assayNames(res))
})

test_that("correct_modprob_SE BetaCorrect: NAs in input remain NA in output", {
    se     <- make_correction_se()
    params <- make_bc_params(50, 2, 2, 50)
    res    <- correct_modprob_SE(se, control_params = params, method = "BetaCorrect")
    raw_nonNA  <- SparseArray::nnawhich(
        SummarizedExperiment::assay(se,  "mod_prob")[["sampleA"]], arr.ind = TRUE)
    corr_nonNA <- SparseArray::nnawhich(
        SummarizedExperiment::assay(res, "mod_prob_corrected")[["sampleA"]], arr.ind = TRUE)
    expect_equal(raw_nonNA, corr_nonNA)
})

test_that("correct_modprob_SE BetaCorrect: corrected values are in [0, 1]", {
    se     <- make_correction_se()
    params <- make_bc_params(50, 2, 2, 50)
    res    <- correct_modprob_SE(se, control_params = params, method = "BetaCorrect")
    idx  <- SparseArray::nnawhich(
        SummarizedExperiment::assay(res, "mod_prob_corrected")[["sampleA"]],
        arr.ind = TRUE)
    vals <- SummarizedExperiment::assay(res, "mod_prob_corrected")[["sampleA"]][idx]
    expect_true(all(vals >= 0 & vals <= 1))
})

test_that("correct_modprob_SE BetaCorrect: errors if sequenceContext missing", {
    se  <- make_correction_se()
    SummarizedExperiment::rowData(se)$sequenceContext <- NULL
    params <- make_bc_params(50, 2, 2, 50)
    expect_error(correct_modprob_SE(se, control_params = params, method = "BetaCorrect"))
})

test_that("correct_modprob_SE BetaCorrect: unknown context uses OTHER without error", {
    se     <- make_correction_se(seqconts = rep("UNKNOWN", 15))
    params <- make_bc_params(50, 2, 2, 50)
    expect_no_error(
        suppressWarnings(
            correct_modprob_SE(se, control_params = params, method = "BetaCorrect")))
})

test_that("correct_modprob_SE BetaCorrect: shifts high/low probs in correct direction", {
    set.seed(1)
    n_pos <- 50; n_frags <- 10
    hi  <- runif(n_pos * n_frags / 2, 0.8, 1.0)
    lo  <- runif(n_pos * n_frags / 2, 0.0, 0.2)
    raw <- c(hi, lo)
    arr <- SparseArray::NaArray(dim = c(n_pos, n_frags), type = "double")
    arr[cbind(rep(seq_len(n_pos), n_frags),
              rep(seq_len(n_frags), each = n_pos))] <- raw
    adf <- S4Vectors::make_zero_col_DFrame(nrow = n_pos)
    adf[["s1"]] <- arr
    se <- SummarizedExperiment::SummarizedExperiment(
        assays  = list(mod_prob = adf),
        rowData = S4Vectors::DataFrame(sequenceContext = rep("GCA", n_pos)))
    params <- make_bc_params(50, 2, 2, 50)
    res    <- correct_modprob_SE(se, control_params = params, method = "BetaCorrect")
    corr   <- as.vector(
        SummarizedExperiment::assay(res, "mod_prob_corrected")[["s1"]][
            cbind(rep(seq_len(n_pos), n_frags),
                  rep(seq_len(n_frags), each = n_pos))])
    expect_gt(mean(corr[raw > 0.8]), mean(raw[raw > 0.8]))
    expect_lt(mean(corr[raw < 0.2]), mean(raw[raw < 0.2]))
})

test_that("correct_modprob_SE BetaCorrect: qnorm_to_raw maps to raw value set", {
    se     <- make_correction_se(seed = 7)
    params <- make_bc_params(50, 2, 2, 50)
    res    <- correct_modprob_SE(se, control_params = params, method = "BetaCorrect",
                                 qnorm_to_raw = TRUE)
    idx <- SparseArray::nnawhich(
        SummarizedExperiment::assay(se, "mod_prob")[["sampleA"]], arr.ind = TRUE)
    raw_vals  <- SummarizedExperiment::assay(se,  "mod_prob")[["sampleA"]][idx]
    corr_vals <- SummarizedExperiment::assay(res, "mod_prob_corrected")[["sampleA"]][idx]
    expect_true(all(corr_vals %in% raw_vals))
})

test_that("correct_modprob_SE BetaCorrect: handles multiple samples independently", {
    set.seed(99)
    n_pos <- 10; n_frags <- 5
    make_arr <- function() {
        arr <- SparseArray::NaArray(dim = c(n_pos, n_frags), type = "double")
        arr[cbind(rep(seq_len(n_pos), n_frags),
                  rep(seq_len(n_frags), each = n_pos))] <- runif(n_pos * n_frags)
        arr
    }
    adf <- S4Vectors::make_zero_col_DFrame(nrow = n_pos)
    adf[["s1"]] <- make_arr(); adf[["s2"]] <- make_arr()
    se <- SummarizedExperiment::SummarizedExperiment(
        assays  = list(mod_prob = adf),
        rowData = S4Vectors::DataFrame(sequenceContext = rep("GCA", n_pos)))
    params <- make_bc_params(50, 2, 2, 50)
    res    <- correct_modprob_SE(se, control_params = params, method = "BetaCorrect")
    corr_s1 <- SummarizedExperiment::assay(res, "mod_prob_corrected")[["s1"]]
    corr_s2 <- SummarizedExperiment::assay(res, "mod_prob_corrected")[["s2"]]
    expect_false(identical(corr_s1, corr_s2))
    expect_equal(ncol(res), 2L)
})


## ---- correct_modprob_SE: BetaUniform ----------------------------------------

test_that("correct_modprob_SE BetaUniform: returns SE with corrected assay added", {
    se     <- make_correction_se()
    params <- make_bu_params()
    res    <- correct_modprob_SE(se, control_params = params, method = "BetaUniform")
    expect_true("mod_prob_corrected" %in% SummarizedExperiment::assayNames(res))
    expect_equal(dim(SummarizedExperiment::assay(res, "mod_prob_corrected")),
                 dim(SummarizedExperiment::assay(res, "mod_prob")))
})

test_that("correct_modprob_SE BetaUniform: corrected values are in [0, 1]", {
    se     <- make_correction_se()
    params <- make_bu_params()
    res    <- correct_modprob_SE(se, control_params = params, method = "BetaUniform")
    idx  <- SparseArray::nnawhich(
        SummarizedExperiment::assay(res, "mod_prob_corrected")[["sampleA"]],
        arr.ind = TRUE)
    vals <- SummarizedExperiment::assay(res, "mod_prob_corrected")[["sampleA"]][idx]
    expect_true(all(vals >= 0 & vals <= 1))
})

test_that("correct_modprob_SE BetaUniform: NAs in input remain NA in output", {
    se     <- make_correction_se()
    params <- make_bu_params()
    res    <- correct_modprob_SE(se, control_params = params, method = "BetaUniform")
    raw_nonNA  <- SparseArray::nnawhich(
        SummarizedExperiment::assay(se,  "mod_prob")[["sampleA"]], arr.ind = TRUE)
    corr_nonNA <- SparseArray::nnawhich(
        SummarizedExperiment::assay(res, "mod_prob_corrected")[["sampleA"]], arr.ind = TRUE)
    expect_equal(raw_nonNA, corr_nonNA)
})

test_that("correct_modprob_SE BetaUniform: monotone (corrected increases with raw, per context)", {
    set.seed(5)
    n_pos <- 100; n_frags <- 1
    raw <- sort(runif(n_pos))
    arr <- SparseArray::NaArray(dim = c(n_pos, n_frags), type = "double")
    arr[cbind(seq_len(n_pos), rep(1L, n_pos))] <- raw
    adf <- S4Vectors::make_zero_col_DFrame(nrow = n_pos)
    adf[["s1"]] <- arr
    se <- SummarizedExperiment::SummarizedExperiment(
        assays  = list(mod_prob = adf),
        rowData = S4Vectors::DataFrame(sequenceContext = rep("GCA", n_pos)))
    params <- make_bu_params()
    res    <- correct_modprob_SE(se, control_params = params, method = "BetaUniform")
    corr   <- SummarizedExperiment::assay(res, "mod_prob_corrected")[["s1"]][
        cbind(seq_len(n_pos), rep(1L, n_pos))]
    expect_true(all(diff(corr) >= -1e-10))
})

test_that("correct_modprob_SE BetaUniform: shifts high/low probs in correct direction", {
    set.seed(2)
    n_pos <- 50; n_frags <- 10
    hi  <- runif(n_pos * n_frags / 2, 0.8, 1.0)
    lo  <- runif(n_pos * n_frags / 2, 0.0, 0.2)
    raw <- c(hi, lo)
    arr <- SparseArray::NaArray(dim = c(n_pos, n_frags), type = "double")
    arr[cbind(rep(seq_len(n_pos), n_frags),
              rep(seq_len(n_frags), each = n_pos))] <- raw
    adf <- S4Vectors::make_zero_col_DFrame(nrow = n_pos)
    adf[["s1"]] <- arr
    se <- SummarizedExperiment::SummarizedExperiment(
        assays  = list(mod_prob = adf),
        rowData = S4Vectors::DataFrame(sequenceContext = rep("GCA", n_pos)))
    params <- make_bu_params()
    res    <- correct_modprob_SE(se, control_params = params, method = "BetaUniform")
    corr   <- as.vector(
        SummarizedExperiment::assay(res, "mod_prob_corrected")[["s1"]][
            cbind(rep(seq_len(n_pos), n_frags),
                  rep(seq_len(n_frags), each = n_pos))])
    expect_gt(mean(corr[raw > 0.8]), mean(raw[raw > 0.8]))
    expect_lt(mean(corr[raw < 0.2]), mean(raw[raw < 0.2]))
})

test_that("correct_modprob_SE BetaUniform: qnorm_to_raw runs without warning", {
    se     <- make_correction_se()
    params <- make_bu_params()
    expect_no_warning(
        correct_modprob_SE(se, control_params = params, method = "BetaUniform",
                           qnorm_to_raw = TRUE))
})

test_that("correct_modprob_SE BetaUniform: errors if missing BU columns", {
    se     <- make_correction_se()
    params <- make_bc_params(50, 2, 2, 50)  # BC params, not BU
    expect_error(
        correct_modprob_SE(se, control_params = params, method = "BetaUniform"))
})
