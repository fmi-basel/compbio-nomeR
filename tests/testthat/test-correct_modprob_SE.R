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

make_shapes <- function(shape1, shape2,
                        ctxts = c("GCA","GCT","GCC","GCG","other")) {
    data.table::data.table(
        seqcont = ctxts,
        n_dat   = 1000L,
        shape1  = shape1,
        shape2  = shape2)
}

## ---- .calc_beta_corrected_mod_prob_dbeta ------------------------------------

test_that(".calc_beta_corrected_mod_prob_dbeta: boundary inputs produce no NaN", {
    res <- footBayes:::.calc_beta_corrected_mod_prob_dbeta(
        mod_prob   = c(0, 1, 0.5),
        neg_shape1 = 2, neg_shape2 = 50,
        pos_shape1 = 50, pos_shape2 = 2)
    expect_false(any(is.nan(res) | is.na(res)))
    expect_true(all(res >= 0 & res <= 1))
})

test_that(".calc_beta_corrected_mod_prob_dbeta: corrects in the right direction", {
    res <- footBayes:::.calc_beta_corrected_mod_prob_dbeta(
        mod_prob   = c(0.05, 0.95),
        neg_shape1 = 2, neg_shape2 = 50,
        pos_shape1 = 50, pos_shape2 = 2)
    expect_lt(res[1], 0.5)
    expect_gt(res[2], 0.5)
    expect_gt(res[2], res[1])
})

test_that(".calc_beta_corrected_mod_prob_dbeta: monotone in mod_prob", {
    probs <- seq(0.05, 0.95, by = 0.1)
    res <- footBayes:::.calc_beta_corrected_mod_prob_dbeta(
        mod_prob   = probs,
        neg_shape1 = 2, neg_shape2 = 50,
        pos_shape1 = 50, pos_shape2 = 2)
    # plogis saturates to exactly 1 at high log-odds, so allow ties at the boundary
    expect_true(all(diff(res) >= 0))
    expect_gt(res[length(res)], res[1])
})

test_that(".calc_beta_corrected_mod_prob_dbeta: mod_prior shifts output", {
    args <- list(mod_prob   = 0.5,
                 neg_shape1 = 5, neg_shape2 = 5,
                 pos_shape1 = 5, pos_shape2 = 5)
    res_low  <- do.call(footBayes:::.calc_beta_corrected_mod_prob_dbeta,
                        c(args, mod_prior = 0.1))
    res_high <- do.call(footBayes:::.calc_beta_corrected_mod_prob_dbeta,
                        c(args, mod_prior = 0.9))
    expect_lt(res_low, res_high)
})

test_that(".calc_beta_corrected_mod_prob_dbeta: simultaneous underflow falls back to prior", {
    mod_prior <- 0.3
    # Beta(500, 0.1) concentrates all mass so near 1 that dbeta(0.5, ...) underflows
    # for both pos and neg shapes → NaN path → output should equal plogis(log-prior-odds)
    res <- footBayes:::.calc_beta_corrected_mod_prob_dbeta(
        mod_prob   = 0.5,
        neg_shape1 = 500, neg_shape2 = 0.1,
        pos_shape1 = 500, pos_shape2 = 0.1,
        mod_prior  = mod_prior)
    expected <- stats::plogis(log(mod_prior) - log(1 - mod_prior))
    expect_equal(res, expected, tolerance = 1e-9)
})

## ---- correct_modprob_SE -----------------------------------------------------

test_that("correct_modprob_SE: returns SE with corrected assay added", {
    se  <- make_correction_se()
    neg <- make_shapes(2, 50)
    pos <- make_shapes(50, 2)
    res <- correct_modprob_SE(se,
                              neg_control_shapes = neg,
                              pos_control_shapes = pos)
    expect_true("mod_prob_corrected" %in% SummarizedExperiment::assayNames(res))
    expect_equal(dim(SummarizedExperiment::assay(res, "mod_prob_corrected")),
                 dim(SummarizedExperiment::assay(res, "mod_prob")))
})

test_that("correct_modprob_SE: respects corrected_assayName argument", {
    se  <- make_correction_se()
    neg <- make_shapes(2, 50)
    pos <- make_shapes(50, 2)
    res <- correct_modprob_SE(se,
                              neg_control_shapes  = neg,
                              pos_control_shapes  = pos,
                              corrected_assayName = "my_corrected")
    expect_true("my_corrected" %in% SummarizedExperiment::assayNames(res))
    expect_false("mod_prob_corrected" %in% SummarizedExperiment::assayNames(res))
})

test_that("correct_modprob_SE: NAs in input remain NA in output", {
    se  <- make_correction_se()
    neg <- make_shapes(2, 50)
    pos <- make_shapes(50, 2)
    res <- correct_modprob_SE(se,
                              neg_control_shapes = neg,
                              pos_control_shapes = pos)
    raw_nonNA  <- SparseArray::nnawhich(
        SummarizedExperiment::assay(se,  "mod_prob")[["sampleA"]], arr.ind = TRUE)
    corr_nonNA <- SparseArray::nnawhich(
        SummarizedExperiment::assay(res, "mod_prob_corrected")[["sampleA"]], arr.ind = TRUE)
    expect_equal(raw_nonNA, corr_nonNA)
})

test_that("correct_modprob_SE: corrected values are in [0, 1]", {
    se  <- make_correction_se()
    neg <- make_shapes(2, 50)
    pos <- make_shapes(50, 2)
    res <- correct_modprob_SE(se,
                              neg_control_shapes = neg,
                              pos_control_shapes = pos)
    idx  <- SparseArray::nnawhich(
        SummarizedExperiment::assay(res, "mod_prob_corrected")[["sampleA"]],
        arr.ind = TRUE)
    vals <- SummarizedExperiment::assay(res, "mod_prob_corrected")[["sampleA"]][idx]
    expect_true(all(vals >= 0 & vals <= 1))
})

test_that("correct_modprob_SE: errors if sequenceContext is missing from rowData", {
    se  <- make_correction_se()
    SummarizedExperiment::rowData(se)$sequenceContext <- NULL
    neg <- make_shapes(2, 50)
    pos <- make_shapes(50, 2)
    expect_error(correct_modprob_SE(se,
                                    neg_control_shapes = neg,
                                    pos_control_shapes = pos))
})

test_that("correct_modprob_SE: unknown seqcontext uses 'other' row without error", {
    se  <- make_correction_se(seqconts = rep("UNKNOWN", 15))
    neg <- make_shapes(2, 50)
    pos <- make_shapes(50, 2)
    expect_no_error(correct_modprob_SE(se,
                                       neg_control_shapes = neg,
                                       pos_control_shapes = pos))
})

test_that("correct_modprob_SE: shifts high/low probs in correct direction", {
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

    neg <- make_shapes(2, 50)
    pos <- make_shapes(50, 2)
    res <- correct_modprob_SE(se,
                              neg_control_shapes = neg,
                              pos_control_shapes = pos)

    corr <- as.vector(
        SummarizedExperiment::assay(res, "mod_prob_corrected")[["s1"]][
            cbind(rep(seq_len(n_pos), n_frags),
                  rep(seq_len(n_frags), each = n_pos))])

    expect_gt(mean(corr[raw > 0.8]), mean(raw[raw > 0.8]))
    expect_lt(mean(corr[raw < 0.2]), mean(raw[raw < 0.2]))
})

test_that("correct_modprob_SE: qnorm_to_raw maps to raw value set", {
    se  <- make_correction_se(seed = 7)
    neg <- make_shapes(2, 50)
    pos <- make_shapes(50, 2)
    res <- correct_modprob_SE(se,
                              neg_control_shapes = neg,
                              pos_control_shapes = pos,
                              qnorm_to_raw       = TRUE)

    idx <- SparseArray::nnawhich(
        SummarizedExperiment::assay(se, "mod_prob")[["sampleA"]],
        arr.ind = TRUE)
    raw_vals  <- SummarizedExperiment::assay(se,  "mod_prob")[["sampleA"]][idx]
    corr_vals <- SummarizedExperiment::assay(res, "mod_prob_corrected")[["sampleA"]][idx]
    # rank-lookup maps each corrected value to one of the original raw values
    expect_true(all(corr_vals %in% raw_vals))
})

test_that("correct_modprob_SE: handles multiple samples independently", {
    set.seed(99)
    n_pos <- 10; n_frags <- 5
    make_arr <- function() {
        arr <- SparseArray::NaArray(dim = c(n_pos, n_frags), type = "double")
        arr[cbind(rep(seq_len(n_pos), n_frags),
                  rep(seq_len(n_frags), each = n_pos))] <- runif(n_pos * n_frags)
        arr
    }
    adf <- S4Vectors::make_zero_col_DFrame(nrow = n_pos)
    adf[["s1"]] <- make_arr()
    adf[["s2"]] <- make_arr()
    se <- SummarizedExperiment::SummarizedExperiment(
        assays  = list(mod_prob = adf),
        rowData = S4Vectors::DataFrame(sequenceContext = rep("GCA", n_pos)))

    neg <- make_shapes(2, 50)
    pos <- make_shapes(50, 2)
    res <- correct_modprob_SE(se,
                              neg_control_shapes = neg,
                              pos_control_shapes = pos)

    corr_s1 <- SummarizedExperiment::assay(res, "mod_prob_corrected")[["s1"]]
    corr_s2 <- SummarizedExperiment::assay(res, "mod_prob_corrected")[["s2"]]
    expect_false(identical(corr_s1, corr_s2))
    expect_equal(ncol(res), 2L)
})
