## shared SE test fixtures ----------------------------------------------------
.se_ftp_models <- function() {
    list("Nucl" = list("PROTECT_PROB" = rep(0.99, 120),
                       "COVER_PRIOR"  = 0.6,
                       "NAME"         = "Nucl",
                       "GROUP"        = "Nucl"),
         "TF"   = list("PROTECT_PROB" = rep(0.99, 30),
                       "COVER_PRIOR"  = 0.01,
                       "NAME"         = "TF",
                       "GROUP"        = "TF"))
}

test_that("predict_footprints_SE works", {
    ## load data
    dlist <- readRDS(test_path("testdata/test-predict_footprints_SE_data.rds"))

    ### test with SE
    ftp_models <- list("Nucl" = list("PROTECT_PROB" = rep(0.99, 120),
                                     "COVER_PRIOR" = 0.6,
                                     "NAME" = "Nucl",
                                     "GROUP" = "Nucl"),

                       "TF" = list("PROTECT_PROB" = rep(0.99, 30),
                                   "COVER_PRIOR" = 0.01,
                                   "NAME" = "TF",
                                   "GROUP" = "TF"))
    ## test Viterbi

    ftp_pred <- predict_footprints_SE(se = dlist$test_se,
                                      footprint_models = ftp_models,
                                      bgprotectprob = 0.01,
                                      bgcoverprior = 0.59,
                                      ftpConfigMethod = "Viterbi",
                                      keepStartProb=TRUE,
                                      ncpu = 1,
                                      profile=T)
    ## remove timing metadata
    mtdat <- metadata(ftp_pred)
    mtdat$timings <- NULL
    metadata(ftp_pred) <- mtdat

    expect_equal(ftp_pred, dlist$exp_output_viterbi)

    ## test Posterior-Viterbi
    ftp_pred <- predict_footprints_SE(se = dlist$test_se,
                                      footprint_models = ftp_models,
                                      bgprotectprob = 0.01,
                                      bgcoverprior = 0.59,
                                      ftpConfigMethod = "PV",
                                      keepStartProb=TRUE,
                                      ncpu = 1)
    ## remove timing metadata
    mtdat <- metadata(ftp_pred)
    mtdat$timings <- NULL
    metadata(ftp_pred) <- mtdat
    expect_equal(ftp_pred, dlist$exp_output_PV)

    ## test returnAs="data.table"
    ftp_pred <- predict_footprints_SE(se = dlist$test_se,
                                      footprint_models = ftp_models,
                                      bgprotectprob = 0.01,
                                      bgcoverprior = 0.59,
                                      ftpConfigMethod = "PV",
                                      returnAs = "data.table",
                                      keepStartProb=TRUE,
                                      ncpu = 1)
    expect_equal(ftp_pred, dlist$exp_output_PV_dt)
})


test_that("predict_footprints_SE PosteriorDecoding method returns valid output", {
    dlist <- readRDS(test_path("testdata/test-predict_footprints_SE_data.rds"))
    ftp_models <- .se_ftp_models()

    out <- predict_footprints_SE(se               = dlist$test_se,
                                 footprint_models  = ftp_models,
                                 bgprotectprob     = 0.01,
                                 bgcoverprior      = 0.59,
                                 ftpConfigMethod   = "PosteriorDecoding",
                                 ncpu              = 1)

    ## same assay structure as other methods
    expect_true(all(c("Nucl_coverProb_nomeR", "TF_coverProb_nomeR",
                      "background_coverProb_nomeR") %in% assayNames(out)))

    ## footprint configs present in colData
    expect_true(all(c("Nucl_nomeR", "TF_nomeR") %in% names(colData(out))))

    ## scores in [0,1]
    all_scores <- unlist(lapply(colData(out)$Nucl_nomeR$s1,
                                function(ir) S4Vectors::mcols(ir)$score))
    expect_true(all(all_scores >= 0 - .Machine$double.eps^0.5 &
                        all_scores <= 1 + .Machine$double.eps^0.5))
})


test_that("predict_footprints_SE ncpu=1 and ncpu=2 produce identical results", {
    dlist <- readRDS(test_path("testdata/test-predict_footprints_SE_data.rds"))
    ftp_models <- .se_ftp_models()

    run <- function(ncpu)
        predict_footprints_SE(se               = dlist$test_se,
                              footprint_models  = ftp_models,
                              bgprotectprob     = 0.01,
                              bgcoverprior      = 0.59,
                              ftpConfigMethod   = "PV",
                              keepStartProb     = TRUE,
                              ncpu              = ncpu)

    out1 <- run(1)
    out2 <- run(2)

    ## strip timings before comparison
    strip_timings <- function(se) {
        m <- metadata(se); m$timings <- NULL; metadata(se) <- m; se
    }
    expect_equal(strip_timings(out1), strip_timings(out2))
})


test_that("predict_footprints_SE keepStartProb=FALSE omits startProb assays", {
    dlist <- readRDS(test_path("testdata/test-predict_footprints_SE_data.rds"))
    ftp_models <- .se_ftp_models()

    out <- predict_footprints_SE(se               = dlist$test_se,
                                 footprint_models  = ftp_models,
                                 bgprotectprob     = 0.01,
                                 bgcoverprior      = 0.59,
                                 keepStartProb     = FALSE,
                                 ncpu              = 1)

    expect_false(any(grepl("_startProb_nomeR", assayNames(out))))
    expect_true(all(c("Nucl_coverProb_nomeR",
                      "TF_coverProb_nomeR",
                      "background_coverProb_nomeR") %in% assayNames(out)))
})


test_that("predict_footprints_SE aggrByGroup=FALSE reports per-name assays and colData", {
    dlist <- readRDS(test_path("testdata/test-predict_footprints_SE_data.rds"))
    ftp_models <- .se_ftp_models()

    out <- predict_footprints_SE(se               = dlist$test_se,
                                 footprint_models  = ftp_models,
                                 bgprotectprob     = 0.01,
                                 bgcoverprior      = 0.59,
                                 aggrByGroup       = FALSE,
                                 keepStartProb     = FALSE,
                                 ncpu              = 1)

    ## assay names use model NAMEs ("Nucl", "TF"), same here since NAME==GROUP
    expect_true(all(c("Nucl_coverProb_nomeR", "TF_coverProb_nomeR") %in%
                        assayNames(out)))

    ## colData columns named by ftp_name
    expect_true(all(c("Nucl_nomeR", "TF_nomeR") %in% names(colData(out))))
})


test_that("predict_footprints_SE cover probabilities sum to 1 at each non-NA position", {
    dlist <- readRDS(test_path("testdata/test-predict_footprints_SE_data.rds"))
    ftp_models <- .se_ftp_models()

    out <- predict_footprints_SE(se               = dlist$test_se,
                                 footprint_models  = ftp_models,
                                 bgprotectprob     = 0.01,
                                 bgcoverprior      = 0.59,
                                 ncpu              = 1)

    cover_assays <- grep("_coverProb_nomeR", assayNames(out), value = TRUE)

    ## for each sample column, sum the cover-prob assays at every row position
    for (sI in seq_len(ncol(out))) {
        cover_mat <- do.call(cbind, lapply(cover_assays, function(a) {
            as.numeric(assay(out, a)[[sI]])
        }))
        row_sums <- rowSums(cover_mat)
        non_na_rows <- which(!is.na(row_sums))
        expect_true(all(abs(row_sums[non_na_rows] - 1) < 1e-8))
    }
})


test_that("predict_footprints_SE errors when all fragments are filtered out", {
    dlist <- readRDS(test_path("testdata/test-predict_footprints_SE_data.rds"))
    ftp_models <- .se_ftp_models()

    expect_error(
        predict_footprints_SE(se                = dlist$test_se,
                              footprint_models   = ftp_models,
                              bgprotectprob      = 0.01,
                              bgcoverprior       = 0.59,
                              min_frag_data_len  = 100000L,
                              ncpu               = 1),
        regexp = "No fragments"
    )
})


test_that("predict_footprints_SE verbose=TRUE runs without error", {
    dlist <- readRDS(test_path("testdata/test-predict_footprints_SE_data.rds"))
    ftp_models <- .se_ftp_models()

    expect_no_error(
        predict_footprints_SE(se               = dlist$test_se,
                              footprint_models  = ftp_models,
                              bgprotectprob     = 0.01,
                              bgcoverprior      = 0.59,
                              verbose           = TRUE,
                              ncpu              = 1)
    )
})

