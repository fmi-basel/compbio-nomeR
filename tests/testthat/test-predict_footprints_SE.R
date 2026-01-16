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
    ## check Viterbi

    ftp_pred <- predict_footprints_SE(se = dlist$test_se,
                                      footprint_models = ftp_models,
                                      bgprotectprob = 0.01,
                                      bgcoverprior = 0.59,
                                      ftpConfigMethod = "Viterbi",
                                      ncpu = 1,
                                      profile=T)
    ## remove timing metadata
    mtdat <- metadata(ftp_pred)
    mtdat$timings <- NULL
    metadata(ftp_pred) <- mtdat

    expect_equal(ftp_pred, dlist$exp_output_viterbi)

    ## check Posterior-Viterbi
    ftp_pred <- predict_footprints_SE(se = dlist$test_se,
                                      footprint_models = ftp_models,
                                      bgprotectprob = 0.01,
                                      bgcoverprior = 0.59,
                                      ftpConfigMethod = "PV",
                                      ncpu = 1)
    ## remove timing metadata
    mtdat <- metadata(ftp_pred)
    mtdat$timings <- NULL
    metadata(ftp_pred) <- mtdat
    expect_equal(ftp_pred, dlist$exp_output_PV)


})
