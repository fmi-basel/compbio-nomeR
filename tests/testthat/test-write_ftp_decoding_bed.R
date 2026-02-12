test_that("write_ftp_decoding_bed works", {
    ## load data
    dlist <- readRDS(test_path("testdata/test-predict_footprints_SE_data.rds"))


    ftp_models <- list("Nucl" = list("PROTECT_PROB" = rep(0.99, 120),
                                     "COVER_PRIOR" = 0.6,
                                     "NAME" = "Nucl",
                                     "GROUP" = "Nucl"),

                       "TF" = list("PROTECT_PROB" = rep(0.99, 30),
                                   "COVER_PRIOR" = 0.01,
                                   "NAME" = "TF",
                                   "GROUP" = "TF"))

    ftp_pred <- predict_footprints_SE(se = dlist$test_se,
                                      footprint_models = ftp_models,
                                      bgprotectprob = 0.01,
                                      bgcoverprior = 0.59,
                                      ftpConfigMethod = "PV",
                                      returnAs = "data.table",
                                      ncpu = 1,
                                      profile=F)
    nuclconf <- ftp_pred$FOOTPRINT_CONF[ftp_group == "Nucl"]

    tmpdir <- tempdir()
    outfile <- file.path(tmpdir, "nucl_ftp.bed")

    write_ftp_decoding_bed(ftp_decode_dt = nuclconf,
                           file = outfile)
    expect_true(file.exists(outfile))
})
