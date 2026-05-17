## shared test fixtures -------------------------------------------------------
.make_rmatr <- function(nr = 50, nc = 50, seed = 3346) {
    set.seed(seed)
    matrix(data = as.integer(rnorm(nc * nr) >= 0.5), ncol = nc, nrow = nr)
}
.make_single_ftp_models <- function(ft.len = 15, ft.pr = 0.5) {
    list(list("PROTECT_PROB" = rep(0.99, ft.len),
              "COVER_PRIOR"  = ft.pr,
              "NAME"         = "FOOTPRINT"))
}

test_that("wrong parameters for predict_footprints are handled correctly",{
    expect_error(predict_footprints(data = "aaa"))
    expect_error(predict_footprints(data = matrix()))

    ## check whether ncpu handled correctly
    ## create dummy random data
    set.seed(3346)
    nc <- 50
    nr <- 50
    rmatr <- matrix(data = as.integer(rnorm(nc * nr) >= 0.5),
                    ncol = nc, nrow = nr)

    ## create dummy footprints
    bg.pr <- 0.5
    ft.pr <- 1 - bg.pr
    ft.len <- 15

    ## creating a list of binding models for nomeR
    ftp.models <- list(list("PROTECT_PROB" = rep(0.99, ft.len),
                            "COVER_PRIOR" = ft.pr,
                            "NAME" = "FOOTPRINT"))

    expect_error(nomeR.out <- predict_footprints(data = rmatr,
                                                 footprint_models = ftp.models,
                                                 bgprotectprob = 0.05,
                                                 bgcoverprior = bg.pr,
                                                 ncpu = -1))

})


test_that("predict_footprints returns correct object",{
    ## create dummy random data
    set.seed(3346)
    nc <- 50
    nr <- 50
    rmatr <- matrix(data = as.integer(rnorm(nc * nr) >= 0.5),
                    ncol = nc, nrow = nr)

    ## create dummy footprints
    bg.pr <- 0.5
    ft.pr <- 1 - bg.pr
    ft.len <- 15

    ## creating a list of binding models for nomeR
    ftp.models <- list(list("PROTECT_PROB" = rep(0.99, ft.len),
                            "COVER_PRIOR" = ft.pr,
                            "NAME" = "FOOTPRINT"))

    nomeR.out <- predict_footprints(data = rmatr,
                                    footprint_models = ftp.models,
                                    bgprotectprob = 0.05,
                                    bgcoverprior = bg.pr,
                                    keepStartProb=TRUE,
                                    ncpu = 1L)

    ## check whether slots exist
    expect_true(all(c("START_PROB", "COVER_PROB", "FOOTPRINT_CONF") %in% names(nomeR.out)))

    ## check whether all required seq exist
    expect_true(all(as.character(seq_len(nr)) %in% nomeR.out[["START_PROB"]][["seq"]]) &
                    all(as.character(seq_len(nr)) %in% nomeR.out[["COVER_PROB"]][["seq"]])
    )

    ## check all pos exist
    expect_true(all(1:nc %in% nomeR.out[["START_PROB"]][["pos"]]) &
                    all(1:nc %in% nomeR.out[["COVER_PROB"]][["pos"]])
    )

    ## check whether FOOTPRINT and background exist
    expect_true(all(c("FOOTPRINT", "background") %in% colnames(nomeR.out[["START_PROB"]])) &
                    all(c("FOOTPRINT", "background") %in% colnames(nomeR.out[["COVER_PROB"]]))
    )

    ## check that we do not have incorrect probs
    expect_false(any(nomeR.out[["START_PROB"]][, c("FOOTPRINT", "background")] < 0 - .Machine$double.eps ^ 0.5) | any(nomeR.out[["START_PROB"]][, c("FOOTPRINT", "background")] > 1 + .Machine$double.eps ^ 0.5) |
                     any(nomeR.out[["COVER_PROB"]][, c("FOOTPRINT", "background")] < 0 - .Machine$double.eps ^ 0.5) | any(nomeR.out[["COVER_PROB"]][, c("FOOTPRINT", "background")] > 1 + .Machine$double.eps ^ 0.5))

    ## check whether sum of start probs does not exceed 1
    expect_false(any(rowSums(nomeR.out[["START_PROB"]][, c("FOOTPRINT", "background")]) > 1 + .Machine$double.eps ^ 0.5))

    ## check if sum of cover probs sum up to 1
    cover.prob.rowsum <- rowSums(nomeR.out[["COVER_PROB"]][, c("FOOTPRINT", "background")])
    expect_true(all(abs(cover.prob.rowsum - 1) < 1.0e-8))

})


test_that("predict_footprints returns expected probabilities and ftp configuration",{

    ## load data
    dlist <- readRDS(test_path("testdata/test-predict_footprints_data.rds"))


    ## calculate for all footprints aggregated by group
    testinsil <- predict_footprints(data=dlist$test_dat_mat,
                                    footprint_models = dlist$ftp_models,
                                    bgprotectprob = 0.05304034,
                                    bgcoverprior = 0.4822005,
                                    aggrByGroup = TRUE,
                                    ftpConfigMethod = "Viterbi",
                                    keepStartProb=TRUE,
                                    ncpu = 1L)
    ## check start probs

    expect_equal(testinsil$START_PROB,dlist$exp_output$START_PROB)
    ## check cover probs
    expect_equal(testinsil$COVER_PROB,dlist$exp_output$COVER_PROB)

    ## check whether Viterbi config is correct
    map_conf <- subset(testinsil$FOOTPRINT_CONF,ftp_name != "background")
    map_conf <- map_conf[order(map_conf$start),]
    row.names(map_conf) <- NULL
    exp_conf <- data.frame(seq = 1,
                           start=c(151,273,492),
                           width=c(50,150,150),
                           ftp_name = c("ftp1--50","ftp2--150","ftp2--150"),
                           ftp_group = c("ftp1","ftp2","ftp2")
                           #score = c(0.9421603, 0.9496543, 0.8856803)
                           )
    expect_equal(map_conf[,colnames(exp_conf)],exp_conf)


    ## check if Posterior-Viterbi is correct
    testinsil <- predict_footprints(data=dlist$test_dat_mat,
                                    footprint_models = dlist$ftp_models,
                                    bgprotectprob = 0.05304034,
                                    bgcoverprior = 0.4822005,
                                    aggrByGroup = TRUE,
                                    ftpConfigMethod = "PV",
                                    keepStartProb=TRUE,
                                    ncpu = 1L)
    map_conf <- subset(testinsil$FOOTPRINT_CONF,ftp_name != "background")
    map_conf <- map_conf[order(map_conf$start),]
    row.names(map_conf) <- NULL

    expect_equal(map_conf[,colnames(exp_conf)],exp_conf)

})


test_that("predict_footprints PosteriorDecoding method returns valid output", {
    rmatr <- .make_rmatr()
    ftp.models <- .make_single_ftp_models()

    out <- predict_footprints(data            = rmatr,
                              footprint_models = ftp.models,
                              bgprotectprob    = 0.05,
                              bgcoverprior     = 0.5,
                              ftpConfigMethod  = "PosteriorDecoding",
                              ncpu             = 1L)

    expect_named(out, c("COVER_PROB", "FOOTPRINT_CONF"))
    expect_true(all(c("FOOTPRINT", "background") %in% colnames(out$COVER_PROB)))

    conf <- out$FOOTPRINT_CONF
    expect_true(all(c("seq", "start", "width", "ftp_name", "ftp_group", "score") %in%
                        colnames(conf)))
    expect_true(all(conf$score >= 0 - .Machine$double.eps^0.5 &
                        conf$score <= 1 + .Machine$double.eps^0.5))

    ## cover probs still sum to 1 with PD
    cover_sum <- rowSums(out$COVER_PROB[, c("FOOTPRINT", "background")])
    expect_true(all(abs(cover_sum - 1) < 1e-8))
})


test_that("predict_footprints ncpu=1 and ncpu=2 produce identical results", {
    dlist <- readRDS(test_path("testdata/test-predict_footprints_data.rds"))

    run <- function(ncpu)
        predict_footprints(data             = dlist$test_dat_mat,
                           footprint_models  = dlist$ftp_models,
                           bgprotectprob     = 0.05304034,
                           bgcoverprior      = 0.4822005,
                           aggrByGroup       = TRUE,
                           ftpConfigMethod   = "PV",
                           keepStartProb     = TRUE,
                           ncpu              = ncpu)

    out1 <- run(1L)
    out2 <- run(2L)

    expect_equal(out1$START_PROB,     out2$START_PROB)
    expect_equal(out1$COVER_PROB,     out2$COVER_PROB)
    expect_equal(out1$FOOTPRINT_CONF, out2$FOOTPRINT_CONF)
})


test_that("predict_footprints keepStartProb=FALSE omits START_PROB", {
    rmatr <- .make_rmatr()
    ftp.models <- .make_single_ftp_models()

    out <- predict_footprints(data             = rmatr,
                              footprint_models  = ftp.models,
                              bgprotectprob     = 0.05,
                              bgcoverprior      = 0.5,
                              keepStartProb     = FALSE,
                              ncpu              = 1L)

    expect_named(out, c("COVER_PROB", "FOOTPRINT_CONF"))
    expect_false("START_PROB" %in% names(out))
})


test_that("predict_footprints accepts list input and matches matrix input", {
    rmatr <- .make_rmatr()
    ftp.models <- .make_single_ftp_models()

    lst <- lapply(seq_len(nrow(rmatr)), function(i) rmatr[i, ])

    out_mat <- predict_footprints(data             = rmatr,
                                  footprint_models  = ftp.models,
                                  bgprotectprob     = 0.05,
                                  bgcoverprior      = 0.5,
                                  keepStartProb     = TRUE,
                                  ncpu              = 1L)
    out_lst <- predict_footprints(data             = lst,
                                  footprint_models  = ftp.models,
                                  bgprotectprob     = 0.05,
                                  bgcoverprior      = 0.5,
                                  keepStartProb     = TRUE,
                                  ncpu              = 1L)

    expect_equal(out_lst$COVER_PROB,     out_mat$COVER_PROB)
    expect_equal(out_lst$START_PROB,     out_mat$START_PROB)
    expect_equal(out_lst$FOOTPRINT_CONF, out_mat$FOOTPRINT_CONF)
})


test_that("predict_footprints aggrByGroup=FALSE reports per-name columns", {
    dlist <- readRDS(test_path("testdata/test-predict_footprints_data.rds"))
    model_names <- vapply(dlist$ftp_models, `[[`, character(1), "NAME")

    out <- predict_footprints(data             = dlist$test_dat_mat,
                              footprint_models  = dlist$ftp_models,
                              bgprotectprob     = 0.05304034,
                              bgcoverprior      = 0.4822005,
                              aggrByGroup       = FALSE,
                              keepStartProb     = TRUE,
                              ncpu              = 1L)

    ## output columns are footprint NAMEs, not GROUPs
    expect_true(all(model_names %in% colnames(out$COVER_PROB)))
    expect_true(all(model_names %in% colnames(out$START_PROB)))

    ## cover probs still sum to 1
    cover_sum <- rowSums(out$COVER_PROB[, c(model_names, "background")])
    expect_true(all(abs(cover_sum - 1) < 1e-8))
})


test_that("predict_footprints works for a single-fragment input", {
    set.seed(42)
    single_row <- matrix(as.integer(rnorm(60) >= 0.5), nrow = 1)
    ftp.models <- .make_single_ftp_models(ft.len = 15)

    out <- predict_footprints(data             = single_row,
                              footprint_models  = ftp.models,
                              bgprotectprob     = 0.05,
                              bgcoverprior      = 0.5,
                              ncpu              = 1L)

    expect_equal(nrow(out$COVER_PROB), ncol(single_row))
    cover_sum <- rowSums(out$COVER_PROB[, c("FOOTPRINT", "background")])
    expect_true(all(abs(cover_sum - 1) < 1e-8))
})


test_that("predict_footprints verbose=TRUE runs without error", {
    rmatr <- .make_rmatr(nr = 5)
    ftp.models <- .make_single_ftp_models()

    expect_no_error(
        predict_footprints(data             = rmatr,
                           footprint_models  = ftp.models,
                           bgprotectprob     = 0.05,
                           bgcoverprior      = 0.5,
                           verbose           = TRUE,
                           ncpu              = 1L)
    )
})
