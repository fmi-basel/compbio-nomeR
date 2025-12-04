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


test_that("predict_footprints returns expected probabilities and MAP configuration",{

    ## load data
    dlist <- readRDS(test_path("testdata/test-predict_footprints_data.rds"))

    ## calculate for all footprints aggregated by group
    testinsil <- predict_footprints(data=dlist$test_dat_mat,
                                    footprint_models = dlist$ftp_models,
                                    bgprotectprob = 0.05304034,
                                    bgcoverprior = 0.4822005,
                                    aggrByGroup = TRUE,
                                    ftpConfigMethod = "Viterbi",
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
                                    ncpu = 1L)
    map_conf <- subset(testinsil$FOOTPRINT_CONF,ftp_name != "background")
    map_conf <- map_conf[order(map_conf$start),]
    row.names(map_conf) <- NULL

    expect_equal(map_conf[,colnames(exp_conf)],exp_conf)

})


