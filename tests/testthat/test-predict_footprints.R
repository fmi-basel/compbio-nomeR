test_that("wrong parameters for predict_footprints are handled correctly",{
  expect_error(predict_footprints(data = "aaa"))
  expect_error(predict_footprints(data = matrix()))
  
  ## check whether ncpu handled correctly
  ## create dummy random data
  set.seed(3346)
  nc <- 50
  nr <- 50
  rmatr <- matrix(data = as.integer(rnorm(nc * nr) >= 0.5),
                  ncol = nc,nrow=nr)
  
  ## create dummy footprints
  bg.pr <- 0.5
  ft.pr <- 1-bg.pr
  ft.len <- 15
  
  ## creating a list of binding models for nomeR
  ftp.models <- list(list("PROTECT_PROB" = rep(0.99,ft.len),
                          "COVER_PRIOR" = ft.pr,
                          "NAME" = "FOOTPRINT"))
  
  
  # 
  # expect_warning(nomeR.out <- predict_footprints(data=rmatr,
  #                                           footprint_models = ftp.models,
  #                                           bgprotectprob = 0.05,
  #                                           bgcoverprior = bg.pr,
  #                                           ncpu = 0))
  # expect_warning(nomeR.out <- predict_footprints(data=rmatr,
  #                                           footprint_models = ftp.models,
  #                                           bgprotectprob = 0.05,
  #                                           bgcoverprior = bg.pr,
  #                                           ncpu = Inf))
  
  expect_error(nomeR.out <- predict_footprints(data=rmatr,
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
                  ncol = nc,nrow=nr)
  
  ## create dummy footprints
  bg.pr <- 0.5
  ft.pr <- 1-bg.pr
  ft.len <- 15
  
  ## creating a list of binding models for nomeR
  ftp.models <- list(list("PROTECT_PROB" = rep(0.99,ft.len),
                          "COVER_PRIOR" = ft.pr,
                          "NAME" = "FOOTPRINT"))
  
  nomeR.out <- predict_footprints(data=rmatr,
                                  footprint_models = ftp.models,
                                  bgprotectprob = 0.05,
                                  bgcoverprior = bg.pr,
                                  report_prediction_in_flanks = T,
                                  ncpu = 1L)
  
  ## check whether slots exist
  expect_true(all(c("START_PROB", "COVER_PROB","SUMMARY") %in% names(nomeR.out)))
  
  ## check whether all required seq exist
  expect_true(all(as.character(1:nr) %in% nomeR.out[["START_PROB"]][["seq"]]) &
                all(as.character(1:nr) %in% nomeR.out[["COVER_PROB"]][["seq"]])
  )
  
  ## check all pos exist
  expect_true(all(1:nc %in% nomeR.out[["START_PROB"]][["pos"]]) &
                all(1:nc %in% nomeR.out[["COVER_PROB"]][["pos"]])
  )
  
  ## check whether FOOTPRINT and background exist
  expect_true(all(c("FOOTPRINT","background") %in% colnames(nomeR.out[["START_PROB"]])) &
                all(c("FOOTPRINT","background") %in% colnames(nomeR.out[["COVER_PROB"]])) &
                all(c("FOOTPRINT","background") %in% colnames(nomeR.out[["SUMMARY"]]))
  )
  
  ## check that we do not have incorrect probs
  expect_false(any(nomeR.out[["START_PROB"]][,c("FOOTPRINT","background")] < 0 - .Machine$double.eps^0.5) | any(nomeR.out[["START_PROB"]][,c("FOOTPRINT","background")] > 1 + .Machine$double.eps^0.5) |
                 any(nomeR.out[["COVER_PROB"]][,c("FOOTPRINT","background")] < 0 - .Machine$double.eps^0.5) | any(nomeR.out[["COVER_PROB"]][,c("FOOTPRINT","background")] > 1 + .Machine$double.eps^0.5))
  
  ## check whether sum of start probs does not exceed 1
  expect_false(any(rowSums(nomeR.out[["START_PROB"]][,c("FOOTPRINT","background")]) > 1 + .Machine$double.eps^0.5))
  
  ## check if sum of cover probs sum up to 1
  cover.prob.rowsum <- rowSums(nomeR.out[["COVER_PROB"]][,c("FOOTPRINT","background")])
  expect_true(all(abs(cover.prob.rowsum - 1) < 1.0e-8))
  # a <- sapply(cover.prob.rowsum,
  #             function(x){
  #               
  #               expect_equal(x,1,tolerance = 1.0e-8)
  #             })
  
  
  
})



