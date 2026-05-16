test_that("get_ctable_from_SE works", {
	## load data
	dlist <- readRDS(test_path("testdata/test-predict_footprints_SE_data.rds"))
	ctbl <- get_ctable_from_SE(dlist$test_se)
	expect_equal(ctbl,dlist$exp_ctable)

})
