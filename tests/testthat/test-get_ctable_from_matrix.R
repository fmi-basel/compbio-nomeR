test_that("get_ctable_from_matrix works", {
	dlist <- readRDS(test_path("testdata/test-predict_footprints_data.rds"))
	ctbl <- get_ctable_from_matrix(dlist$test_dat_mat)
	expect_equal(ctbl,dlist$exp_ctable)

})
