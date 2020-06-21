test_that("init_ftp_abundance_for_inference works", {
  expect_error(init_ftp_abundance_for_inference(-c(200,rep(1,199))))
  expect_type(init_ftp_abundance_for_inference(),"list")
  
})
