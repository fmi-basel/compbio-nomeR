test_that("init_ftp_model_params_for_hmc works", {
  expect_error(init_ftp_model_params_for_hmc(-c(200,rep(1,199))))
  expect_type(init_ftp_model_params_for_hmc(),"list")
  
})
