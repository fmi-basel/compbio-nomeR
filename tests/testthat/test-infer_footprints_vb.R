test_that("infer_footprints_stan_vb works", {

	## simple data with two ftp of 5 and 10 bps
	ftp_5_10_data <- data.frame("S" = 1:15,
															"N00" = c(1626964,1508381,
																				1420066,1336897,
																				1258679,1185045,
																				1136784,1090029,
																				1045066,1001670,
																				959927,983020,
																				1000410,1012326,
																				1019633),
															"N01" = c(0,113856,
																				197657,276539,
																				350664,420399,
																				464873,507988,
																				549466,589487,
																				627997,601629,
																				580965,565776,
																				555217),
															"N10" = c(0, 113921,
																				197854,276862,
																				351147,421042,
																				465716,509005,
																				550645,590837,
																				629533,603328,
																				582814,567737,
																				557220),
															"N11" = c(873036,758842,
																				674423,594702,
																				519510,448514,
																				402627,357978,
																				314823,273006,
																				232543,257023,
																				275811,289161,
																				297930))


	expect_error(sampling_res <- infer_footprints_vb(ftp_5_10_data,
																									 ftp_lengths = -1:20,
																									 iter=1,
																									 adapt_iter = 1))
	expect_error(sampling_res <- infer_footprints_vb(ftp_5_10_data[-1,],
																												ftp_lengths = -1:20,
																												iter=1,
																												adapt_iter = 1))
	expect_error(sampling_res <- infer_footprints_vb(ftp_5_10_data[,1:3],
																									 infer_footprints_vb))

	## run 1 iteration
	vb_res <- infer_footprints_vb(cooc_ctable = ftp_5_10_data,
											ftp_lengths = 2:15,
											iter = 1,
											adapt_iter = 1)

	## check if return is stanfit object
	expect_s4_class(vb_res,class = "stanfit")

})
