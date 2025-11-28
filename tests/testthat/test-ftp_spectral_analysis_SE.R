test_that("ftp_spectral_analysis_SE works", {
    ## load data
    dlist <- readRDS(test_path("testdata/test-predict_footprints_SE_data.rds"))
    expect_warning(tst_fsa <- ftp_spectral_analysis_SE(dlist$test_se,
    																									 iter = 1,
    																									 adapt_iter = 1,
    																									 verbose = FALSE,
    																									 max_nruns = 1,
    																									 max_pareto_k = Inf,
    																									 output_samples = 1))


    ## test that all columns are present
    expect_true(all(c('sample',
                      'modbase',
                      'n_reads',
                      'readInfo',
                      'pairStats',
                      'VB_success',
                      'pareto_k',
                      'bg_emis_mean',
                      'bg_emis_sd',
                      'bg_emis_2.5perc',
                      'bg_emis_50perc',
                      'bg_emis_97.5perc',
                      'ftp_emis_mean',
                      'ftp_emis_sd',
                      'ftp_emis_2.5perc',
                      'ftp_emis_50perc',
                      'ftp_emis_97.5perc',
                      'bg_coverage_mean',
                      'bg_coverage_sd',
                      'bg_coverage_2.5perc',
                      'bg_coverage_50perc',
                      'bg_coverage_97.5perc',
                      'ftp_spectrum') %in%
                        colnames(tst_fsa)))
})
