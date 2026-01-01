#' Footprint spectral analysis for a SummarizedExperiment containing
#' single-molecule footprinting (SMF) data
#'
#' @description
#' Performs footprint spectral analysis on a \code{SummarizedExperiment}
#' object containing single-molecule footprinting (SMF) data. For each sample,
#' the function extracts pair-state statistics, computes co-occurrence tables,
#' and infers footprint spectra and emission probabilities using the specified
#' Bayesian inference method.
#'
#' @inheritParams predict_footprints_SE
#' @inheritParams get_ctable_from_SE
#' @inheritParams infer_footprints_vb
#'
#' @return
#' A \code{DataFrame} (from \code{colData(se)}) augmented with additional
#' columns containing:
#' \itemize{
#'   \item pair-state summary statistics,
#'   \item inferred footprint spectra,
#'   \item estimated emission probabilities.
#' }
#' Each row corresponds to a sample in the input \code{SummarizedExperiment}.
#'
#' @importFrom SummarizedExperiment SummarizedExperiment rowRanges colData
#'     colData<-
#' @importFrom SparseArray NaArray
#' @importFrom GenomicRanges GPos match seqnames start end strand seqinfo
#' @importFrom IRanges subsetByOverlaps IRanges IRangesList
#' @importFrom S4Vectors DataFrame SimpleList metadata metadata<-
#'     make_zero_col_DFrame
#' @import data.table
#'
#' @export

ftp_spectral_analysis_SE <- function(
        se,
        assayName = "mod_prob",
        threshMod = 0.5,
        threshUnmod = threshMod,
        min_frag_data_len = 50L,
        min_frag_data_dens = 0.05,
        max_spacing = 200,
        ftp_lengths = 20:200,
        ftp_prior_cover = NULL,
        bg_prior_cover = 0.5,
        total_cnt_prior_dirich = NULL,
        ftp_bg_model = c("informative_prior", "bg_fixed", "ftp_bg_fixed"),
        bg_model_params = list(
            bg_protect_prob_fixed = 0.05,
            bg_protect_min = 0.01,
            bg_protect_max = 0.2,
            bg_protect_mean = 0.05,
            bg_protect_totcount = 100
        ),
        ftp_model_params = list(
            ftp_protect_prob_fixed = 0.95,
            ftp_protect_min = 0.8,
            ftp_protect_max = 0.99,
            ftp_protect_mean = 0.95,
            ftp_protect_totcount = 100
        ),
        max_nruns = 3,
        max_pareto_k = 10,
        ncpu = 1L,
        verbose = FALSE,
        iter = 15000,
        tol_rel_obj = 1e-8,
        output_samples = 2000,
        grad_samples = 1,
        algorithm = "meanfield",
        ...) {
    ftp_bg_model <- match.arg(ftp_bg_model)

    ### validate ncpu
    assert_int(x = ncpu, lower = 0, na.ok = TRUE)
    avail_ncpu <- parallel::detectCores()
    if (is.na(avail_ncpu)) {
        .warning_timestamp(
            "Could not detect number of available cpu. Setting ncpu to 1L."
        )
        ncpu <- 1L
    } else if (ncpu > avail_ncpu || ncpu == 0) {
        .warning_timestamp(c(
            "Number of ncpu is 0 or exceeds number of ",
            "available cpu. Setting ncpu to number of ",
            "available cpus."
        ))
        ncpu <- avail_ncpu
    }

    ## get pair statistics for each sample
    .message_timestamp("Collecting statistics of pair states")

    infDFout <- colData(se)
    ctables_list <- get_ctable_from_SE(se = se,
                                       assayName = assayName,
                                       threshUnmod = threshUnmod,
                                       threshMod = threshMod,
                                       min_frag_data_len = min_frag_data_len,
                                       min_frag_data_dens = min_frag_data_dens,
                                       max_spacing = max_spacing,
                                       aggrSamples = FALSE,
                                       ncpu = ncpu,
                                       verbose = verbose
    )
    infDFout[["pairStats"]] <- ctables_list

    ## perform parameter inference for each sample
    if (verbose) {
        .message_timestamp("Performing footprint spectral analysis")
    }

    infDF <- do.call(rbind, lapply(
        seq_len(nrow(infDFout)),
        function(idx) {
            .message_timestamp(paste0("Inference for ", infDFout$sample[idx]))
            vb_res <- infer_footprints_vb(ctables_list[[idx]],
                                          ftp_lengths,
                                          ftp_prior_cover,
                                          bg_prior_cover,
                                          total_cnt_prior_dirich,
                                          ftp_bg_model,
                                          bg_model_params,
                                          ftp_model_params,
                                          max_nruns,
                                          max_pareto_k,
                                          output_samples = output_samples,
                                          iter = iter,
                                          grad_samples = grad_samples,
                                          tol_rel_obj = tol_rel_obj,
                                          algorithm = algorithm,
                                          refresh = ifelse(verbose, 100, 0),
                                          ...
            )

            DFout <- DataFrame(
                VB_success = FALSE,
                pareto_k = NA,
                bg_emis_mean = NA,
                bg_emis_sd = NA,
                bg_emis_2.5perc = NA,
                bg_emis_50perc = NA,
                bg_emis_97.5perc = NA,
                ftp_emis_mean = NA,
                ftp_emis_sd = NA,
                ftp_emis_2.5perc = NA,
                ftp_emis_50perc = NA,
                ftp_emis_97.5perc = NA,
                bg_coverage_mean = NA,
                bg_coverage_sd = NA,
                bg_coverage_2.5perc = NA,
                bg_coverage_50perc = NA,
                bg_coverage_97.5perc = NA,
                ftp_spectrum = I(list(NULL))
            )

            ## get summary and populate the DataFrame
            if (!is.null(vb_res)) {
                ftpsumm <- get_ftp_inference_summary(vb_res)
                ftp_spec <- ftpsumm$ESTIMATES$ftp_abundance_estimates[, c(
                    "ftp_length", "mean", "sd", "2.5%", "50%", "97.5%")]
                DFout$VB_success <- TRUE
                DFout$pareto_k <- vb_res@sim$diagnostics$psis$pareto_k

                ## add inf results for BG emission probs
                if (ftp_bg_model == "informative_prior") {
                    DFout$bg_emis_mean <-
                        ftpsumm$ESTIMATES$bg_protect_prob_estimate[1, "mean"]
                    DFout$bg_emis_sd <-
                        ftpsumm$ESTIMATES$bg_protect_prob_estimate[1, "sd"]
                    DFout$bg_emis_2.5perc <-
                        ftpsumm$ESTIMATES$bg_protect_prob_estimate[1, "2.5%"]
                    DFout$bg_emis_50perc <-
                        ftpsumm$ESTIMATES$bg_protect_prob_estimate[1, "50%"]
                    DFout$bg_emis_97.5perc <-
                        ftpsumm$ESTIMATES$bg_protect_prob_estimate[1, "97.5%"]
                } else {
                    DFout$bg_emis_mean <-
                        bg_model_params[["bg_protect_prob_fixed"]]
                }

                ## add inf results for FTP emission probs
                if (ftp_bg_model %in% c("informative_prior","bg_fixed")) {
                    DFout$ftp_emis_mean <-
                        ftpsumm$ESTIMATES$ftp_protect_prob_estimate[1, "mean"]
                    DFout$ftp_emis_sd <-
                        ftpsumm$ESTIMATES$ftp_protect_prob_estimate[1, "sd"]
                    DFout$ftp_emis_2.5perc <-
                        ftpsumm$ESTIMATES$ftp_protect_prob_estimate[1, "2.5%"]
                    DFout$ftp_emis_50perc <-
                        ftpsumm$ESTIMATES$ftp_protect_prob_estimate[1, "50%"]
                    DFout$ftp_emis_97.5perc <-
                        ftpsumm$ESTIMATES$ftp_protect_prob_estimate[1, "97.5%"]
                } else {
                    DFout$ftp_emis_mean <-
                        ftp_model_params[["ftp_protect_prob_fixed"]]
                }

                ## add inf results for background
                DFout$bg_coverage_mean <-
                    subset(ftp_spec, ftp_length == 1)[, "mean"]
                DFout$bg_coverage_sd <-
                    subset(ftp_spec, ftp_length == 1)[, "sd"]
                DFout$bg_coverage_2.5perc <-
                    subset(ftp_spec, ftp_length == 1)[, "2.5%"]
                DFout$bg_coverage_50perc <-
                    subset(ftp_spec, ftp_length == 1)[, "50%"]
                DFout$bg_coverage_97.5perc <-
                    subset(ftp_spec, ftp_length == 1)[, "97.5%"]

                DFout$ftp_spectrum[[1]] <- subset(ftp_spec, ftp_length > 1)
            }

            return(DFout)
        }
    ))

    infDFout <- cbind(
        infDFout,
        infDF
    )
    return(infDFout)
}
