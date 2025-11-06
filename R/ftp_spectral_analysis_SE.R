#' Footprint spectral analysis for SummarizedExperiment object containing
#' single-molecule footprinting (SMF) data
#'
#'
#' @param se A \code{\link[SummarizedExperiment]{SummarizedExperiment}} object
#'     containing read-level data returned by \code{\link[footprintR]{readModBam}} function
#'     containing modification probabilities.
#' @param assayName Character scalar describing the name of tje assay in \code{se} containing
#'     read-level data.
#' @param threshUnmod,threshMod Numeric scalars used to classify observations
#'     as modified (modification probability >= threshMod), unmodified
#'     (modification probability < threshUnmod) or unknown (otherwise).
#' @param max_spacing \code{integer} that defines maximum distance between positions for aggregating
#'     frequencies of combinations 00, 01, 10, 11 at distances up to \code{max_spacing} observed in
#'     the SMF dataset.
#' @param ftp_lengths A numeric vector representing the lengths of footprints
#'     for which abundance is being analyzed. This parameter allows users to
#'     input a vector of footprint lengths of interest for further analysis,
#'     excluding the length of 1, which is reserved for background.
#' @param ftp_prior_cover A numeric vector representing the expected coverages
#'     for footprints with lengths corresponding to the values provided in the
#'     ftp_lengths parameter. This parameter allows users to specify the
#'     expected coverage for each footprint length, influencing the Dirichlet
#'     distribution used as the prior distribution in the Bayesian model.
#'     Higher values indicate a higher expected coverage for the corresponding
#'     footprint length, affecting the model's prior assumptions. Users can
#'     adjust this parameter to reflect their prior knowledge or assumptions
#'     about the coverage of specific footprint lengths in the dataset. Each
#'     value in the vector must be between 0 and 1. The sum of the values in
#'     ftp_prior_cover and bg_prior_cover should equal 1. If the sum is not 1,
#'     the values are scaled accordingly, and a warning is issued.
#' @param bg_prior_cover A numeric value representing the expected fraction of
#'     unprotected positions, or coverage of background, in the dataset. This
#'     parameter influences the Dirichlet distribution used as the prior
#'     distribution in the Bayesian model. Higher values indicate a higher
#'     proportion of background coverage, affecting the model's prior
#'     assumptions. Users can adjust this parameter to reflect their prior
#'     knowledge or assumptions about the background coverage in the dataset.
#'     It must be a value between 0 and 1.
#' @param total_cnt_prior_dirich A numeric value representing the total count
#'     parameterization of the prior Dirichlet distribution used in the
#'     Bayesian model. This parameter is related to the bg_prior_cover and
#'     ftp_prior_cover parameters and reflects the total count of observations.
#'     In Bayesian statistics, the Dirichlet distribution is often
#'     parameterized using mean and total count, where the total count
#'     influences the spread of the distribution. Higher values of
#'     total_cnt_prior_dirich result in a narrower distribution, implying
#'     stronger prior beliefs, while lower values lead to a wider distribution,
#'     indicating weaker prior beliefs. Users can adjust this parameter to
#'     reflect their confidence in the prior assumptions encoded by
#'     bg_prior_cover and ftp_prior_cover.
#' @param ftp_bg_model Type of model used for inference.
#'     \describe{
#'     \item{"informative_prior"}{Inference is performed on parameters
#'     bg_protect_prob, ftp_protect_prob, and footprint abundances.}
#'     \item{"bg_fixed"}{bg_protect_prob is fixed and determined by
#'     bg_model_params[["bg_protect_prob_fixed"]], while inference is conducted
#'     on ftp_protect_prob and footprint abundances.}
#'     \item{"ftp_bg_fixed"}{Both bg_protect_prob and ftp_protect_prob are
#'     fixed, defined by corresponding values in bg_model_params and
#'     ftp_model_params, respectively. Inference is solely focused on footprint
#'     abundances.}
#' }
#' @param bg_model_params A list containing parameters for the background
#'     model, which must contain the following elements:
#'     \describe{
#'     \item{bg_protect_prob_fixed}{Constant value for the model parameter
#'     "bg_protect_prob" used in "bg_fixed" and "ftp_bg_fixed" models and
#'     ignored if ftp_bg_model is "informative_prior".}
#'     \item{bg_protect_min}{Minimum allowed value for the background emission
#'     probability of the protected state (i.e., 1's in the data). This value
#'     is ignored if ftp_bg_model is "bg_fixed" or "ftp_bg_fixed".}
#'     \item{bg_protect_max}{Maximum allowed value for the background emission
#'     probability of the protected state (i.e., 1's in the data). This value
#'     is ignored if ftp_bg_model is "bg_fixed" or "ftp_bg_fixed".}
#'     \item{bg_protect_mean}{Mean of the prior Beta distribution for the model
#'     parameter "bg_protect_prob". This parameter corresponds to shape
#'     parameters of the Beta distribution as mean=alpha/(alpha+beta). This
#'     value is ignored if ftp_bg_model is "bg_fixed" or "ftp_bg_fixed".}
#'     \item{bg_protect_totcount}{Total count parameter for the prior beta
#'     distribution for the model parameter "bg_protect_prob". This parameter
#'     corresponds to shape parameters of the Beta distribution as
#'     tot_count=alpha+beta and influences the spread of the distribution.
#'     This value is ignored if ftp_bg_model is "bg_fixed" or "ftp_bg_fixed".}
#' }
#' @param ftp_model_params A list containing parameters for the footprint
#'     model, which must contain the following elements:
#'     \describe{
#'     \item{ftp_protect_prob_fixed}{Constant value for the model parameter
#'     "ftp_protect_prob" used in "ftp_bg_fixed" models and ignored if
#'     ftp_bg_model is "informative_prior" or "bg_fixed".}
#'     \item{ftp_protect_min}{Minimum allowed value for the footprint emission
#'     probability of the protected state (i.e., 1's in the data). This value
#'     is ignored if ftp_bg_model is "ftp_bg_fixed".}
#'     \item{ftp_protect_max}{Maximum allowed value for the footprint emission
#'     probability of the protected state (i.e., 1's in the data). This value
#'     is ignored if ftp_bg_model is "ftp_bg_fixed".}
#'     \item{ftp_protect_mean}{Alpha parameter for the prior beta distribution
#'     for the model parameter "ftp_protect_prob". This parameter corresponds
#'     to shape parameters of the Beta distribution as mean=alpha/(alpha+beta).
#'     This value is ignored if ftp_bg_model is "ftp_bg_fixed".}
#'     \item{ftp_protect_totcount}{Beta parameter for the prior beta
#'     distribution for the model parameter "ftp_protect_prob". This parameter
#'     corresponds to shape parameters of the Beta distribution as
#'     tot_count=alpha+beta and influences the spread of the distribution.
#'     This value is ignored if ftp_bg_model is "ftp_bg_fixed".}
#' }
#' @param max_nruns Maximum number of trials to run stan function
#'     \code{\link[rstan]{vb}}. Sometimes, due to bad initial point or other
#'     reasons this function fails to converge. \code{max_nruns} controls
#'     maximum number of attempts for inference.
#' @param max_pareto_k maximum pareto_k returned by \code{\link[rstan]{vb}}.
#'     If it exceeds \code{max_pareto_k} the function will run again until
#'     max_nruns attempts have been done.
#' @param ncpu number of threads to use.
#' @param verbose verbose mode for bug fixing.
#'
#' @return \code{DataFrame} object from \code{colData} of \code{se} with additional columns
#'     containing pair state statistics, inferred footprint spectra and emission probabilities
#'     for each sample in code{se}.
#' @importFrom SummarizedExperiment SummarizedExperiment rowRanges colData colData<-
#' @importFrom SparseArray NaArray
#' @importFrom GenomicRanges GPos match seqnames start end strand seqinfo
#' @importFrom IRanges subsetByOverlaps IRanges IRangesList
#' @importFrom S4Vectors DataFrame SimpleList metadata metadata<- make_zero_col_DFrame
#' @import data.table
#'
#' @export
ftp_spectral_analysis_SE <- function(se,
                                     assayName = "mod_prob",
                                     threshUnmod = 0.5,
                                     threshMod = 0.5,
                                     max_spacing = 200,
                                     ftp_lengths = 20:200,
                                     ftp_prior_cover = NULL,
                                     bg_prior_cover = 0.5,
                                     total_cnt_prior_dirich = NULL,
                                     ftp_bg_model = c("informative_prior", "bg_fixed", "ftp_bg_fixed"),
                                     bg_model_params = list(
                                       bg_protect_prob_fixed = 0.05, bg_protect_min = 0.01,
                                       bg_protect_max = 0.2, bg_protect_mean = 0.05, bg_protect_totcount = 100
                                     ),
                                     ftp_model_params = list(
                                       ftp_protect_prob_fixed = 0.95, ftp_protect_min = 0.8,
                                       ftp_protect_max = 0.99, ftp_protect_mean = 0.95, ftp_protect_totcount = 100
                                     ),
                                     max_nruns = 7,
                                     max_pareto_k = 10,
                                     ncpu = 1L,
                                     verbose = FALSE,
                                     ...) {
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
  ctables_list <- get_ctable_from_SE(se,
    assayName,
    threshUnmod,
    threshMod,
    max_spacing,
    aggrSamples = F,
    ncpu,
    verbose
  )
  infDFout[["pairStats"]] <- ctables_list

  ## perform parameter inference for each sample
  if (verbose) {
    .warning_timestamp("Performing footprint spectral analysis")
  }

  ## check input parameters for VB

  dots <- list(...)
  if ("iter" %in% names(dots)) {
    iter <- dots$iter
  } else {
    iter <- 5000
  }
  if ("tol_rel_obj" %in% names(dots)) {
    tol_rel_obj <- dots$tol_rel_obj
  } else {
    tol_rel_obj <- 1e-8
  }
  if ("output_samples" %in% names(dots)) {
    output_samples <- dots$output_samples
  } else {
    output_samples <- 4000
  }
  if ("grad_samples" %in% names(dots)) {
    grad_samples <- dots$grad_samples
  } else {
    grad_samples <- 1
  }
  if ("algorithm" %in% names(dots)) {
    algorithm <- dots$algorithm
  } else {
    algorithm <- "meanfield"
  }


  infDF <- do.call(rbind, lapply(
    1:nrow(infDFout),
    function(idx) {
      .message_timestamp(paste0("Inference for ", infDFout$sample[idx]))
      if (verbose) {
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
          ...
        )
      } else {
        invisible(utils::capture.output(
          suppressWarnings(suppressMessages(vb_res <- infer_footprints_vb(ctables_list[[idx]],
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
            ...
          )))
        ))
      }

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

      ## get summary and populate a DataFrame
      if (!is.null(vb_res)) {
        ftpsumm <- get_ftp_inference_summary(vb_res)
        ftp_spec <- ftpsumm$ESTIMATES$ftp_abundance_estimates[, c("ftp_length", "mean", "sd", "2.5%", "50%", "97.5%")]
        DFout$VB_success <- TRUE
        DFout$pareto_k <- vb_res@sim$diagnostics$psis$pareto_k

        ## add inf results for BG emission probs
        if (!is.null(ftpsumm$ESTIMATES$bg_protect_prob_estimate$mean) & !is.na(ftpsumm$ESTIMATES$bg_protect_prob_estimate$mean)) {
          DFout$bg_emis_mean <- ftpsumm$ESTIMATES$bg_protect_prob_estimate[1, "mean"]
          DFout$bg_emis_sd <- ftpsumm$ESTIMATES$bg_protect_prob_estimate[1, "sd"]
          DFout$bg_emis_2.5perc <- ftpsumm$ESTIMATES$bg_protect_prob_estimate[1, "2.5%"]
          DFout$bg_emis_50perc <- ftpsumm$ESTIMATES$bg_protect_prob_estimate[1, "50%"]
          DFout$bg_emis_97.5perc <- ftpsumm$ESTIMATES$bg_protect_prob_estimate[1, "97.5%"]
        }

        ## add inf results for FTP emission probs
        if (!is.null(ftpsumm$ESTIMATES$ftp_protect_prob_estimate$mean) & !is.na(ftpsumm$ESTIMATES$ftp_protect_prob_estimate$mean)) {
          DFout$ftp_emis_mean <- ftpsumm$ESTIMATES$ftp_protect_prob_estimate[1, "mean"]
          DFout$ftp_emis_sd <- ftpsumm$ESTIMATES$ftp_protect_prob_estimate[1, "sd"]
          DFout$ftp_emis_2.5perc <- ftpsumm$ESTIMATES$ftp_protect_prob_estimate[1, "2.5%"]
          DFout$ftp_emis_50perc <- ftpsumm$ESTIMATES$ftp_protect_prob_estimate[1, "50%"]
          DFout$ftp_emis_97.5perc <- ftpsumm$ESTIMATES$ftp_protect_prob_estimate[1, "97.5%"]
        }
        ## add inf results for background
        DFout$bg_coverage_mean <- subset(ftp_spec, ftp_length == 1)[, "mean"]
        DFout$bg_coverage_sd <- subset(ftp_spec, ftp_length == 1)[, "sd"]
        DFout$bg_coverage_2.5perc <- subset(ftp_spec, ftp_length == 1)[, "2.5%"]
        DFout$bg_coverage_50perc <- subset(ftp_spec, ftp_length == 1)[, "50%"]
        DFout$bg_coverage_97.5perc <- subset(ftp_spec, ftp_length == 1)[, "97.5%"]

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
