
#' Bayesian posterior correction using Beta distributions (BetaCorrect)
#'
#' @noRd
#' @keywords internal
#' @importFrom stats dbeta plogis
.calc_beta_corrected_mod_prob_dbeta <- function(mod_prob,
                                                 alpha_neg,
                                                 beta_neg,
                                                 alpha_pos,
                                                 beta_pos,
                                                 mod_prior = 0.5,
                                                 eps = .Machine$double.eps) {
    mod_prob_safe <- pmax(pmin(mod_prob, 1 - eps), eps)
    log_f_pos <- dbeta(mod_prob_safe, shape1 = alpha_pos, shape2 = beta_pos, log = TRUE)
    log_f_neg <- dbeta(mod_prob_safe, shape1 = alpha_neg, shape2 = beta_neg, log = TRUE)
    mod_prior_safe    <- pmax(pmin(mod_prior, 1 - eps), eps)
    log_prior_lo      <- log(mod_prior_safe) - log(1 - mod_prior_safe)
    log_lr            <- log_f_pos - log_f_neg
    log_lr[is.nan(log_lr)] <- 0
    plogis(log_lr + log_prior_lo)
}


#' Standard quantile normalisation
#' @noRd
#' @keywords internal
.quantile_normalise <- function(p_z1, y_raw) {
    n          <- length(p_z1)
    ref_sorted <- sort(y_raw)
    ranks      <- rank(p_z1, ties.method = "average")
    ref_idx    <- pmax(1L, pmin(n, round(ranks / n * length(ref_sorted))))
    ref_sorted[ref_idx]
}


#' Resolve sequence context to row index in control_params
#'
#' Returns a 1-based integer index into `control_params` for every row in `se`.
#' Falls back to the "OTHER" row when an exact match is not found.
#' Rows with no match at all get `NA_integer_`.
#' @noRd
#' @keywords internal
.resolve_context_idx <- function(seqcont_per_row, control_params) {
    idx       <- match(seqcont_per_row, control_params$seqcont)
    other_idx <- match("OTHER", control_params$seqcont)
    if (!is.na(other_idx)) idx[is.na(idx)] <- other_idx
    idx
}


#' Warn about rows with no matching context model
#' @noRd
#' @keywords internal
.warn_no_model <- function(seqcont_per_row, param_idx_per_row) {
    no_model <- is.na(param_idx_per_row)
    if (any(no_model)) {
        missing_ctx <- unique(seqcont_per_row[no_model])
        cli::cli_warn(
            "{sum(no_model)} row(s) across {length(missing_ctx)} context(s) \\
             have no model and will retain raw values: \\
             {paste(missing_ctx, collapse = ', ')}")
    }
    no_model
}


# ----------------------------------------------------------------
# Per-sample correction helpers
# ----------------------------------------------------------------

#' BetaCorrect correction for one sample column
#' @noRd
#' @keywords internal
.correct_sample_bc <- function(dat_na_matrix,
                                alpha_neg_vec,
                                beta_neg_vec,
                                alpha_pos_vec,
                                beta_pos_vec,
                                mod_prior    = 0.5,
                                eps          = .Machine$double.eps,
                                qnorm_to_raw = FALSE) {
    nonna_idx <- nnawhich(dat_na_matrix, arr.ind = TRUE)
    corrected_namat <- NaArray(dim = dim(dat_na_matrix),
                               dimnames = dimnames(dat_na_matrix),
                               type = "double")
    if (nrow(nonna_idx) == 0L) return(corrected_namat)
    corrected_namat[nonna_idx] <- .calc_beta_corrected_mod_prob_dbeta(
        mod_prob  = dat_na_matrix[nonna_idx],
        alpha_neg = alpha_neg_vec[nonna_idx[, 1]],
        beta_neg  = beta_neg_vec[nonna_idx[, 1]],
        alpha_pos = alpha_pos_vec[nonna_idx[, 1]],
        beta_pos  = beta_pos_vec[nonna_idx[, 1]],
        mod_prior = mod_prior,
        eps       = eps)
    if (qnorm_to_raw)
        corrected_namat[nonna_idx] <- .quantile_normalise(
            p_z1  = corrected_namat[nonna_idx],
            y_raw = dat_na_matrix[nonna_idx])
    corrected_namat
}


#' BetaUniform correction for one sample column
#' @noRd
#' @keywords internal
.correct_sample_bu <- function(dat_na_matrix,
                                param_idx_per_row,
                                control_params,
                                pi_pos       = 0.5,
                                isotonic     = TRUE,
                                qnorm_to_raw = FALSE) {
    nonna_idx <- nnawhich(dat_na_matrix, arr.ind = TRUE)
    corrected_namat <- NaArray(dim = dim(dat_na_matrix),
                               dimnames = dimnames(dat_na_matrix),
                               type = "double")
    if (nrow(nonna_idx) == 0L) return(corrected_namat)

    raw_p         <- dat_na_matrix[nonna_idx]
    obs_param_idx <- param_idx_per_row[nonna_idx[, 1]]

    corrected_vals <- .correct_posterior_cpp(
        raw_p     = raw_p,
        param_idx = obs_param_idx,
        alpha_pos = control_params$alpha_pos,
        beta_pos  = control_params$beta_pos,
        eps_pos   = control_params$eps_pos,
        alpha_neg = control_params$alpha_neg,
        beta_neg  = control_params$beta_neg,
        eps_neg   = control_params$eps_neg,
        pi_pos    = pi_pos,
        isotonic  = isotonic)

    no_model_obs <- is.na(obs_param_idx)
    if (any(no_model_obs)) corrected_vals[no_model_obs] <- raw_p[no_model_obs]

    corrected_namat[nonna_idx] <- corrected_vals
    if (qnorm_to_raw)
        corrected_namat[nonna_idx] <- .quantile_normalise(
            p_z1  = corrected_namat[nonna_idx],
            y_raw = raw_p)
    corrected_namat
}


# ----------------------------------------------------------------
# Public function
# ----------------------------------------------------------------

#' Correct modification probabilities for sequence-context bias
#'
#' @description
#' Applies per-position Bayesian correction to raw modification probabilities
#' using fitted control models stratified by sequence context.
#'
#' Two methods are supported:
#' \describe{
#'   \item{\code{"BetaCorrect"}}{Bayesian posterior using Beta distributions fitted
#'     by method of moments to positive and negative controls. Parameters are
#'     obtained from \code{\link{get_SeqContext_control_beta_shapes_SE}}.}
#'   \item{\code{"BetaUniform"}}{Bayesian posterior using Beta-Uniform mixture
#'     models fitted by constrained EM+Newton-Raphson. Includes built-in
#'     pool-adjacent-violators isotonic regression to enforce monotonicity
#'     between raw and corrected probabilities. Parameters are obtained from
#'     \code{\link{fit_SeqContext_BetaUnif_params_SE}}.}
#' }
#'
#' @param control_params A \code{data.table} returned by either
#'   \code{\link{get_SeqContext_control_beta_shapes_SE}} (for
#'   \code{method = "BetaCorrect"}) or
#'   \code{\link{fit_SeqContext_BetaUnif_params_SE}} (for
#'   \code{method = "BetaUniform"}).
#' @param method Character scalar. Correction method: \code{"BetaCorrect"} or
#'   \code{"BetaUniform"}.
#' @param corrected_assayName Character scalar. Name of the additional assay
#'   added to the returned \code{se}.
#' @param mod_prior Numeric in (0, 1). Prior probability of a modified
#'   (accessible) base. Used as \code{pi_pos} in the Bayesian posterior.
#'   Default: \code{0.5}.
#' @param eps Numeric. Clamping tolerance for boundary probabilities (only used
#'   by \code{method = "BetaCorrect"}). Default: \code{.Machine$double.eps}.
#' @param isotonic Logical. If \code{TRUE} (default), apply pool-adjacent-
#'   violators (PAV) isotonic regression per sequence context after computing
#'   the Bayesian posterior, enforcing monotonicity between raw and corrected
#'   probabilities. Set to \code{FALSE} to return the raw posterior without
#'   isotonic smoothing. Only applicable when \code{method = "BetaUniform"}.
#' @param qnorm_to_raw Logical. If \code{TRUE}, apply quantile normalisation so
#'   the corrected probabilities match the marginal distribution of the raw
#'   probabilities. Default: \code{TRUE}.
#'
#' @inheritParams predict_footprints_SE
#'
#' @returns The input \code{se} with an additional assay (named by
#'   \code{corrected_assayName}) containing the corrected modification
#'   probabilities.
#' @export
#' @importFrom SummarizedExperiment rowData assay assay<- assayNames
#' @importFrom SparseArray NaArray nnawhich
#' @importFrom S4Vectors make_zero_col_DFrame
#' @import data.table
correct_modprob_SE <- function(se,
                               control_params,
                               method              = c("BetaCorrect", "BetaUniform"),
                               assayName           = "mod_prob",
                               corrected_assayName = "mod_prob_corrected",
                               mod_prior           = 0.5,
                               eps                 = .Machine$double.eps,
                               isotonic            = TRUE,
                               qnorm_to_raw        = TRUE) {

    method <- match.arg(method)

    stopifnot("sequenceContext" %in% colnames(rowData(se)))
    stopifnot(is.data.frame(control_params) || is.data.table(control_params))
    stopifnot("seqcont" %in% names(control_params))

    if (method == "BetaCorrect") {
        req_cols <- c("alpha_pos", "beta_pos", "alpha_neg", "beta_neg")
        missing  <- setdiff(req_cols, names(control_params))
        if (length(missing))
            cli::cli_abort("control_params is missing columns for method 'BetaCorrect': {paste(missing, collapse=', ')}")
    } else {
        req_cols <- c("alpha_pos", "beta_pos", "eps_pos", "alpha_neg", "beta_neg", "eps_neg")
        missing  <- setdiff(req_cols, names(control_params))
        if (length(missing))
            cli::cli_abort("control_params is missing columns for method 'BetaUniform': {paste(missing, collapse=', ')}")
    }

    seqcont_per_row  <- as.character(rowData(se)[, "sequenceContext"])
    param_idx_per_row <- .resolve_context_idx(seqcont_per_row, control_params)
    no_model          <- .warn_no_model(seqcont_per_row, param_idx_per_row)

    mod_prob_assays <- assay(se, assayName)

    assayMat <- make_zero_col_DFrame(nrow = nrow(se))

    if (method == "BetaCorrect") {
        alpha_neg_per_row <- control_params$alpha_neg[param_idx_per_row]
        beta_neg_per_row  <- control_params$beta_neg[param_idx_per_row]
        alpha_pos_per_row <- control_params$alpha_pos[param_idx_per_row]
        beta_pos_per_row  <- control_params$beta_pos[param_idx_per_row]

        for (sI in seq_len(ncol(se))) {
            cd <- .correct_sample_bc(
                dat_na_matrix = mod_prob_assays[, sI],
                alpha_neg_vec = alpha_neg_per_row,
                beta_neg_vec  = beta_neg_per_row,
                alpha_pos_vec = alpha_pos_per_row,
                beta_pos_vec  = beta_pos_per_row,
                mod_prior     = mod_prior,
                eps           = eps,
                qnorm_to_raw  = qnorm_to_raw)
            assayMat[[sI]] <- cd
        }

    } else {
        for (sI in seq_len(ncol(se))) {
            cd <- .correct_sample_bu(
                dat_na_matrix     = mod_prob_assays[, sI],
                param_idx_per_row = param_idx_per_row,
                control_params    = control_params,
                pi_pos            = mod_prior,
                isotonic          = isotonic,
                qnorm_to_raw      = qnorm_to_raw)
            assayMat[[sI]] <- cd
        }
    }

    colnames(assayMat) <- colnames(se)
    assay(se, corrected_assayName, withDimnames = FALSE) <- assayMat
    se
}
