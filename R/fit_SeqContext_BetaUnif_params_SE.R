
#' Fit Beta-Uniform mixture models for positive and negative controls
#'
#' Fits per-context Beta-Uniform mixture parameters using joint constrained EM
#' (Expectation-Maximisation + Newton-Raphson), with a monotone-likelihood-
#' ratio-preserving (MLRP) projection at each step. The constraint
#' \code{alpha_pos >= alpha_neg} and \code{beta_pos <= beta_neg} guarantees
#' that the resulting posterior correction is monotone in the raw probability.
#'
#' Sequence contexts with fewer than \code{min_obs} finite observations in
#' **either** control are pooled into a single `"OTHER"` category before
#' fitting.
#'
#' @param neg_control_sampleName Character scalar. Sample name for the negative
#'   control (SMF experiment without MTase treatment).
#' @param pos_control_sampleName Character scalar. Sample name for the positive
#'   control (SMF experiment with MTase treatment on naked DNA).
#' @param min_obs Integer. Minimum number of finite observations required in
#'   **both** controls for a context to receive its own model. Default: \code{100}.
#' @param eps_max_pos Numeric. Maximum uniform weight for the positive control
#'   mixture. Default: \code{0.15}.
#' @param eps_max_neg Numeric. Maximum uniform weight for the negative control
#'   mixture. Default: \code{0.25}. Higher values prevent the negative-control
#'   density from collapsing to zero near the accessible end.
#' @param ncpu Integer. Number of parallel workers for context fitting via
#'   \code{parallel::mclapply} (Unix/macOS only). Default: \code{1L}.
#' @inheritParams predict_footprints_SE
#'
#' @returns A \code{data.table} with one row per fitted sequence context and
#'   columns:
#'   \describe{
#'     \item{seqcont}{Context label (`"OTHER"` for pooled sparse contexts).}
#'     \item{alpha_pos, beta_pos, eps_pos}{Beta-Uniform parameters for the positive control.}
#'     \item{alpha_neg, beta_neg, eps_neg}{Beta-Uniform parameters for the negative control.}
#'     \item{n_pos, n_neg}{Number of observations used for fitting.}
#'     \item{converged}{Logical: did the EM converge?}
#'     \item{mlrp_ok}{Logical: is the MLRP constraint satisfied after fitting?}
#'   }
#' @export
#' @importFrom SummarizedExperiment assay rowData
#' @importFrom SparseArray nnawhich
#' @import data.table
#' @import parallel
fit_SeqContext_BetaUnif_params_SE <- function(se,
                                               neg_control_sampleName,
                                               pos_control_sampleName,
                                               assayName   = "mod_prob",
                                               min_obs     = 100L,
                                               eps_max_pos = 0.15,
                                               eps_max_neg = 0.25,
                                               ncpu        = 1L) {

    stopifnot("sequenceContext" %in% colnames(rowData(se)))
    stopifnot(all(c(neg_control_sampleName, pos_control_sampleName) %in% colnames(se)))
    min_obs <- as.integer(min_obs)
    ncpu    <- as.integer(ncpu)

    seqcont <- as.character(rowData(se)[, "sequenceContext"])

    neg_data <- assay(se, assayName)[, neg_control_sampleName]
    pos_data <- assay(se, assayName)[, pos_control_sampleName]

    neg_nonna_idx <- nnawhich(neg_data, arr.ind = TRUE)
    pos_nonna_idx <- nnawhich(pos_data, arr.ind = TRUE)

    neg_ctx <- seqcont[neg_nonna_idx[, 1]]
    pos_ctx <- seqcont[pos_nonna_idx[, 1]]

    row_fit_ctx <- .joint_fit_ctx(seqcont, neg_ctx, pos_ctx, min_obs)
    neg_fit_ctx <- row_fit_ctx[neg_nonna_idx[, 1]]
    pos_fit_ctx <- row_fit_ctx[pos_nonna_idx[, 1]]

    neg_vals  <- neg_data[neg_nonna_idx]
    pos_vals  <- pos_data[pos_nonna_idx]

    neg_split <- split(neg_vals, neg_fit_ctx)
    pos_split <- split(pos_vals, pos_fit_ctx)

    fit_ctxs <- union(names(neg_split), names(pos_split))

    .fit_one <- function(ctx) {
        nv <- neg_split[[ctx]]; if (is.null(nv)) nv <- numeric(0)
        pv <- pos_split[[ctx]]; if (is.null(pv)) pv <- numeric(0)
        if (length(nv) < 2L || length(pv) < 2L) {
            return(data.table(
                seqcont   = ctx,
                alpha_pos = NA_real_, beta_pos = NA_real_, eps_pos = NA_real_,
                alpha_neg = NA_real_, beta_neg = NA_real_, eps_neg = NA_real_,
                n_pos     = length(pv), n_neg = length(nv),
                converged = FALSE, mlrp_ok = FALSE))
        }
        fit <- .fit_one_context_cpp(pv, nv, eps_max_pos, eps_max_neg)
        data.table(
            seqcont   = ctx,
            alpha_pos = fit$alpha_pos, beta_pos = fit$beta_pos, eps_pos = fit$eps_pos,
            alpha_neg = fit$alpha_neg, beta_neg = fit$beta_neg, eps_neg = fit$eps_neg,
            n_pos     = length(pv),    n_neg    = length(nv),
            converged = fit$converged, mlrp_ok  = fit$mlrp_ok
        )
    }

    use_parallel <- ncpu > 1L && .Platform$OS.type != "windows"
    rows <- if (use_parallel)
        parallel::mclapply(fit_ctxs, .fit_one, mc.cores = ncpu)
    else
        lapply(fit_ctxs, .fit_one)

    rbindlist(rows)
}
