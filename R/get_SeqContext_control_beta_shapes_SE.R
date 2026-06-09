
#' Fit Beta distribution shapes per sequence context using method of moments
#'
#' @param x Numeric vector of modification probabilities.
#' @returns A list with `n_dat`, `alpha`, `beta`.
#' @noRd
#' @keywords internal
#' @importFrom stats var
.fit_beta_mom <- function(x) {
    m  <- mean(x)
    v  <- var(x)
    common <- m * (1 - m) / v - 1
    list(n_dat = length(x), alpha = m * common, beta = (1 - m) * common)
}


#' Assign joint fit-context labels to SE rows
#'
#' Returns the fit-context label for every row in `se`: contexts with at least
#' `min_obs` finite observations in **both** controls keep their own name;
#' all others are mapped to `"OTHER"`.
#'
#' @return A character vector of length `nrow(se)`.
#' @noRd
#' @keywords internal
.joint_fit_ctx <- function(seqcont, neg_ctx, pos_ctx, min_obs) {
    all_ctx <- union(unique(neg_ctx), unique(pos_ctx))
    neg_n   <- tabulate(match(neg_ctx, all_ctx), nbins = length(all_ctx))
    pos_n   <- tabulate(match(pos_ctx, all_ctx), nbins = length(all_ctx))
    names(neg_n) <- names(pos_n) <- all_ctx
    keep_ctx <- all_ctx[neg_n >= min_obs & pos_n >= min_obs]
    ifelse(seqcont %in% keep_ctx, seqcont, "OTHER")
}


#' Fit Beta distribution shapes for positive and negative controls
#'
#' Fits Beta distribution shape parameters (via method of moments) for each
#' sequence context using modification probabilities from a positive and a
#' negative control sample stored in a \code{SummarizedExperiment}.
#'
#' Sequence contexts with fewer than \code{min_obs} finite observations in
#' **either** control are pooled into a single `"OTHER"` category before
#' fitting, so that all downstream correction has a fallback.
#'
#' @param neg_control_sampleName Character scalar. Sample name for the negative
#'   control (SMF experiment without MTase treatment).
#' @param pos_control_sampleName Character scalar. Sample name for the positive
#'   control (SMF experiment with MTase treatment on naked DNA).
#' @param min_obs Integer. Minimum number of finite observations required in
#'   **both** controls for a context to receive its own fitted parameters.
#'   Contexts below this threshold in either control are collapsed into
#'   `"OTHER"`. Default: \code{100}.
#' @inheritParams predict_footprints_SE
#'
#' @returns A \code{data.table} with one row per fitted sequence context and
#'   columns:
#'   \describe{
#'     \item{seqcont}{Context label (`"OTHER"` for pooled sparse contexts).}
#'     \item{n_pos, n_neg}{Number of observations used for fitting.}
#'     \item{alpha_pos, beta_pos}{Beta shape parameters for the positive control.}
#'     \item{alpha_neg, beta_neg}{Beta shape parameters for the negative control.}
#'   }
#' @export
#' @importFrom SummarizedExperiment assay rowData
#' @importFrom SparseArray nnawhich
#' @import data.table
get_SeqContext_control_beta_shapes_SE <- function(se,
                                                   neg_control_sampleName,
                                                   pos_control_sampleName,
                                                   assayName  = "mod_prob",
                                                   min_obs    = 100L) {

    stopifnot("sequenceContext" %in% colnames(rowData(se)))
    stopifnot(all(c(neg_control_sampleName, pos_control_sampleName) %in% colnames(se)))
    min_obs <- as.integer(min_obs)

    seqcont <- as.character(rowData(se)[, "sequenceContext"])

    neg_data <- assay(se, assayName)[, neg_control_sampleName]
    pos_data <- assay(se, assayName)[, pos_control_sampleName]

    neg_nonna_idx <- nnawhich(neg_data, arr.ind = TRUE)
    pos_nonna_idx <- nnawhich(pos_data, arr.ind = TRUE)

    neg_ctx <- seqcont[neg_nonna_idx[, 1]]
    pos_ctx <- seqcont[pos_nonna_idx[, 1]]

    row_fit_ctx  <- .joint_fit_ctx(seqcont, neg_ctx, pos_ctx, min_obs)
    neg_fit_ctx  <- row_fit_ctx[neg_nonna_idx[, 1]]
    pos_fit_ctx  <- row_fit_ctx[pos_nonna_idx[, 1]]

    neg_vals <- neg_data[neg_nonna_idx]
    pos_vals <- pos_data[pos_nonna_idx]

    neg_split <- split(neg_vals, neg_fit_ctx)
    pos_split <- split(pos_vals, pos_fit_ctx)

    fit_ctxs <- union(names(neg_split), names(pos_split))

    rbindlist(lapply(fit_ctxs, function(ctx) {
        nv <- neg_split[[ctx]]
        pv <- pos_split[[ctx]]
        neg_sh <- if (length(nv) >= 2L) .fit_beta_mom(nv) else list(n_dat=0L, alpha=1, beta=1)
        pos_sh <- if (length(pv) >= 2L) .fit_beta_mom(pv) else list(n_dat=0L, alpha=1, beta=1)
        data.table(
            seqcont   = ctx,
            n_pos     = as.integer(pos_sh$n_dat),
            alpha_pos = pos_sh$alpha,
            beta_pos  = pos_sh$beta,
            n_neg     = as.integer(neg_sh$n_dat),
            alpha_neg = neg_sh$alpha,
            beta_neg  = neg_sh$beta
        )
    }))
}
