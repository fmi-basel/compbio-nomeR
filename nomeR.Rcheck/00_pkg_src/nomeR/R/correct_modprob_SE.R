
#' Make correction of modification probability
#'
#' @param mod_prob modification probability
#' @param neg_shape1,neg_shape2  shapes of beta distribution from negative control
#' @param pos_shape1,pos_shape2 shapes of beta distribution from positive control
#' @param pos_prior prior probability of positive (modidified) base
#'
#' @returns corrected modification probability
#' @noRd
#' @keywords internal
#' @importFrom stats dbeta plogis
.calc_beta_corrected_mod_prob_dbeta <- function (mod_prob,
                                                 neg_shape1,
                                                 neg_shape2,
                                                 pos_shape1,
                                                 pos_shape2,
                                                 mod_prior = 0.5,
                                                 eps = .Machine$double.eps)
{

    # --- Guard: clamp y away from exact 0/1 ---
    # dbeta(log=TRUE) is finite for y in (0,1) but undefined at boundaries
    mod_prob_safe <- pmax(pmin(mod_prob, 1 - eps), eps)

    boundary_flag <- (mod_prob != mod_prob_safe)      # track which values were clamped

    # --- Log densities ---
    log_f_pos <- dbeta(mod_prob_safe,
                       shape1 = pos_shape1,
                       shape2 = pos_shape2,
                       log    = TRUE)
    log_f_neg <- dbeta(mod_prob_safe,
                       shape1 = neg_shape1,
                       shape2 = neg_shape2,
                       log    = TRUE)

    # --- Prior log-odds ---
    # From per-position theta (clamp theta away from 0/1 for log stability)
    mod_prior_safe    <- pmax(pmin(mod_prior, 1 - eps), eps)
    log_prior_lo  <- log(mod_prior_safe) - log(1 - mod_prior_safe)

    # --- Log posterior odds = log-likelihood ratio + log prior odds ---
    # log_lr = log_f_pos - log_f_neg
    # Cases:
    #   finite  - finite  = finite  → normal
    #   -Inf    - finite  = -Inf    → p(Z=1) = 0
    #   finite  - (-Inf)  = +Inf    → p(Z=1) = 1
    #   -Inf    - (-Inf)  = NaN     → simultaneous underflow → flag
    log_lr            <- log_f_pos - log_f_neg
    simultaneous_flag <- is.nan(log_lr)

    # Resolve NaN: when both densities underflow, fall back to prior
    # (the observation gives no information — trust the prior)?plogis

    log_lr[simultaneous_flag] <- 0

    log_posterior_odds <- log_lr + log_prior_lo

    # --- Sigmoid for final probability (numerically stable) ---
    # sigmoid(x) = 1/(1+exp(-x)), computed carefully for large |x|
    mod_prob_correct <- plogis(log_posterior_odds)   # plogis IS the sigmoid, already stable
    return(mod_prob_correct)
}

#' Standard quantile normalisation
#'
#' @param p_z1   Numeric vector: posterior p(Z=1) (to be normalised)
#' @param y_raw  Numeric vector: reference distribution (raw scores)
#' @return Numeric vector with distribution matched to y_raw
.quantile_normalise <- function(p_z1, y_raw) {
    n        <- length(p_z1)
    # Empirical quantile function of reference
    ref_sorted <- sort(y_raw)
    # Map each posterior to its rank-equivalent raw score
    ranks    <- rank(p_z1, ties.method = "average")
    ref_idx  <- pmax(1L, pmin(n, round(ranks / n * length(ref_sorted))))
    ref_sorted[ref_idx]
}



.correct_sample <- function(dat_na_matrix,
                            neg_shape1_vec,
                            neg_shape2_vec,
                            pos_shape1_vec,
                            pos_shape2_vec,
                            mod_prior = 0.5,
                            eps = .Machine$double.eps,
                            qnorm_to_raw = FALSE){

    nonna_idx <- nnawhich(dat_na_matrix, arr.ind = TRUE)
    ## row; column

    ## create NaArray
    corrected_namat <- NaArray(dim = dim(dat_na_matrix),
                               dimnames = dimnames(dat_na_matrix),
                               type = "double")
    ## add data
    corrected_namat[nonna_idx] <- .calc_beta_corrected_mod_prob_dbeta(mod_prob = dat_na_matrix[nonna_idx],
                                                                      neg_shape1 = neg_shape1_vec[nonna_idx[,1]],
                                                                      neg_shape2 = neg_shape2_vec[nonna_idx[,1]],
                                                                      pos_shape1 = pos_shape1_vec[nonna_idx[,1]],
                                                                      pos_shape2 = pos_shape2_vec[nonna_idx[,1]],
                                                                      mod_prior = mod_prior,
                                                                      eps=eps)
    if(qnorm_to_raw){
        corrected_namat[nonna_idx] <- .quantile_normalise(p_z1 = corrected_namat[nonna_idx],
                                                          y_raw = dat_na_matrix[nonna_idx])
    }

    return(corrected_namat)

}

#' Correct modification probabilities for sequence-context bias
#'
#' @description
#' Applies a per-position Bayesian correction to raw modification probabilities
#' using Beta distributions fitted to negative-control (unmodified) and
#' positive-control (fully modified) samples, stratified by sequence context.
#'
#' @param corrected_assayName Character scalar. Name of the additional assay
#'   added to the returned \code{se} that will contain the corrected
#'   read-level modification probabilities.
#' @param neg_control_shapes A \code{data.table} containing Beta distribution
#'   shape parameters for the negative control experiment (no MTase treatment),
#'   as returned by \code{\link{get_SeqContext_control_beta_shapes_SE}}.
#' @param pos_control_shapes A \code{data.table} containing Beta distribution
#'   shape parameters for the positive control experiment (MTase treatment of
#'   naked DNA), as returned by \code{\link{get_SeqContext_control_beta_shapes_SE}}.
#' @param mod_prior prior probability of a modified base.
#' @param qnorm_to_raw if \code{TRUE} perform quantile normalization of the corrected
#'   probabilities to match the distribution of the raw probabilities.
#'
#' @inheritParams predict_footprints_SE
#'
#' @returns The input \code{se} with an additional assay (named by
#'   \code{corrected_assayName}) containing the corrected modification
#'   probabilities.
#' @export
#' @importFrom SummarizedExperiment SummarizedExperiment rowData colData
#'     assay assayNames assay<-
#' @importFrom SparseArray NaArray nnawhich
#' @importFrom GenomicRanges GPos match seqnames start end
#' @importFrom IRanges subsetByOverlaps
#' @importFrom S4Vectors DataFrame SimpleList metadata
#' @importFrom rlang .data
#' @import data.table
correct_modprob_SE <- function(se,
                               assayName = "mod_prob",
                               corrected_assayName = "mod_prob_corrected",
                               neg_control_shapes,
                               pos_control_shapes,
                               mod_prior = 0.5,
                               eps=.Machine$double.eps,
                               qnorm_to_raw = FALSE){
    #browser()
    stopifnot("sequenceContext" %in% colnames(rowData(se)))
    mod_prob_assays <- assay(se, assayName)


    ## look up sequence contexts for each row in SE in the  neg and pos control tables
    negctrl_other_idx <- which(neg_control_shapes$seqcont == "other")
    posctrl_other_idx <- which(pos_control_shapes$seqcont == "other")
    seqcont_per_row <- rowData(se)[,"sequenceContext"]
    seqcont_per_row <- data.frame(sequenceContext = seqcont_per_row,
                                  rowIndex = 1:length(seqcont_per_row)
    ) %>%
        mutate(neg_control_rowindex = match(.data$sequenceContext,neg_control_shapes$seqcont,
                                            nomatch = negctrl_other_idx
        ),
        pos_control_rowindex = match(.data$sequenceContext,pos_control_shapes$seqcont
                                     ,nomatch = posctrl_other_idx
        ))

    ## extract shapes for pos and neg beta distributions for each row in SE
    neg_shape1_per_row <- neg_control_shapes[["shape1"]][seqcont_per_row$neg_control_rowindex]
    neg_shape2_per_row <- neg_control_shapes[["shape2"]][seqcont_per_row$neg_control_rowindex]
    pos_shape1_per_row <- pos_control_shapes[["shape1"]][seqcont_per_row$pos_control_rowindex]
    pos_shape2_per_row <- pos_control_shapes[["shape2"]][seqcont_per_row$pos_control_rowindex]

    ## correct mod probs for each sample
    assayMat <- make_zero_col_DFrame(nrow = nrow(se))
    for(sI in seq_len(ncol(se))) {
        corrected_data <- .correct_sample(dat_na_matrix = mod_prob_assays[,sI],
                                          neg_shape1_vec = neg_shape1_per_row,
                                          neg_shape2_vec = neg_shape2_per_row,
                                          pos_shape1_vec = pos_shape1_per_row,
                                          pos_shape2_vec = pos_shape2_per_row,
                                          mod_prior = mod_prior,
                                          eps=eps,
                                          qnorm_to_raw = qnorm_to_raw)
        assayMat[[sI]] <- corrected_data
    }
    colnames(assayMat) <- colnames(se)

    ## add mod_prob_corrected assay to se
    assay(se,corrected_assayName,withDimnames=FALSE) <- assayMat
    return(se)
}
