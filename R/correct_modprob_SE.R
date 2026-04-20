
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
#' @importFrom stats pbeta
.calc_beta_corrected_mod_prob <- function(mod_prob,
                                         neg_shape1,
                                         neg_shape2,
                                         pos_shape1,
                                         pos_shape2,
                                         pos_prior = 0.5){

    posbeta <- pos_prior * pbeta(q = mod_prob,
                                 shape1 = pos_shape1,
                                 shape2 = pos_shape2)
    negbeta <- (1 - pos_prior) * (1 - pbeta(q = mod_prob,
                                            shape1 = neg_shape1,
                                            shape2 = neg_shape2))
    mod_prob_correct <- posbeta/(posbeta + negbeta)
    return(mod_prob_correct)
}

.correct_sample <- function(dat_na_matrix,
                            seqcontext,
                            neg_control_shapes,
                            pos_control_shapes,
                            pos_prior = 0.5){

    nonna_idx <- nnawhich(dat_na_matrix, arr.ind = TRUE)


    ## create NAArray
    corrected_namat <- NaArray(dim = dim(dat_na_matrix),
                     dimnames = dimnames(dat_na_matrix),
                     type = "double")
    ## add data
    corrected_namat[nonna_idx] <- .calc_beta_corrected_mod_prob(mod_prob = dat_na_matrix[as.matrix(nonna_idx)],
                                                                neg_shape1 = neg_control_shapes$shape1[seqcontext[,"neg_control_rowindex"]],
                                                                neg_shape2 = neg_control_shapes$shape2[seqcontext[,"neg_control_rowindex"]],
                                                                pos_shape1 = pos_control_shapes$shape1[seqcontext[,"pos_control_rowindex"]],
                                                                pos_shape2 = pos_control_shapes$shape2[seqcontext[,"pos_control_rowindex"]],
                                                                pos_prior = pos_prior)
    return(corrected_namat)

}

#' Correct modification probabilities
#'

#' @param corrected_assayName Character scalar specifying the name of the additional
#'   assay in returned \code{se} that will contain read-level corrected modification probabilities.
#' @param neg_control_shapes \code{data.table}'s containing shapes of Beta distributions
#'   for negative control experiments. Can be obtained using \code{\link{get_SeqContext_control_beta_shapes_SE}}.
#' @param pos_control_shapes \code{data.table}'s containing shapes of Beta distributions
#'   for positive control experiments. Can be obtained using \code{\link{get_SeqContext_control_beta_shapes_SE}}.
#'
#' @inheritParams predict_footprints_SE
#'
#' @returns \code{se} with additional assay (with name defined by \code{corrected_assayName}) containing
#'   corrected modification probabilities.
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
                                pos_control_shapes){
    stopifnot("sequenceContext" %in% colnames(rowData(se)))
    mod_prob_assays <- assay(se, assayName)


    ## make correpsondence between sequence context and
    ## vocabulary in neg and pos shapes


    negctrl_other_idx <- which(neg_control_shapes$seqcont == "other")
    posctrl_other_idx <- which(pos_control_shapes$seqcont == "other")
    seqcont <- rowData(se)[,"sequenceContext"]
    seqcont <- data.frame(sequenceContext = seqcont,
                          rowIndex = 1:length(seqcont)
                          ) %>%
        mutate(neg_control_rowindex = match(.data$sequenceContext,neg_control_shapes$seqcont,nomatch = negctrl_other_idx),
               pos_control_rowindex = match(.data$sequenceContext,neg_control_shapes$seqcont,nomatch = posctrl_other_idx))


    ## correct mod probs for each sample
    assayMat <- make_zero_col_DFrame(nrow = nrow(se))
    for(sI in seq_len(ncol(se))) {
        corrected_data <- .correct_sample(dat_na_matrix = mod_prob_assays[,sI],
                                          seqcontext = seqcont,
                                          neg_control_shapes = neg_control_shapes,
                                          pos_control_shapes = pos_control_shapes)
        assayMat[[sI]] <- corrected_data
    }
    colnames(assayMat) <- colnames(se)

    ## add mod_prob_corrected assay to se
    assay(se,corrected_assayName) <- assayMat
    return(se)
}
