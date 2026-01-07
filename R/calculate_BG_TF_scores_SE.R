#' Perform log-ratio transformations to predicted posterior coverages
#'
#' @inheritParams predict_footprints_SE
#' @param bgAssayNames vector of assay names that contain predicted coverage
#'   probabilities for background (accessible state)
#' @param tfAssayNames vector of assay names that contain predicted coverage
#'   probabilities for transcription factors
#' @param nuclAssayNames vector of assay names that contain predicted coverage
#'   probabilities for nucleosomes
#' @param psc pseudo-count to avoid log(0)
#'
#' @returns \code{SummarizeExperiment} object with BG and TF scores stored
#' as assays
#'
#' @importFrom SummarizedExperiment SummarizedExperiment assay assayNames
#'     colnames
#' @importFrom SparseArray NaArray
#' @importFrom checkmate assert_vector assert_number
#' @export
#'
calculate_BG_TF_scores_SE <- function(se,
                                      bgAssayNames = c("background_coverProb_nomeR"),
                                      tfAssayNames = c("TF_coverProb_nomeR"),
                                      nuclAssayNames = c("Nucl_coverProb_nomeR"),
                                      mod_probAssayName = "mod_prob",
                                      psc = 0.1) {

    ## check assay names
    assert_vector(x = bgAssayNames, any.missing = FALSE, min.len = 1,
                  unique = TRUE, null.ok = FALSE)
    assert_vector(x = tfAssayNames, any.missing = FALSE, min.len = 1,
                  unique = TRUE, null.ok = FALSE)
    assert_vector(x = nuclAssayNames, any.missing = FALSE, min.len = 1,
                  unique = TRUE, null.ok = FALSE)
    assert_number(x = psc, lower = 0, finite = TRUE)

    stopifnot(all(c(bgAssayNames, tfAssayNames, nuclAssayNames) %in%
                      assayNames(se)))

    ## extract/aggregate bg assays
    sample_names <- colnames(se)
    bg_scores_DF <- make_zero_col_DFrame(nrow = nrow(se))
    tf_scores_DF <- make_zero_col_DFrame(nrow = nrow(se))

    for (snm in sample_names) {
        ## aggregate across assays
        bg_smpl <- Reduce("+", lapply(bgAssayNames, function(as) {
            assay(se, as)[[snm]]
        }))

        tf_smpl <- Reduce("+", lapply(tfAssayNames, function(as) {
            assay(se, as)[[snm]]
        }))

        nucl_smpl <- Reduce("+", lapply(nuclAssayNames, function(as) {
            assay(se, as)[[snm]]
        }))

        ## add pseudo-count
        bg_smpl <- bg_smpl + psc
        tf_smpl <- tf_smpl + psc
        nucl_smpl <- nucl_smpl + psc

        ## calculate log transformations
        bg_score_smpl <- sqrt(2/3) * log(bg_smpl/sqrt(tf_smpl * nucl_smpl))
        tf_score_smpl <- sqrt(1/2) * log(tf_smpl/nucl_smpl)

        bg_scores_DF[[snm]] <- bg_score_smpl
        tf_scores_DF[[snm]] <- tf_score_smpl
    }
    if (mod_probAssayName %in% assayNames(se)) {
        assayList <- list("mod_prob" = assay(se, mod_probAssayName),
                          "BG_score_nomeR" = bg_scores_DF,
                          "TF_score_nomeR" = tf_scores_DF)
        names(assayList) <- mod_probAssayName
    } else{
        assayList <- list("BG_score_nomeR" = bg_scores_DF,
                          "TF_score_nomeR" = tf_scores_DF)
    }

    seOut <- SummarizedExperiment(
        assays = assayList,
        rowRanges = rowRanges(se),
        colData = colData(se),
        metadata = metadata(se)
    )
    return(seOut)
}
