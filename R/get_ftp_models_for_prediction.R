#' Create parameters for \code{\link{predict_footprints}} and \code{\link{predict_footprints_SE}}
#' from a footprint inference summary
#'
#' @description
#' Utility function to generate a list of parameters required by
#' \code{\link{predict_footprints}} or \code{\link{predict_footprints_SE}}
#' using a footprint inference summary obtained from
#' \code{\link{get_ftp_inference_summary}}. These parameters include
#' footprint models, background protection probabilities, and coverage priors.
#'
#' @param infer_summary A \code{list} returned by
#'   \code{\link{get_ftp_inference_summary}}, which must contain
#'   \code{"ESTIMATES"} and \code{"FTP_SUGGEST"}.
#'   If you wish to override the suggested footprints, you can replace
#'   the \code{"FTP_SUGGEST"} element with your own \code{n x 2} matrix.
#'
#' @return A \code{list} containing the following elements:
#' \describe{
#'   \item{\code{FTP_MODELS}}{Footprint models suitable for the
#'     \code{footprint_models} parameter in \code{\link{predict_footprints}}.}
#'   \item{\code{bgprotectprob}}{Numeric value of background protection
#'     probability, corresponding to the \code{bgprotectprob} parameter in
#'     \code{\link{predict_footprints}}.}
#'   \item{\code{bgcoverprior}}{Numeric value of background coverage prior,
#'     corresponding to the \code{bgcoverprior} parameter in
#'     \code{\link{predict_footprints}}.}
#' }
#'
#' @export

# @examples
get_ftp_models_for_prediction <- function(infer_summary) {

    if (is.null(infer_summary[["ESTIMATES"]]) ||
        is.null(infer_summary[["FTP_SUGGEST"]])) {
        stop("Incorrect input. infer_summary must be a list containing ",
             "ESTIMATES and FTP_SUGGEST. See ?nomeR::get_inference_summary")
    }

    if (is.null(row.names(infer_summary[["FTP_SUGGEST"]]))) {
        stop("Can't find 'row.names' of matrix 'FTP_SUGGEST'. Please use ",
             "footprint names as 'row.names'.")
    }
    infer_ftp_abund_probs <-
        infer_summary[["ESTIMATES"]][["ftp_abundance_estimates"]]
    if (inherits(infer_summary[["ESTIMATES"]][["ftp_protect_prob_estimate"]],
                 "data.frame")) {
        ftp_protect_prob <-
            infer_summary[["ESTIMATES"]][["ftp_protect_prob_estimate"]][1,
                                                                        "mean"]
    } else {
        ftp_protect_prob <-
            infer_summary[["ESTIMATES"]][["ftp_protect_prob_estimate"]]
    }

    if (inherits(infer_summary[["ESTIMATES"]][["bg_protect_prob_estimate"]],
                 "data.frame")) {
        bg_protect_prob <-
            infer_summary[["ESTIMATES"]][["bg_protect_prob_estimate"]][1,
                                                                       "mean"]
    } else {
        bg_protect_prob <-
            infer_summary[["ESTIMATES"]][["bg_protect_prob_estimate"]]
    }

    if (is.na(bg_protect_prob)) {
        warning("bg_protect_prob in infer_summary is NA. Likely inference ",
                "was done with background_model = \"fixed\". Use parameter ",
                "\"bg_protect_prob_fixed\" which was used for inference as ",
                "parameter bgprotectprob in predict_footprints.")
    }
    ftp_suggest <- infer_summary[["FTP_SUGGEST"]]

    ftp_lengths <- sapply(row.names(ftp_suggest),
                          function(ridx) {
                              seq(from = ftp_suggest[ridx, 1],
                                  to = ftp_suggest[ridx, 2],
                                  by = 1)
                          }, simplify = FALSE, USE.NAMES = TRUE)
    ## select required footprints
    req_infer_ftp_lengths <-
        infer_ftp_abund_probs[infer_ftp_abund_probs[, "ftp_length"] %in%
                                  c(1, unlist(ftp_lengths)), ]
    req_infer_ftp_lengths[,"mean"] <- req_infer_ftp_lengths[, "mean"] /
        sum(req_infer_ftp_lengths[, "mean"])

    ## create models
    ftp_models <- do.call(
        c,
        lapply(
            names(ftp_lengths),
            function(nm) {
                lapply(
                    ftp_lengths[[nm]],
                    function(flen) {
                        list("PROTECT_PROB" = rep(ftp_protect_prob,flen),
                             "COVER_PRIOR" =
                                 req_infer_ftp_lengths[
                                     req_infer_ftp_lengths[, "ftp_length"] ==
                                         flen,
                                     "mean"],
                             "NAME" = paste0(nm, "--", flen))
                    })
            }))

    return(list(
        "FTP_MODELS" = ftp_models,
        "bgprotectprob" = bg_protect_prob,
        "bgcoverprior" =
            req_infer_ftp_lengths[req_infer_ftp_lengths[, "ftp_length"] == 1,
                                  "mean"]))
}
