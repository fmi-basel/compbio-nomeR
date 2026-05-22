#' Create co-occurrence count tables (N00, N01, N10, N11) from an SMF data matrix
#'
#' @description
#' Computes expected co-occurrence counts for accessibility state pairs
#' (\code{N00}, \code{N01}, \code{N10}, \code{N11}) at varying spacings within
#' a matrix of single-molecule footprinting (SMF) modification probabilities.
#' Rows correspond to individual molecules (or fragments) and columns correspond
#' to genomic positions. Each cell holds a modification probability in \code{[0, 1]},
#' where high values indicate accessible (methylated) positions and low values
#' indicate protected (unmethylated) positions. \code{NA} encodes missing data.
#' Expected counts are computed as sums of products of probabilities:
#' e.g. \code{N00 += p_i * p_j} for each valid pair at spacing \code{s}.
#'
#' @param data A numeric matrix of SMF modification probabilities for a single ROI.
#'   Values must be in \code{[0, 1]} or \code{NA}.
#' @param max_spacing Integer specifying the maximum spacing (in bases)
#'     between positions for which co-occurrence counts are computed.
#' @param ncpu Number of CPU cores to use for parallel computation.
#' @param verbose Logical. If \code{TRUE}, prints progress and diagnostic
#'     messages for debugging.
#'
#' @return
#' A \code{data.frame} containing expected co-occurrence counts for \code{N00},
#' \code{N01}, \code{N10}, and \code{N11} at each spacing from 1 up to
#' \code{max_spacing}.
#'
#' @export
get_ctable_from_matrix <- function(data,
                                   max_spacing = 200L,
                                   ncpu = 1L,
                                   verbose = FALSE) {

    ### validate data
    data <- validate_prepare_listOrMat(data)
    out_data <- count_spacing_freq_cpp(
        data[["nonNA_data"]][, "fidx_glob"], ## unique fragment ID or index
        data[["nonNA_data"]][, "fragpos"],   ## position within fragment, 1 - based
        data[["nonNA_data"]][, "mod_prob"],  ## modification probability in [0,1]
        max_spacing,
        ncpu,
        verbose)
    return(out_data)
}
