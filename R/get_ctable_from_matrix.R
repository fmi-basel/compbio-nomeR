#' Create co-occurrence count tables (00, 01, 10, 11, etc.) from a binary SMF data matrix
#'
#' @description
#' Computes co-occurrence frequencies of binary accessibility states
#' (`00`, `01`, `10`, `11`) at varying spacings within a matrix representing
#' a single-molecule footprinting (SMF) region of interest (ROI).
#' Rows correspond to individual molecules (or fragments) and columns
#' correspond to genomic positions. In the matrix, `1` represents a
#' protected position and `0` represents an accessible position.
#'
#' @param data A binary matrix of SMF data for a single ROI.
#' @param max_spacing Integer specifying the maximum spacing (in bases)
#'     between positions for which co-occurrence frequencies are computed.
#' @param ncpu Number of CPU cores to use for parallel computation.
#' @param verbose Logical. If \code{TRUE}, prints progress and diagnostic
#'     messages for debugging.
#'
#' @return
#' A \code{data.frame} containing co-occurrence frequencies for `00`, `01`,
#' `10`, and `11` at each spacing from 1 up to \code{max_spacing}.
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
        data[["nonNA_data"]][, "protect"],   ## binary protection data, 0 - accessible, 1 - protected
        max_spacing,
        ncpu,
        verbose)
    return(out_data)
}
