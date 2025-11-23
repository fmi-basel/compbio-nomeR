#' Create count table for 00, 01, 10, 11, etc occurrences SMF data matrix
#'
#' @param data A matrix containing binary SMF data for a ROI. 1s represent 
#'     protected and 0s represent accessible positions.
#' @param max_spacing Maximum spacing between positions.
#' @param ncpu number of cores to use.
#' @param verbose verbose mode for bug fixing.
#'
#' @return data.frame containing frequencies for 0,0; 0,1; etc for each spacing
#' between 1 (total number of 0s and 1 in the data) and \code{max_spacing}.
#'
#' @export
#'
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
