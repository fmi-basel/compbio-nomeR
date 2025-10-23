#' Create count table for 00, 01, 10, 11, etc occurences in NOMe-Seq data
#'
#' @param data A matrix containing NOMe-Seq data for an amplicon.
#' @param max_spacing Maximum spacing between positions.
#' @param ncpu number of cores to use.
#'
#' @return data.frame containing frequencies for 0,0; 0,1; etc for each spacing
#' between 1 (total number of 0s and 1 in the data) and \code{max_spacing}.
#'
#' @export
#'
get_ctable_from_matrix <- function(data, max_spacing,
                                   ncpu = 1L) {

    ### validate data
    data <- validate_data_for_predict2(data)
    out_data <- count_spacing_freq_cpp(data[["nonNA_data"]][,"fidx_glob"], ## unique fragment ID or index
    																	 data[["nonNA_data"]][,"fragpos"],      ## position within fragment, 0 - based
    																	 data[["nonNA_data"]][,"protect"],   ## binary protection data, 0 - accessible, 1 - protected
                                       max_spacing)
    return(as.data.frame(out_data))
}
