#' Create count table for 00, 01, 10, 11, etc occurences in NOMe-Seq data
#'
#' @param data A matrix containing NOMe-Seq data for an amplicon.
#' @param max_spacing Maximum spacing between positions.
#' @param max_ftplen A priori maximum length of possible footprint in the data.
#'     If \code{max_ftplen} is not 0, then the sequences are extended by 
#'     \code{max_ftplen} on the left and on the right and filled with NA.
#' @param ncpu number of cores to use.
#'
#' @return data.frame containing frequencies for 0,0; 0,1; etc for each spacing 
#' between 1 (total number of 0s and 1 in the data) and \code{max_spacing}.
#'
#' @export
#'
get_ctable_from_matrix <- function(data, max_spacing, max_ftplen = 0L,
                                   ncpu = 1L) {
    
    ### validate data
    data <- validate_data_for_predict(data)
    out_data <- count_spacing_freq_cpp(data[["data_list"]],
                                       data[["fragnames"]],
                                       max_spacing,
                                       max_ftplen)
    return(as.data.frame(out_data))
}
