




#' Truncate protected positions from left and right sides
#' 
#' This function for each row in a matrix containing NOMe-Seq data fills with NA all protected positions from left and right sides until it meets \code{trunc_until_nbg} unprotected positions.
#' The main aim of this function is to get read of incomplete fottprints which are located at the edges of sequences fragment.
#'
#' @param data_matr matrix containing NOMe-Seq data for an amplicon.
#' @param trunc_until_nbg number of unprotected positions to stop truncation.
#'
#' @return same matrix as `data_matr` but flanking protected positions are replaced by `NA`.
#' @export
#'
truncate_footprint_overhangs <- function(data_matr,
                                         trunc_until_nbg = 1L){
  invert_data <- 1 - data_matr
  invert_data[is.na(invert_data)] <- 0
  
  trunc_data <- do.call(rbind,lapply(1:nrow(data_matr),
                                     function(irow){
                                       row_dat <- data_matr[irow,]
                                       cols_to_truncate <- cumsum(invert_data[irow,]) < trunc_until_nbg | rev(cumsum(rev(invert_data[irow,]))) < trunc_until_nbg
                                       if(any(cols_to_truncate)){
                                         row_dat[cols_to_truncate] <- NA
                                       }
                                       return(row_dat)
                                     }))
  return(trunc_data)
}


