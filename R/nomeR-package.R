## R functions for running nomeR

#' nomeR: A package that uses a probabilistic model for interpreting NOME-Seq data
#'
#' The package uses a probabilistic model for calculating posterior probabilities that a certain position 
#' in a NOME-Seq read is covered by a certain footprint.
#' 
#' @section nomeR functions:
# \code{create_data_list} - function for preprocessing a NOME-Seq matrix containing information about
# "protected" and "unprotected" positions in each read. Rows must represent reads (fragments) and columns must represent positions in the amplicon.
#' 
#'   
#' \code{run_nomeR} - function which takes as input 1) a matrix with NOME-Seq data where rows represent reads (or fragments)
#' and columns represent positions in the amplicon, 2) models for a set of user-defined footprints,
#' and using a probabilistic model it returns a data.frame with posterior probabilities for each position and fragment to be covered by each footprint.
#' @docType package
#' @name nomeR
#' 
#' 
#' 
#' 
NULL


## usethis namespace: start
#' @useDynLib nomeR, .registration = TRUE
## usethis namespace: end
NULL


## usethis namespace: start
#' @importFrom Rcpp sourceCpp
## usethis namespace: end
NULL













