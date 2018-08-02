## R functions for running rnome



#' Function for running probabilistic model of interpretation of NOMe-Seq data
#'
#' @param data A list containing preprocessed NOMe-Seq data returned from \code{create_data_list}.
#' @param binding_models A list containing binding models for proteins. 
#' Each element must have 3 slops:
#' PROTECT_PROB - a numeric vector with probability to find protected 'C' within the footprint
#' PRIOR - prior (numeric) probability reflecting how often we expect to find a footprint
#' NAME - name (character) of a model, e.g. "Nucleosome"
#' @param bgprotectprob probability to find a protected 'C' in open (in reality uprotected) regions
#' @param bgprior prior probability for a free (unprotected, or background) position
#' @param bound_fit_tol fitting tolerance for estimating initial values of partition sums
#' @param bound_min_fderiv_val This value used in Haley's numerical method for solving equations and represent a 
#' minimal value of denominator in first derivative of a function
#' @param bound_max_steps Maximum number of iterations for algorithm which estimates boundary values for partition sums
#' @param run_priorEM If TRUE the function runs Expectation Maximization algorithm for fitting prior probabilities of binding objects
#' @param priorEM_fit_tol Fitting tolerance for running prior EM 
#' @param priorEM_max_steps Maximum number of iterations in prior EM
#'
#' @return
#' @export
#'
#' @examples
rnome <- function(data, binding_models, bgprotectprob, bgprior, bound_fit_tol, bound_min_fderiv_val, bound_max_steps, run_priorEM, priorEM_fit_tol, priorEM_max_steps) {
  out.list <- run_cpp_rnome(data, binding_models, bgprotectprob, bgprior, bound_fit_tol, bound_min_fderiv_val, bound_max_steps, run_priorEM, priorEM_fit_tol, priorEM_max_steps)
  
  lapply(out.list,as.data.frame)
  
}


#' Preprocesses the data for rnome
#'
#' @param matr A matrix containing NOMe-Seq data. Rows represent reads (or fragments), 
#' columns represent position in an amplicon. NOTE that \code{data} must contain only positions of selected GpCs where '0' represent unprotected "C", 
#' '1' represent protected "C" and all other positions must be filled with "NA".
#' @return
#' @export
#'
#' @examples
create_data_list <- function(matr){
  stopifnot(is.matrix(matr))
  if(is.null(row.names(matr))){
    row.names(matr) <- 1:nrow(matr)
  }
  
  dat.out <- lapply(row.names(matr),function(nm){
    seq <- matr[nm,]
    seq[is.na(seq)] <- 2
    seq <- paste(seq,collapse = "")
    list("DATA" = seq,
         "NAME" = nm)
  })
  dat.out
}