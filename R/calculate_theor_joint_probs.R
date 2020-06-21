


#' Function for calculating theoretical joint probabilities at certain distance
#' 
#' Using assumption that all non-background footprints have the same emission probabilities we calculate 
#' theoretical probabilities 
#' \code{P(0,0|S)}, \code{P(0,1|S)}, \code{P(1,0|S)}, \code{P(1,1|S)} for all distances \code{S} up to \code{max_spacing}.
#' \code{P(0,0|S)} etc, are probabilities to find \code{0} and \code{0} separated by distance \code{S}.
#' @param ftp_cover_priors coverage probabilities for each footprint. 
#' IMPORTANT! This must be a vector where \code{ftp_cover_priors[1]} is coverage probability of background and we assume that indexes in the vector correspond to the footprint lengths.
#' @param bg_protect_prob emission probability for background for \code{1}.
#' @param footprint_protect_prob emission probability for footprints for \code{1}.
#' @param max_spacing maximum distance
#'
#' @return \code{data.frame} containing columns: 
#' \code{S} - distance where \code{S = 1} corresponds to the same position and probabilities are marginal probabilities to observe 0 and 1; 
#' \code{P00}, \code{P01}, \code{P10}, \code{P11} - calculated theoretical \code{P(0,0|S)}, \code{P(0,1|S)}, \code{P(1,0|S)}, \code{P(1,1|S)}.
#' @export
#'
calculate_theor_joint_probs <- function(ftp_cover_priors, # here vector of priors also represent lengths, namely ith element of the vector has length i+1, e.g. ftp_cover_priors[0] is a prior for bg with length 1. Make sure that R function passes correct vector with priors
                                  bg_protect_prob,
                                  footprint_protect_prob,
                                  max_spacing){
  
  # if(any(ftp_cover_priors < 0) | abs(sum(ftp_cover_priors) - 1) > .Machine$double.eps){
  #   stop("ftp_cover_priors must be vector of non-negative numbers which sum up to 1")
  # }
  
  theor_joint_prob <- calculate_theor_joint_prob_cpp(ftp_cover_priors,
                                                     bg_protect_prob,
                                                     footprint_protect_prob,
                                                     max_spacing)
  
  return(as.data.frame(theor_joint_prob))
  
}


