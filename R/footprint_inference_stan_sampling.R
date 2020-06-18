
#' Bayesian inference of footprint abundance using STAN HMC algorithm
#'
#' This function uses Hamiltonian Monte Carlo and NUTS algorithms from rstan package for bayesian inference of footprint abundance using sampling from posterior distribution.
#'
#' @param joint_freq_table \code{data.frame} which must contain columns \code{S}, \code{N00},\code{N01},\code{N10},\code{N11}, 
#' where \code{S} represent spacings and \code{N00},\code{N01},\code{N10},\code{N11} are observed counts for 00, 01, 10 and 11 at spacing \code{S}.
#' This table can be obtained using function nomeR::count_joint_frequencies(...)
#' @param footprint_prior_diralphas vector containing parameters for dirichlet distribution which is used to model footprint abundances.
#' Default vector assumes that around 50\% of data is free of any footprints and all other footprint up-to 200bp have uninformative prior.
#' 
#' @param bg_protect_prob expected probability to observe 1 at a free position. In other words it is emission probability of background for 1.
#' @param ftp_protect_min minimum allowed value for footprint emission probability for 1.
#' @param ftp_protect_max maximum allowed value for footprint emission probability for 1.
#' @param ftp_protect_mean expected mean for footprint emission probability for 1.
#' @param ftp_protect_var expected variance for footprint emission probability for 1.
#' @param nchains number of markov chains to run.
#' @param ncpu number of CPUs to use.
#' @param rstan_control a list of tuning parameters for stan sampler algorithm.
#' @param ... parameters for `rstan::sampling` function. Please check \code{?rstan::sampling}.
#' 
#' Emission probability for a footprint to observe 1 is one parameters 
#' apart from footprint abundances which will be inferred from the data.
#' It is modelled using truncated beta distribution which is parametrized by mean (ftp_protect_mean) and variance (ftp_protect_var).
#' Inequality ftp_protect_var < ftp_protect_mean * (1 - ftp_protect_mean) must hold for correct parameterization of beta distribution.
#'
#' @return functions returns S4 class \code{rstan::stanfit} representing the fitted results. Pleae check \code{?rstan::stanfit}. 
#' @export
#'
#' @examples
footprint_inference_stan_sampling <- function(joint_freq_table,
                                              footprint_prior_diralphas = c(200,rep(1,199)),
                                              
                                              bg_protect_prob = 0.01,
                                              
                                              ftp_protect_min = 0.6,
                                              ftp_protect_max = 0.9999,
                                              ftp_protect_mean = 0.95,
                                              ftp_protect_var = 0.01,
                                              nchains = 4,
                                              ncpu = 1L,
                                              rstan_control = list(max_treedepth = 15,
                                                                  adapt_delta = 0.9),
                                              ...
                                              ){
  
  ## check if required columns in the table
  if(!all(c("S","N00","N01","N10","N11") %in% colnames(joint_freq_table))){
    stop("Incorrect joint_freq_table! joint_freq_table must contain columns S, N00, N01, N10, N11")
  }
  
  ## remove rows without data
  joint_freq_table <- joint_freq_table[rowSums(joint_freq_table[,c("N00","N01","N10","N11")]) > 0,,drop=F]
  
  ## check if anything left
  if(nrow(joint_freq_table) == 0){
    stop("Incorrect joint_freq_table! joint_freq_table does not contain non-zero counts")
  }
  ## check that there is no duplicated S
  if(any(duplicated(joint_freq_table[,"S"]))){
    stop("Incorrect joint_freq_table! joint_freq_table contains duplicated entries in S column")
  }
  
  ## check if S=1 is in the table
  if(all(joint_freq_table[,"S"] != 1)){
    stop("Incorrect joint_freq_table! Could not find a row with S = 1 which should represent total number of observed 0s and 1s in the data.")
  }
  
  ## sort table
  joint_freq_table <- joint_freq_table[order(joint_freq_table[,"S"]),]
  
  ## check footprint_prior_diralphas
  if(length(footprint_prior_diralphas) == 0){
    stop("footprint_prior_diralphas seems to be empty")
  }
  if(any(footprint_prior_diralphas <= 0)){
    stop("Incorrect footprint_prior_diralphas! non-positive values are prohibited in footprint_prior_diralphas.")
  }
  
  
  ## create data list
  stan_inputdata <- list("n_ftp" = length(footprint_prior_diralphas),
                         "ftp_prior_alphas" = footprint_prior_diralphas,
                         "n_spac" = nrow(joint_freq_table),
                         "spacings" = joint_freq_table[,"S"],
                         "spacing_counts" = joint_freq_table[,c("N00","N01","N10","N11")],
                         
                         "bg_protect_prob" = bg_protect_prob,
                         
                         "ftp_protect_min" = ftp_protect_min,
                         "ftp_protect_max" = ftp_protect_max,
                         "ftp_protect_mean" = ftp_protect_mean,
                         "ftp_protect_var" = ftp_protect_var)
  
  ## get initial values for sampling
  stan_initvals <- init_ftp_model_params_for_hmc(dir_alpha = footprint_prior_diralphas,
                                                 nchains = nchains,
                                                 ftp_protect_min = ftp_protect_min,
                                                 ftp_protect_max = ftp_protect_max,
                                                 ftp_protect_mean = ftp_protect_mean,
                                                 ftp_protect_var = ftp_protect_var)
  
  
  ## hmc sampling
  stanfit_out <- rstan::sampling(stanmodels$footprint_inference_model_v1,
                                 data = stan_inputdata,
                                 chains = nchains,
                                 init = stan_initvals,
                                 cores = ncpu,
                                 ...)
  
  return(stanfit_out)
}








