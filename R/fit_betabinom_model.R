#' Estimate parameters for beta-binomial distribution
#' 
#' This function takes protection statistics as an input (see description for \code{protect_stat_table}) and infers parameters for beta-binomial distribution which fits the data.
#'
#' @param protect_stat_table \code{matrix} with columns TotalGCH - total number of GCHs in a fragment, ProtectedGCH - number of protected GCHs in a fragments, 
#' Nfrags - number of fragments with combination of TotalGCH and ProtectedGCH.
#' @param method method for parameter inference that calls corresponding \code{\link[rstan]{rstan}} function. 
#' \code{optimizing} returns point estimate corresponding to maximum in posterior denstiy.
#' \code{sampling} returns Hamiltonian Monte Carlo draws from posterior distribution.
#' \code{vb} return draws from ADVI approximation of posterior distribution (not tested yet).
#' @param nchains number of markov chains for \code{\link[rstan]{sampling}}. 
#' @param ncpu number of cpu to use for \code{\link[rstan]{sampling}}. 
#' @param ... parameters for \code{\link[rstan]{optimizing}}

#'
#' @return a list containing following slots:
#' \code{rstan_output_object} \code{\link[rstan]{stanfit}} object returned by function \code{\link[rstan]{sampling}} or \code{\link[rstan]{vb}}, or \code{list} returned by \code{\link[rstan]{optimizing}}.
#' \code{betabinom_params} list containing infered parameters for beta-binomial distribution.
#' The list \code{betabinom_params} contains 2 values:
#' \describe{
#' \item{alpha}{alpha (shape1 in \code{\link[stats]{rbeta}}) parameter for beta distribution which can be used as parameter \code{bg_protect_alpha} for prior distribution of background noise or parameter \code{ftp_protect_alpha} for prior distribution of footprint noise}
#' \item{beta}{beta (shape2 in \code{\link[stats]{rbeta}}) parameter for beta distribution which can be used as parameter \code{bg_protect_beta} for prior distribution of background noise or parameter \code{ftp_protect_beta} for prior distribution of footprint noise}
#' }
#' 
#' 
#' @export
#'
fit_betabinom_model <- function(protect_stat_table,
                                method = c("optimizing",
                                           "sampling",
                                           "vb"),
                                nchains = 4,
                                ncpu = nchains,
                                ...) {
  method <- match.arg(method)
  
  ## check if required columns in the input matrix are present
  if(any(!c("TotalGCH","ProtectedGCH","Nfrags") %in% colnames(protect_stat_table)))
    stop("protect_stat_table must be a matrix with columns \"TotalGCH\",\"ProtectedGCH\",\"Nfrags\"")
  
  
  betabinom_model_input_data <- list(n_table_entries = nrow(protect_stat_table),
                                     n_total_pos = protect_stat_table[,"TotalGCH"],
                                     n_protect_pos = protect_stat_table[,"ProtectedGCH"],
                                     n_frequencies = protect_stat_table[,"Nfrags"])
  
  ## get random initial values which are not too far off to prevent algorithms fail
  bbinit <- get_mixmodel_init_values(nchains = nchains,model = "betabinom")
  
  
  
  if(method == "optimizing"){
    betabinom_opt_params <- rstan::optimizing(object = stanmodels[["background_betabinom_model"]],
                                              data = betabinom_model_input_data,
                                              #init="random",
                                              init = bbinit[[1]],
                                              ...)
    
    ## extract fitted parameters
    alpha_betabinom <- betabinom_opt_params$par["alpha_betabinom"]
    beta_betabinom <- betabinom_opt_params$par["beta_betabinom"]
    mean_beta <- betabinom_opt_params$par["mean_beta"]
    var_beta <- betabinom_opt_params$par["var_beta"]
    
    
    rstan_output_object <- betabinom_opt_params
    
  } else if(method == "sampling") {
    betabinom_sampling_params <- rstan::sampling(object = stanmodels[["background_betabinom_model"]],
                                                 data = betabinom_model_input_data,
                                                 chains = nchains,
                                                 #init="random",
                                                 init = bbinit,
                                                 cores = ncpu,
                                                 ...)
    ## extract fitted parameters
    sampl_summ <- rstan::summary(betabinom_sampling_params)$summary
    alpha_betabinom <- sampl_summ["alpha_betabinom","mean"]
    beta_betabinom <- sampl_summ["beta_betabinom","mean"]
    mean_beta <- sampl_summ["mean_beta","mean"]
    var_beta <- sampl_summ["var_beta","mean"]
    
    rstan_output_object <- betabinom_sampling_params
    
  } else if(method == "vb"){
    betabinom_vb_params <- rstan::vb(object = stanmodels[["background_betabinom_model"]],
                                     data = betabinom_model_input_data,
                                     #init="random",
                                     init = bbinit[[1]],
                                     ...)
    ## extract fitted parameters
    sampl_summ <- rstan::summary(betabinom_vb_params)$summary
    alpha_betabinom <- sampl_summ["alpha_betabinom","mean"]
    beta_betabinom <- sampl_summ["beta_betabinom","mean"]
    mean_beta <- sampl_summ["mean_beta","mean"]
    var_beta <- sampl_summ["var_beta","mean"]
    
    rstan_output_object <- betabinom_vb_params
    
    
  }
  
  return(list("rstan_output_object" = rstan_output_object,
              "betabinom_params" = list("alpha" = alpha_betabinom,
                                        "beta" = beta_betabinom,
                                        "mean_beta" = mean_beta,
                                        "var_beta" = var_beta)))
}