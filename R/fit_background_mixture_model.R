
get_mixmodel_init_values <- function(nchains,
                                     model = c("mixture","betabinom"),
                                     mixmodel_indata){
  
  model <- match.arg(model)
  
  if(model == "mixture"){
    init_val <- lapply(1:nchains,
                       function(chain){
                         list(theta = 0.99,
                              # bg_protect_prob = stats::rbeta(1,
                              #                         mixmodel_indata$alpha_distrbetaprior_bg_protect,
                              #                         mixmodel_indata$beta_distrbetaprior_bg_protect),
                              bg_protect_prob = 0.01,
                              # alpha_betabinom = stats::runif(1,0,5),
                              # beta_betabinom = stats::runif(1,0,5)
                              alpha_betabinom = 1,
                              beta_betabinom = 1
                              
                              )
                       })
  } else if(model == "betabinom"){
    init_val <- lapply(1:nchains,
                       function(chain){
                         list(alpha_betabinom = stats::runif(1,1e-20,1000),
                              beta_betabinom = stats::runif(1,1e-20,1000))
                       })
  }
}



#' Fit background model
#' 
#' 
#' Using NOMe-seq data for lambda DNA this function fits mixture model which has binomial and beta-binomial components.
#' After the model is fit parameters for background prior distribution is constructed snd returned,
#' namely parameters \code{alpha} and \code{beta} for beta prior distribution for parameter \code{bg_protect_prob} for the footprint inference.  
#'
#'
#' @param lambda_DNA_nome_data NOMe-Seq data for lambda DNA. 
#' lambda_DNA_nome_data must be either \code{matrix} or \code{list}.
#' @param method method for parameter inference that calls corresponding \code{rstan} function. 
#' \code{optimizing} returns point estimate corresponding to maximum in posterior denstiy.
#' \code{sampling} returns Hamiltonian Monte Carlo draws from posterior distribution.
#' \code{vb} return draws from ADVI approximation of posterior distribution (not tested yet).
#' @param alpha_distrbetaprior_mixcoeff alpha (shape1 in \code{rbeta}) parameter for beta prior distribution for mixing coefficient.
#' @param beta_distrbetaprior_mixcoeff beta (shape2 in \code{rbeta}) parameter for beta prior distribution for mixing coefficient.
#' @param alpha_distrbetaprior_bg_protect alpha (shape1 in \code{rbeta}) parameter for beta prior distribution for \code{bg_protect_prob}.
#' @param beta_distrbetaprior_bg_protect beta (shape2 in \code{rbeta}) parameter for beta prior distribution for \code{bg_protect_prob}.
#' @param nchains number of markov chains for \code{rstan::sampling}. Currently has no effect
#' @param ncpu number of cpu to use for \code{rstan::sampling}. Currently has no effect
#' @param ... parameters for \code{rstan::optimizing}

#'
#' @return a list containing following slots:
#' \code{rstan_output_object} \code{rstan} objected returned by function \code{optimizing}.
#' \code{background_model_params} list containing parameters for function which performs footprint inference, i.e. \code{nomeR::infer_footprints_stan_sampling}.
#' The list \code{background_model_params} contains 3 values:
#' \code{bg_protect_prob_fixed} inferred value for background probability to emit protected base, i.e. \code{bg_protect_prob}.
#' bg_protect_prior_alpha alpha (shape1 in \code{rbeta}) parameter for beta distribution which will be used as a prior distribution for \code{bg_protect_prob} in function for footprint inference.
#' bg_protect_prior_beta beta (shape2 in \code{rbeta}) parameter for beta distribution which will be used as a prior distribution for \code{bg_protect_prob} in function for footprint inference.
#' 
#' 
#' @export
#'
fit_background_mixture_model <- function(lambda_DNA_nome_data,
                                         method = c("optimizing",
                                                    "sampling",
                                                    "vb"),
                                         ## prior distribution of mixture coefficient is beta. below are params for it
                                         alpha_distrbetaprior_mixcoeff = 1,
                                         beta_distrbetaprior_mixcoeff = 0.001,
                                         
                                         ## prior distribution of bg_protect_prob is beta. below are params for it
                                         alpha_distrbetaprior_bg_protect = 0.1,
                                         beta_distrbetaprior_bg_protect = 1,
                                         nchains = 4,
                                         ncpu = nchains,
                                         ...) {
  method <- match.arg(method)
  
  if(method != "optimizing"){
    stop(paste0("Currently only method=\"optimizing\" is supported! Please change method=\"",method,"\" to method=\"optimizing\""))
  }
  
  ## depending on the type of lambda_DNA_nome_data get aggregate statistic
  if(is.matrix(lambda_DNA_nome_data)){
    ## count total number of nonNA positions
    n_total_pos <- rowSums(!is.na(lambda_DNA_nome_data))
    ## total protected counts
    n_protect_pos <- rowSums(lambda_DNA_nome_data,na.rm = T)

  } else if(is.list(lambda_DNA_nome_data)){
    
    ## count total number of nonNA positions
    n_total_pos <- sapply(lambda_DNA_nome_data,
                          function(x){
                            sum(!is.na(x))
                          },simplify = T)
    
    ## total protected counts
    n_protect_pos <- sapply(lambda_DNA_nome_data,
                            function(x){
                              sum(x,na.rm = T)
                            },simplify = T)
  } else {
    stop("lambda_DNA_nome_data must be matrix of list of vectors")
  }
  
  
  
  ## aggregate stats
  aggr.stats <- as.data.frame(table(n_total_pos,n_protect_pos),
                              stringsAsFactors = F)
  aggr.stats <- aggr.stats[aggr.stats$Freq > 0,]
  
  mixture_model_input_data <- list(n_table_entries = nrow(aggr.stats),
                                   n_total_pos = as.integer(as.character(aggr.stats$n_total_pos)),
                                   n_protect_pos = as.integer(as.character(aggr.stats$n_protect_pos)),
                                   n_frequencies = as.integer(aggr.stats$Freq),
                                   
                                   ## prior distribution of mixture coefficient is beta. below are params for it
                                   alpha_distrbetaprior_mixcoeff = alpha_distrbetaprior_mixcoeff,
                                   beta_distrbetaprior_mixcoeff = beta_distrbetaprior_mixcoeff,
                                   
                                   ## prior distribution of bg_protect_prob is beta. below are params for it
                                   alpha_distrbetaprior_bg_protect = alpha_distrbetaprior_bg_protect,
                                   beta_distrbetaprior_bg_protect = beta_distrbetaprior_bg_protect
  )
  
  ## get initial values for model
  mixmodel_init_vals <- get_mixmodel_init_values(nchains = nchains,
                                                 model = "mixture",
                                                 mixmodel_indata = mixture_model_input_data)
  
  if(method == "optimizing"){
    
    mixturemodel_opt_params <- rstan::optimizing(object = stanmodels[["background_mixture_model"]],
                                                 data = mixture_model_input_data,
                                                 init = mixmodel_init_vals[[1]],
                                                 ...)
    
    ## extract fitted parameters
    theta <- mixturemodel_opt_params$par["theta"]
    bg_protect_prob <- mixturemodel_opt_params$par["bg_protect_prob"]
    alpha_betabinom <- mixturemodel_opt_params$par["alpha_betabinom"]
    beta_betabinom <- mixturemodel_opt_params$par["beta_betabinom"]
    rstan_output_object <- mixturemodel_opt_params
    
    ## warn if theta is too low or bg_protect_prob too high
    if(theta <= 0.1 | bg_protect_prob >= 0.4){
      warning(paste0("WARNING: MAP estimates for theta and bg_protect_prob seem unrealistic!\n",
                     "theta=",theta,"; expect >= 0.5\n",
                     "bg_protect_prob=",bg_protect_prob,"; expect <= 0.4\n",
                     "Please check rstan_output_object$return_code in output object. If it is not 0 the algorithm failed!"))
    }
    
    
  } else if(method == "sampling") {
    mixturemodel_sampling_params <- rstan::sampling(object = stanmodels[["background_mixture_model"]],
                                                    data = mixture_model_input_data,
                                                    chains = nchains,
                                                    init = mixmodel_init_vals,
                                                    cores = ncpu,
                                                    ...)
    
  } else if(method == "vb"){
    mixturemodel_vb_params <- rstan::vb(object = stanmodels[["background_mixture_model"]],
                                                     data = mixture_model_input_data,
                                                     init = mixmodel_init_vals[[1]],
                                                     ...)
  } else {
    stop("method must be either of \"optimizing\", \"sampling\", \"vb\"")
  }
  
  ### using fitted parameters for the mixture model calculate posterior of belonging to binom or betabinom model
  binom_component <- theta * stats::dbinom(x = mixture_model_input_data[["n_protect_pos"]],
                                           size = mixture_model_input_data[["n_total_pos"]],
                                           prob = bg_protect_prob)
  
  betabinom_component <- (1-theta) * extraDistr::dbbinom(x = mixture_model_input_data[["n_protect_pos"]],
                                                         size = mixture_model_input_data[["n_total_pos"]],
                                                         alpha = alpha_betabinom,
                                                         beta = beta_betabinom)
  
  posterior_binom <- binom_component/(binom_component + betabinom_component)
  
  ## using posterior calculate parameters for beta prior for footprint inference
  bg_protect_prior_alpha <- sum(posterior_binom * mixture_model_input_data[["n_protect_pos"]] * mixture_model_input_data[["n_frequencies"]])
  bg_protect_prior_beta <- sum(posterior_binom * mixture_model_input_data[["n_total_pos"]] * mixture_model_input_data[["n_frequencies"]]) - 
    bg_protect_prior_alpha
  
  bg_protect_prior_alpha <- bg_protect_prior_alpha + alpha_distrbetaprior_bg_protect
  bg_protect_prior_beta <- bg_protect_prior_beta + beta_distrbetaprior_bg_protect
  
  bg_protect_prob_fixed <- bg_protect_prior_alpha/(bg_protect_prior_alpha + bg_protect_prior_beta)
  
  return(list("rstan_output_object" = rstan_output_object,
              "background_model_params" = list("bg_protect_prob_fixed" = bg_protect_prob_fixed, 
                                               "bg_protect_prior_alpha" = bg_protect_prior_alpha,
                                               "bg_protect_prior_beta" = bg_protect_prior_beta)))
}