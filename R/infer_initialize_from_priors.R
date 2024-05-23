## this function set initial values for inference from prior distributions

.infer_initialize_from_prior <- function(dir_alpha = c(200,rep(1,199)),
                                         nchains = 4,
                                         ftp_background_model = c("ftp_bg_fixed","bg_fixed","informative_prior"),
                                         background_model_params = list("bg_protect_prob_fixed" = 0.01,
                                                                        "bg_protect_min" = 0.0,
                                                                        "bg_protect_max" = 0.4,
                                                                        "bg_protect_alpha" = 675,
                                                                        "bg_protect_beta" = 57126),
                                         
                                         ftp_model_params = list("ftp_protect_prob_fixed" = 0.99,
                                                                 "ftp_protect_min" = 0.6,
                                                                 "ftp_protect_max" = 1,
                                                                 "ftp_protect_alpha" = 1,
                                                                 "ftp_protect_beta" = 0.01),
                                         delta_from_max_min = 0.001){
  ## check footprint_prior_diralphas
  if(length(dir_alpha) == 0){
    stop("dir_alpha seems to be empty")
  }
  if(any(dir_alpha <= 0)){
    stop("Incorrect dir_alpha! non-positive values are prohibited in footprint_prior_diralphas.")
  }
  
  ftp_background_model <- match.arg(ftp_background_model)
  
  init_vals <- lapply(1:nchains,function(x){
    ftp_cover_probs <- as.vector(extraDistr::rdirichlet(1,dir_alpha))
    ## substitute 0s if any by small number because rstan fails because of log(0)
    ftp_cover_probs[ftp_cover_probs == 0] <- .Machine$double.eps/2
    ftp_cover_probs <- ftp_cover_probs/sum(ftp_cover_probs)
    
    
    if(ftp_background_model == "informative_prior"){
      bg_protect_prob <- truncdist::rtrunc(1,spec="beta",
                                           a = background_model_params[["bg_protect_min"]] + delta_from_max_min,
                                           b = background_model_params[["bg_protect_max"]] - delta_from_max_min,
                                           shape1 = background_model_params[["bg_protect_alpha"]],
                                           shape2 = background_model_params[["bg_protect_beta"]])
      footprint_protect_prob <- truncdist::rtrunc(1,spec="beta",
                                                  a = ftp_model_params[["ftp_protect_min"]] + delta_from_max_min,
                                                  b = ftp_model_params[["ftp_protect_max"]] - delta_from_max_min,
                                                  shape1 = ftp_model_params[["ftp_protect_alpha"]],
                                                  shape2 = ftp_model_params[["ftp_protect_beta"]])
      chain_init_out <- list("ftp_cover_probs" = ftp_cover_probs,
                             "bg_protect_prob" = bg_protect_prob,
                             "footprint_protect_prob" = footprint_protect_prob)
    } else if(ftp_background_model == "bg_fixed"){
      footprint_protect_prob <- truncdist::rtrunc(1,spec="beta",
                                                  a = ftp_model_params[["ftp_protect_min"]] + delta_from_max_min,
                                                  b = ftp_model_params[["ftp_protect_max"]] - delta_from_max_min,
                                                  shape1 = ftp_model_params[["ftp_protect_alpha"]],
                                                  shape2 = ftp_model_params[["ftp_protect_beta"]])
      
      chain_init_out <- list("ftp_cover_probs" = ftp_cover_probs,
                             "footprint_protect_prob" = footprint_protect_prob)
    } else if(ftp_background_model == "ftp_bg_fixed"){
      chain_init_out <- list("ftp_cover_probs" = ftp_cover_probs)
    } else{
      stop("Unknown ftp_background_model!")
    }
    
    return(chain_init_out)
    
  })
  
  return(init_vals)
}