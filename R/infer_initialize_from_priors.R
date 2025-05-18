## initialization of parameter values by drawing from prior distribution
#' @keywords internal
#' @noRd
.init_param_from_prior_distr <- function(stan_input,
                                         nchains,
                                         delta_from_max_min =  0.001) {
    
    stan_model_name <- stan_input[["stan_model_name"]]
    stan_inputdata <- stan_input[["stan_inputdata"]]
    ## calculate didchlet alphas
    dirich_alpha <- stan_inputdata[["ftp_prior_cover"]] * 
        stan_inputdata[["total_cnt_prior_dirich"]]
    
    init_vals <- lapply(seq_len(nchains), function(x) {
        ftp_cover_probs <- as.vector(extraDistr::rdirichlet(1,dirich_alpha))
        ## substitute 0s if any by small number because rstan fails because of log(0)
        ftp_cover_probs[ftp_cover_probs == 0] <- .Machine$double.eps / 2
        ftp_cover_probs <- ftp_cover_probs/sum(ftp_cover_probs)
        
        if (stan_model_name == "ftp_inference_informative_prior") {
            bg_protect_alpha <- stan_inputdata[["bg_protect_mean"]] *
                stan_inputdata[["bg_protect_totcount"]]
            bg_protect_beta <- stan_inputdata[["bg_protect_totcount"]] * 
                (1 - stan_inputdata[["bg_protect_mean"]])
            ftp_protect_alpha <- stan_inputdata[["ftp_protect_mean"]] *
                stan_inputdata[["ftp_protect_totcount"]]
            ftp_protect_beta <- stan_inputdata[["ftp_protect_totcount"]] * 
                (1 - stan_inputdata[["ftp_protect_mean"]])
            
            bg_protect_prob <- truncdist::rtrunc(
                1, spec = "beta",
                a = stan_inputdata[["bg_protect_min"]] + delta_from_max_min,
                b = stan_inputdata[["bg_protect_max"]] - delta_from_max_min,
                shape1 = bg_protect_alpha,
                shape2 = bg_protect_beta)
            footprint_protect_prob <- truncdist::rtrunc(
                1, spec = "beta",
                a = stan_inputdata[["ftp_protect_min"]] + delta_from_max_min,
                b = stan_inputdata[["ftp_protect_max"]] - delta_from_max_min,
                shape1 = ftp_protect_alpha,
                shape2 = ftp_protect_beta)
            chain_init_out <- list("ftp_abundances" = ftp_cover_probs,
                                   "bg_protect_prob" = bg_protect_prob,
                                   "ftp_protect_prob" = footprint_protect_prob)
        } else if (stan_model_name == "ftp_inference_bg_fixed") {
            ftp_protect_alpha <- stan_inputdata[["ftp_protect_mean"]] *
                stan_inputdata[["ftp_protect_totcount"]]
            ftp_protect_beta <- stan_inputdata[["ftp_protect_totcount"]] * 
                (1 - stan_inputdata[["ftp_protect_mean"]])
            footprint_protect_prob <- truncdist::rtrunc(
                1, spec = "beta",
                a = stan_inputdata[["ftp_protect_min"]] + delta_from_max_min,
                b = stan_inputdata[["ftp_protect_max"]] - delta_from_max_min,
                shape1 = ftp_protect_alpha,
                shape2 = ftp_protect_beta)
            
            chain_init_out <- list("ftp_abundances" = ftp_cover_probs,
                                   "ftp_protect_prob" = footprint_protect_prob)
        } else if (stan_model_name == "ftp_inference_ftp_bg_fixed") {
            chain_init_out <- list("ftp_abundances" = ftp_cover_probs)
        } else{
            stop("Unknown ftp_background_model!")
        }
        
        return(chain_init_out)
    })
    
    return(init_vals)
}
