## validate and contsruct input for inference
.validate_construct_stan_input <- function(cooc_ctable,
                                           ftp_lengths,
                                           bg_prior_cover,
                                           ftp_prior_cover,
                                           total_cnt_prior_dirich,
                                           ftp_bg_model,
                                           bg_model_params,
                                           ftp_model_params
                                           ){
  #### CHECK ARGUMENTS ####
  coll = checkmate::makeAssertCollection()
  
  ## validate count table
  if(checkmate::check_data_frame(x=cooc_ctable,types = "integerish",min.rows=1)){
    ## check colnames
    checkmate::assert_names(x=colnames(cooc_ctable),
                            must.include = c("S","N00","N01","N10","N11"),
                            add=coll)
    
    
    ## remove rows without data
    cooc_ctable <- cooc_ctable[rowSums(cooc_ctable[,c("N00","N01","N10","N11")]) > 0,,drop=F]
    
    ## check if anything left
    checkmate::assert_data_frame(x=cooc_ctable,
                                 types = "integerish",
                                 min.rows=1,
                                 any.missing = FALSE,
                                 add = coll)
    
    
    ## check vector S
    checkmate::assert_integerish(x=cooc_ctable[,"S"],
                                 lower=1L,
                                 any.missing = F,
                                 min.len=1,
                                 unique = TRUE, 
                                 add=coll)
    
    ## check if S=1 is in the table
    if(all(cooc_ctable[,"S"] != 1)){
      stop("Incorrect cooc_ctable! Could not find a row with S = 1 which should represent total number of observed 0s and 1s in the data.")
    }
    
    ## sort table
    cooc_ctable <- cooc_ctable[order(cooc_ctable[,"S"]),]
  }
  
  ## validate ftp_lengths
  checkmate::assert_integerish(x = ftp_lengths,
                               lower=2,
                               any.missing = F,
                               min.len = 1,
                               unique=T,
                               add = coll)
  
  ## validate bg_prior_cover
  checkmate::assert_number(x = bg_prior_cover,
                           lower=0,upper=1,
                           add=coll)
  
  ## validate ftp_prior_cover
  if(is.null(ftp_prior_cover)){
    ## if ftp_prior_cover is not provided use non-informative prior
    ftp_prior_cover <- rep((1 - bg_prior_cover)/length(ftp_lengths),length(ftp_lengths))
  }
  checkmate::assert_double(x = ftp_prior_cover, lower = 0,upper=1,finite=TRUE,any.missing = F,len = length(ftp_lengths),add = coll)
  #browser()
  ## validate total_cnt_prior_dirich
  if(is.null(total_cnt_prior_dirich)){
    ## if total_cnt_prior_dirich is not provided choose such that footprint dirichlet alpha are all 1, i.e. non-informative dirichlet prior
    total_cnt_prior_dirich <- length(ftp_lengths)/(1 - bg_prior_cover)
  }
  checkmate::assert_number(x = total_cnt_prior_dirich,lower=1,add = coll)
  
  # validate bg_model_params
  if(checkmate::check_list(bg_model_params, types = c("double","integerish"),any.missing = T,all.missing = F)){
    if(ftp_bg_model == "informative_prior"){
      if(checkmate::check_names(x = names(bg_model_params),must.include = c("bg_protect_min","bg_protect_max","bg_protect_mean","bg_protect_totcount"))){
        checkmate::assert(checkmate::check_number(x=bg_model_params[["bg_protect_min"]],lower=0,upper=1),
                          checkmate::check_number(x=bg_model_params[["bg_protect_max"]],lower=0,upper=1),
                          checkmate::check_number(x=bg_model_params[["bg_protect_mean"]],lower=0,upper=1),
                          checkmate::check_number(x=bg_model_params[["bg_protect_totcount"]],lower=0),
                          add=coll,combine="and")
      }
    } else{
      if(checkmate::check_names(x = names(bg_model_params),must.include = c("bg_protect_prob_fixed"))){
        checkmate::assert_number(x=bg_model_params[["bg_protect_prob_fixed"]],lower=0,upper=1,add=coll)
      }
    }
  }
  
  
  ## validate ftp_model_params
  if(checkmate::check_list(ftp_model_params, types = c("double","integerish"),any.missing = T,all.missing = F)){
    if(ftp_bg_model %in% c("informative_prior","bg_fixed")){
      if(checkmate::check_names(x = names(ftp_model_params),must.include = c("ftp_protect_min","ftp_protect_max","ftp_protect_mean","ftp_protect_totcount"))){
        checkmate::assert(checkmate::check_number(x=ftp_model_params[["ftp_protect_min"]],lower=0,upper=1),
                          checkmate::check_number(x=ftp_model_params[["ftp_protect_max"]],lower=0,upper=1),
                          checkmate::check_number(x=ftp_model_params[["ftp_protect_mean"]],lower=0,upper=1),
                          checkmate::check_number(x=ftp_model_params[["ftp_protect_totcount"]],lower=0),
                          add=coll,combine="and")
      }
    } else{
      if(checkmate::check_names(x = names(ftp_model_params),must.include = c("ftp_protect_prob_fixed"))){
        checkmate::assert_number(x=ftp_model_params[["ftp_protect_prob_fixed"]],lower=0,upper=1,add=coll)
      }
    }
  }
  
  
  ## check if sum of bg prior and ftp is 1 and rescale if neccessary
  prior_cover <- c(bg_prior_cover,ftp_prior_cover)
  totcover <- sum(prior_cover)
  if(totcover != 1){
    warning("Summ of prior coverages bg_prior_cover and ftp_prior_cover is ",totcover," but should equal to 1. Rescale prior coverages accordingly.")
    prior_cover <- prior_cover/totcover
  }
  
  ## convert prior cover and tot counts to vector of dirichlet alphas
  ## update. in the v0.9 stan code the required data is prior_cover
  ## dirich_alpha <- prior_cover * total_cnt_prior_dirich
  ## add bg into vector of ftp lengths
  ftp_lengths <- c(1,ftp_lengths)
  
  ## calculate Beta distribution shape parameters
  ## update v0.9. the data does not require shape parameters. they are calculated in the stan code
  # if(ftp_bg_model == "informative_prior"){
  #   bg_protect_alpha <- bg_model_params[["bg_protect_mean"]] * bg_model_params[["bg_protect_totcount"]]
  #   bg_protect_beta <- bg_model_params[["bg_protect_totcount"]] * (1 - bg_model_params[["bg_protect_mean"]])
  # }
  # 
  # if(ftp_bg_model %in% c("informative_prior","bg_fixed")){
  #   ftp_protect_alpha <- ftp_model_params[["ftp_protect_mean"]] * ftp_model_params[["ftp_protect_totcount"]]
  #   ftp_protect_beta <- ftp_model_params[["ftp_protect_totcount"]] * (1 - ftp_model_params[["ftp_protect_mean"]])
  # }
  
  ## prepare input for stan
  ## v0.9 stan code data section
  # data {
  #   // input count coocurrnce table
  #   int n_spac;
  #   array[n_spac] int<lower=1> spacings;
  #   array[n_spac,4] int<lower=0> spacing_counts;
  #   
  #   // footprint hyperparameters
  #   int n_ftp;
  #   array[n_ftp] int ftp_lengths;
  #   vector[n_ftp] ftp_prior_cover; // vector of expected footprint coverages. this vector incudes 1, i.e. background and ith element represent footprint length
  #   real total_cnt_prior_dirich; // total count for dirichlet distribution
  #   
  #   
  #   // background noise hyperparameters
  #   real bg_protect_min; // minimum allowed value
  #   real bg_protect_max; // maximum allowed value
  #   real bg_protect_mean; // expected mean for bg noise
  #   real bg_protect_totcount; // total count for beta distribution
  #   
  #   // footprint noise hyperparameters
  #   real ftp_protect_min;
  #   real ftp_protect_max;
  #   real ftp_protect_mean;
  #   real ftp_protect_totcount;
  # }
  
  
  if(ftp_bg_model == "informative_prior"){
    stan_model_name <- "ftp_inference_informative_prior"
    stan_inputdata <- list("n_ftp" = length(ftp_lengths),
                           "ftp_lengths" = ftp_lengths,
                           "ftp_prior_cover" = prior_cover,
                           "total_cnt_prior_dirich" = total_cnt_prior_dirich,
                           
                           
                           "n_spac" = nrow(cooc_ctable),
                           "spacings" = cooc_ctable[,"S"],
                           "spacing_counts" = cooc_ctable[,c("N00","N01","N10","N11")],
                           
                           "bg_protect_min" = bg_model_params[["bg_protect_min"]],
                           "bg_protect_max" = bg_model_params[["bg_protect_max"]],
                           "bg_protect_mean" = bg_model_params[["bg_protect_mean"]],
                           "bg_protect_totcount" = bg_model_params[["bg_protect_totcount"]],
                           
                           "ftp_protect_min" = ftp_model_params[["ftp_protect_min"]],
                           "ftp_protect_max" = ftp_model_params[["ftp_protect_max"]],
                           "ftp_protect_mean" = ftp_model_params[["ftp_protect_mean"]],
                           "ftp_protect_totcount" = ftp_model_params[["ftp_protect_totcount"]])
  } else if(ftp_bg_model == "bg_fixed"){
    stan_model_name <- "ftp_inference_bg_fixed"
    
    stan_inputdata <- list("n_ftp" = length(ftp_lengths),
                           "ftp_lengths" = ftp_lengths,
                           "ftp_prior_cover" = prior_cover,
                           "total_cnt_prior_dirich" = total_cnt_prior_dirich,
                           
                           
                           "n_spac" = nrow(cooc_ctable),
                           "spacings" = cooc_ctable[,"S"],
                           "spacing_counts" = cooc_ctable[,c("N00","N01","N10","N11")],
                           
                           "bg_protect_prob" = bg_model_params[["bg_protect_prob_fixed"]],
                           
                           "ftp_protect_min" = ftp_model_params[["ftp_protect_min"]],
                           "ftp_protect_max" = ftp_model_params[["ftp_protect_max"]],
                           "ftp_protect_mean" = ftp_model_params[["ftp_protect_mean"]],
                           "ftp_protect_totcount" = ftp_model_params[["ftp_protect_totcount"]])
    
    
    
  } else if(ftp_bg_model == "ftp_bg_fixed"){
    stan_model_name <- "ftp_inference_ftp_bg_fixed"
    stan_inputdata <- list("n_ftp" = length(ftp_lengths),
                           "ftp_lengths" = ftp_lengths,
                           "ftp_prior_cover" = prior_cover,
                           "total_cnt_prior_dirich" = total_cnt_prior_dirich,
                           
                           
                           "n_spac" = nrow(cooc_ctable),
                           "spacings" = cooc_ctable[,"S"],
                           "spacing_counts" = cooc_ctable[,c("N00","N01","N10","N11")],
                           
                           "bg_protect_prob" = bg_model_params[["bg_protect_prob_fixed"]],
                           "ftp_protect_prob" = ftp_model_params[["ftp_protect_prob_fixed"]],
                           
                           "ftp_protect_min" = ftp_model_params[["ftp_protect_min"]],
                           "ftp_protect_max" = ftp_model_params[["ftp_protect_max"]],
                           "ftp_protect_mean" = ftp_model_params[["ftp_protect_mean"]],
                           "ftp_protect_totcount" = ftp_model_params[["ftp_protect_totcount"]])
    
  } else{
    stop("Incorrect ftp_bg_model.")
  }
  
  checkmate::reportAssertions(coll)
  return(list("stan_model_name" = stan_model_name,
              "stan_inputdata" = stan_inputdata))
  
  
}