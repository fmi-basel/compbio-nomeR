## validate input and construct stan data

validate_construct_stan_input <- function(joint_freq_table,
                                          footprint_prior_diralphas,
                                          
                                          ftp_background_model,
                                          background_model_params,
                                          
                                          ftp_model_params){
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
  
  ## check ftp_model_params and ftp_model_params
  if(ftp_background_model == "informative_prior"){
    stan_model_name <- "ftp_inference_background_informative_prior"
    
    ## check bg model parameters 
    if(background_model_params[["bg_protect_min"]] < 0 | background_model_params[["bg_protect_min"]] > 1 |
       background_model_params[["bg_protect_max"]] < 0 | background_model_params[["bg_protect_max"]] > 1 |
       background_model_params[["bg_protect_min"]] > background_model_params[["bg_protect_max"]] | 
       background_model_params[["bg_protect_alpha"]] <= 0 | background_model_params[["bg_protect_beta"]] <= 0){
      
      stop(paste0("Incorrect values in background_model_params for background model parameters.\n",
                  "Expect background_model_params[[\"bg_protect_min\"]] < background_model_params[[\"bg_protect_max\"]]",
                  "background_model_params[[\"bg_protect_min\"]] = ",background_model_params[["bg_protect_min"]],". Expect to be in [0,1].\n",
                  "background_model_params[[\"bg_protect_max\"]] = ",background_model_params[["bg_protect_max"]],". Expect to be in [0,1].\n",
                  "background_model_params[[\"bg_protect_alpha\"]] = ",background_model_params[["bg_protect_alpha"]],". Expect to be > 0.\n",
                  "background_model_params[[\"bg_protect_beta\"]] = ",background_model_params[["bg_protect_beta"]],". Expect to be > 0.\n",
      ))
      
    }
    
    ## check ftp model params
    if(ftp_model_params[["ftp_protect_min"]] < 0 | ftp_model_params[["ftp_protect_min"]] > 1 |
       ftp_model_params[["ftp_protect_max"]] < 0 | ftp_model_params[["ftp_protect_max"]] > 1 |
       ftp_model_params[["ftp_protect_min"]] > ftp_model_params[["ftp_protect_max"]] | 
       ftp_model_params[["ftp_protect_alpha"]] <= 0 | ftp_model_params[["ftp_protect_beta"]] <= 0){
      
      stop(paste0("Incorrect values in ftp_model_params for footprint parameters.\n",
                  "Expect ftp_model_params[[\"ftp_protect_min\"]] < ftp_model_params[[\"ftp_protect_max\"]]",
                  "ftp_model_params[[\"ftp_protect_min\"]] = ",ftp_model_params[["ftp_protect_min"]],". Expect to be in [0,1].\n",
                  "ftp_model_params[[\"ftp_protect_max\"]] = ",ftp_model_params[["ftp_protect_max"]],". Expect to be in [0,1].\n",
                  "ftp_model_params[[\"ftp_protect_alpha\"]] = ",ftp_model_params[["ftp_protect_alpha"]],". Expect to be > 0.\n",
                  "ftp_model_params[[\"ftp_protect_beta\"]] = ",ftp_model_params[["ftp_protect_beta"]],". Expect to be > 0.\n",
      ))
      
    }
    
    
    ## create input data list
    stan_inputdata <- list("n_ftp" = length(footprint_prior_diralphas),
                           "ftp_prior_alphas" = footprint_prior_diralphas,
                           "n_spac" = nrow(joint_freq_table),
                           "spacings" = joint_freq_table[,"S"],
                           "spacing_counts" = joint_freq_table[,c("N00","N01","N10","N11")],
                           
                           "bg_protect_min" = background_model_params[["bg_protect_min"]],
                           "bg_protect_max" = background_model_params[["bg_protect_max"]],
                           "bg_protect_alpha" = background_model_params[["bg_protect_alpha"]],
                           "bg_protect_beta" = background_model_params[["bg_protect_beta"]],
                           
                           "ftp_protect_min" = ftp_model_params[["ftp_protect_min"]],
                           "ftp_protect_max" = ftp_model_params[["ftp_protect_max"]],
                           "ftp_protect_alpha" = ftp_model_params[["ftp_protect_alpha"]],
                           "ftp_protect_beta" = ftp_model_params[["ftp_protect_beta"]])
    
    
    
  } else if(ftp_background_model == "bg_fixed"){
    stan_model_name <- "ftp_inference_background_fixed"
    
    ## check bg parameter
    if(background_model_params[["bg_protect_prob_fixed"]] < 0 | background_model_params[["bg_protect_prob_fixed"]] > 1){
      stop(paste0("bg_protect_prob_fixed in list background_model_params must be in [0,1]!\n",
                  "Current value background_model_params[[\"bg_protect_prob_fixed\"]] = ",background_model_params[["bg_protect_prob_fixed"]]))
    }
    
    ## check ftp model params
    if(ftp_model_params[["ftp_protect_min"]] < 0 | ftp_model_params[["ftp_protect_min"]] > 1 |
       ftp_model_params[["ftp_protect_max"]] < 0 | ftp_model_params[["ftp_protect_max"]] > 1 |
       ftp_model_params[["ftp_protect_min"]] > ftp_model_params[["ftp_protect_max"]] | 
       ftp_model_params[["ftp_protect_alpha"]] <= 0 | ftp_model_params[["ftp_protect_beta"]] <= 0){
      
      stop(paste0("Incorrect values in ftp_model_params for footprint parameters.\n",
                  "Expect ftp_model_params[[\"ftp_protect_min\"]] < ftp_model_params[[\"ftp_protect_max\"]]",
                  "ftp_model_params[[\"ftp_protect_min\"]] = ",ftp_model_params[["ftp_protect_min"]],". Expect to be in [0,1].\n",
                  "ftp_model_params[[\"ftp_protect_max\"]] = ",ftp_model_params[["ftp_protect_max"]],". Expect to be in [0,1].\n",
                  "ftp_model_params[[\"ftp_protect_alpha\"]] = ",ftp_model_params[["ftp_protect_alpha"]],". Expect to be > 0.\n",
                  "ftp_model_params[[\"ftp_protect_beta\"]] = ",ftp_model_params[["ftp_protect_beta"]],". Expect to be > 0.\n",
      ))
      
    }
    
    ## create input data list
    stan_inputdata <- list("n_ftp" = length(footprint_prior_diralphas),
                           "ftp_prior_alphas" = footprint_prior_diralphas,
                           "n_spac" = nrow(joint_freq_table),
                           "spacings" = joint_freq_table[,"S"],
                           "spacing_counts" = joint_freq_table[,c("N00","N01","N10","N11")],
                           
                           "bg_protect_prob" = background_model_params[["bg_protect_prob_fixed"]],
                           
                           "ftp_protect_min" = ftp_model_params[["ftp_protect_min"]],
                           "ftp_protect_max" = ftp_model_params[["ftp_protect_max"]],
                           "ftp_protect_alpha" = ftp_model_params[["ftp_protect_alpha"]],
                           "ftp_protect_beta" = ftp_model_params[["ftp_protect_beta"]])
  } else if(ftp_background_model == "ftp_bg_fixed"){
    stan_model_name <- "ftp_inference_background_footprint_fixed"
    
    ## check bg parameter
    if(background_model_params[["bg_protect_prob_fixed"]] < 0 | background_model_params[["bg_protect_prob_fixed"]] > 1){
      stop(paste0("bg_protect_prob_fixed in list background_model_params must be in [0,1]!\n",
                  "Current value background_model_params[[\"bg_protect_prob_fixed\"]] = ",background_model_params[["bg_protect_prob_fixed"]]))
    }
    
    ## check ftp model parameter
    if(ftp_model_params[["ftp_protect_prob_fixed"]] < 0 | ftp_model_params[["ftp_protect_prob_fixed"]] > 1){
      stop(paste0("ftp_protect_prob_fixed in list ftp_model_params must be in [0,1]!\n",
                  "Current value background_model_params[[\"bg_protect_prob_fixed\"]] = ",ftp_model_params[["ftp_protect_prob_fixed"]]))
    }
    
    ## create input data list
    stan_inputdata <- list("n_ftp" = length(footprint_prior_diralphas),
                           "ftp_prior_alphas" = footprint_prior_diralphas,
                           "n_spac" = nrow(joint_freq_table),
                           "spacings" = joint_freq_table[,"S"],
                           "spacing_counts" = joint_freq_table[,c("N00","N01","N10","N11")],
                           
                           "bg_protect_prob" = background_model_params[["bg_protect_prob_fixed"]],
                           "footprint_protect_prob" = ftp_model_params[["ftp_protect_prob_fixed"]],
                           
                           "ftp_protect_min" = ftp_model_params[["ftp_protect_min"]],
                           "ftp_protect_max" = ftp_model_params[["ftp_protect_max"]],
                           "ftp_protect_alpha" = ftp_model_params[["ftp_protect_alpha"]],
                           "ftp_protect_beta" = ftp_model_params[["ftp_protect_beta"]])
    
    
  } else {
    stop("Unknown ftp_background_model!")
  }
  
  return(list("stan_model_name" = stan_model_name,
              "stan_inputdata" = stan_inputdata))
  
}