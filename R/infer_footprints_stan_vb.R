#' Black-box Variational Bayes for inference of footprint abundance using STAN Automatic Differentiation Variational Inference (ADVI) algorithm
#'
#' @param joint_freq_table \code{data.frame} which must contain columns \code{S}, \code{N00},\code{N01},\code{N10},\code{N11}, 
#' where \code{S} represent spacings and \code{N00},\code{N01},\code{N10},\code{N11} are observed counts for 00, 01, 10 and 11 at spacing \code{S}.
#' This table can be obtained using function nomeR::count_joint_frequencies(...)
#' @param footprint_prior_diralphas vector containing parameters for dirichlet distribution which is used to model footprint abundances.
#' Default vector assumes that around 50\% of data is free of any footprints and all other footprint up-to 200bp have uninformative prior.
#' 
#' @param background_model type of background model.
#' "fixed" background model assumes model parameter "bg_protect_prob" to be constant and equal to slot `bg_protect_prob_fixed` in `background_model_params`.
#' 
#' "informative_prior" background model runs inference for model parameter "bg_protect_prob" using strongly informative prior beta distribution with parameters `bg_protect_alpha` and `bg_protect_beta` in `background_model_params`.
#' @param background_model_params list containing parameters for background model which must contain following elements:
#' `bg_protect_prob_fixed` constant value for model parameter "bg_protect_prob" used in "fixed" background model.
#' 
#' `bg_protect_min` minimum allowed value for background emission probability for 1.
#' `bg_protect_max` maximum allowed value for background emission probability for 1.
#' `bg_protect_alpha` alpha parameter for prior beta distribution for model parameter "bg_protect_prob".
#' `bg_protect_beta` beta parameter for prior beta distribution for model parameter "bg_protect_prob".
#'
#' 
#' 
#' @param ftp_model_params list containing parameters for footprint model which must contain following elements:
#' 
#' `ftp_protect_min` minimum allowed value for footprint emission probability for 1.
#' `ftp_protect_max` maximum allowed value for footprint emission probability for 1.
#' `ftp_protect_mean` expected mean for footprint emission probability for 1.
#' `ftp_protect_var` expected variance for footprint emission probability for 1.
#'
#'
#'
#' @param ... parameters for `rstan::vb` function. Please check `?rstan::vb`.
#'
#' @return S4 class \code{rstan::stanfit} representing the fitted results. Pleae check \code{?rstan::stanfit}. 
#' @export
#'
#' @examples
#' 
#' \dontrun{
#'  
#' ## simple data with two ftp of 5 and 10 bps. 
#' ## The table below is a count table of observed occurrences of
#' ## 00, 01, 10 and 11 at spacings from 1 until 15.
#'  
#' ftp_5_10_data <- data.frame("S" = 1:15,
#'                             "N00" = c(1626964,1508381,
#'                                       1420066,1336897,
#'                                       1258679,1185045,
#'                                       1136784,1090029,
#'                                       1045066,1001670,
#'                                       959927,983020,
#'                                       1000410,1012326,
#'                                       1019633),
#'                             "N01" = c(0,113856,
#'                                       197657,276539,
#'                                       350664,420399,
#'                                       464873,507988,
#'                                       549466,589487,
#'                                       627997,601629,
#'                                       580965,565776,
#'                                       555217),
#'                             "N10" = c(0, 113921,
#'                                       197854,276862,
#'                                       351147,421042,
#'                                       465716,509005,
#'                                       550645,590837,
#'                                       629533,603328,
#'                                       582814,567737,
#'                                       557220),
#'                             "N11" = c(873036,758842,
#'                                       674423,594702,
#'                                       519510,448514,
#'                                       402627,357978,
#'                                       314823,273006,
#'                                       232543,257023,
#'                                       275811,289161,
#'                                       297930))
#' 
#' 
#' ## variational inference for footprints in the data
#' vb_output <- infer_footprints_stan_vb(joint_freq_table = ftp_5_10_data,
#'                                       footprint_prior_diralphas = c(10,rep(1,14)))
#'                                       
#' ## check output in shinyapp for 
#' library(shinystan)
#' launch_shinystan(vb_output)
#' 
#' 
#' 
#' 
#' }
#' 
infer_footprints_stan_vb <- function(joint_freq_table,
                                     footprint_prior_diralphas = c(200,rep(1,199)),
                                     background_model = c("informative_prior","fixed"),
                                     background_model_params = list("bg_protect_prob_fixed" = 0.01,
                                                                    
                                                                    "bg_protect_min" = 0.0,
                                                                    "bg_protect_max" = 0.4,
                                                                    "bg_protect_alpha" = 675,
                                                                    "bg_protect_beta" = 57126),
                                     
                                     ftp_model_params = list("ftp_protect_min" = 0.6,
                                                             "ftp_protect_max" = 1,
                                                             "ftp_protect_alpha" = 1,
                                                             "ftp_protect_beta" = 0.01),
                                     ...
){
  
  
  background_model <- match.arg(background_model)
  
  
  
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
  
  ## check ftp_model_params
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
  
  
  if(background_model == "fixed"){
    stan_model_name <- "ftp_inference_background_fixed"
    
    if(background_model_params[["bg_protect_prob_fixed"]] < 0 | background_model_params[["bg_protect_prob_fixed"]] > 1){
      stop(paste0("bg_protect_prob_fixed in list background_model_params must be in [0,1]!\n",
                  "Current value background_model_params[[\"bg_protect_prob_fixed\"]] = ",background_model_params[["bg_protect_prob_fixed"]]))
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
  } else if(background_model == "informative_prior"){
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
    
    
  } else {
    stop("Unknown background_model!")
  }
  
  
  
  
  ## get initial values for sampling
  
  stan_initvals <- init_ftp_abundance_for_inference(dir_alpha = footprint_prior_diralphas,
                                                    nchains = 1,
                                                    background_model = background_model,
                                                    background_model_params = background_model_params,
                                                    ftp_model_params = ftp_model_params
  )
  
  
  
  
  ## Stan ADVI
  stanVB_out <- rstan::vb(object = stanmodels[[stan_model_name]],
                          data = stan_inputdata,
                          init = stan_initvals[[1]],
                          ...)
  
  return(stanVB_out)
}
