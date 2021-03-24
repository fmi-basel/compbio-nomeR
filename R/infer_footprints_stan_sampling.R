
#' Bayesian inference of footprint abundance using STAN algorithms
#'
#' This function uses \code{rstan} functionality for bayesian inference, in particular Hamiltonian Monte Carlo (HMC) and No-U_Turn Sampler (NUTS) algorithms, to infer footprint abundance in NOMe-Seq data using sampling from posterior distribution.
#'
#' @param joint_freq_table \code{data.frame} which must contain columns \code{S}, \code{N00},\code{N01},\code{N10},\code{N11}, 
#' where \code{S} represent spacings and \code{N00},\code{N01},\code{N10},\code{N11} are observed counts for 00, 01, 10 and 11 at spacing \code{S}.
#' This table can be obtained using function nomeR::count_joint_frequencies(...)
#' @param footprint_prior_diralphas vector containing parameters for dirichlet distribution which is used to model footprint abundances.
#' Default vector assumes that around 50\% of data is free of any footprints and all other footprint up-to 200bp have uninformative prior.
#' 
#' 
#' @param ftp_background_model type of model used for inference.
#' \describe{
#' \item{"informative_prior"}{both parameters \code{bg_protect_prob} and \code{footprint_protect_prob} as well as footprint abundances are subject for inference;}
#' \item{"bg_fixed"}{\code{bg_protect_prob} is fixed value defined by \code{background_model_params[["bg_protect_prob_fixed"]]} and \code{footprint_protect_prob} as well as footprint abundances are subject for inference;}
#' \item{"ftp_bg_fixed"}{\code{bg_protect_prob} and \code{footprint_protect_prob} are fixed values defined by corresponding values in \code{background_model_params} and \code{ftp_model_params} and inference runs only for footprint abundances.}
#' }
#' 
#' @param background_model_params list containing parameters for background model which must contain following elements:
#' \describe{
#' \item{`bg_protect_prob_fixed`}{constant value for model parameter "bg_protect_prob" used in "bg_fixed" and "ftp_bg_fixed" models.}
#' \item{`bg_protect_min`}{minimum allowed value for background emission probability for 1.}
#' \item{`bg_protect_max`}{maximum allowed value for background emission probability for 1.}
#' \item{`bg_protect_alpha`}{alpha parameter for prior beta distribution for model parameter "bg_protect_prob".}
#' \item{`bg_protect_beta`}{beta parameter for prior beta distribution for model parameter "bg_protect_prob".}
#' }
#' 
#' @param ftp_model_params list containing parameters for footprint model which must contain following elements:
#' \describe{
#' \item{`ftp_protect_prob_fixed`}{constant value for model parameter "footprint_protect_prob" used in "bg_fixed" and "ftp_bg_fixed" models.}
#' \item{`ftp_protect_min`}{minimum allowed value for footprint emission probability for 1.}
#' \item{`ftp_protect_max`}{maximum allowed value for footprint emission probability for 1.}
#' \item{`ftp_protect_alpha`}{alpha parameter for prior beta distribution for model parameter "footprint_protect_prob".}
#' \item{`ftp_protect_beta`}{beta parameter for prior beta distribution for model parameter "footprint_protect_prob".}
#' }
#'
#'
#'
#' @param nchains number of markov chains to run.
#' @param ncpu number of CPUs to use.
#' @param rstan_control a list of tuning parameters for stan sampler algorithm.

#' @param ... parameters for \code{\link[rstan]{sampling}} function. Please check \code{\link[rstan]{sampling}}.

#'
#'
#' @return S4 class \code{\link[rstan]{stanfit-class}} representing the inference results. Please check \code{\link[rstan]{sampling}}. 
#' @export
#'
#' @examples
#' 
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
#' hmc_output <- infer_footprints_stan_sampling(joint_freq_table = ftp_5_10_data,
#'                                       footprint_prior_diralphas = c(10,rep(1,14)))
#'                                       
#' ## check output in shinyapp for 
#' library(shinystan)
#' launch_shinystan(hmc_output)
#' 
#' 
#' 
#' 
#' }
#' 
#' 
infer_footprints_stan_sampling <- function(joint_freq_table,
                                           footprint_prior_diralphas = c(200,rep(1,199)),
                                           
                                           ftp_background_model = c("informative_prior","bg_fixed","ftp_bg_fixed"),
                                           background_model_params = list("bg_protect_prob_fixed" = 0.01,
                                                                          
                                                                          "bg_protect_min" = 1e-3,
                                                                          "bg_protect_max" = 0.4,
                                                                          "bg_protect_alpha" = 675,
                                                                          "bg_protect_beta" = 57126),
                                           
                                           ftp_model_params = list("ftp_protect_min" = 0.6,
                                                                   "ftp_protect_max" = 0.999,
                                                                   "ftp_protect_alpha" = 1,
                                                                   "ftp_protect_beta" = 0.01),
                                           nchains = 4,
                                           ncpu = nchains,
                                           rstan_control = list(max_treedepth = 15,
                                                                adapt_delta = 0.9),
                                           ...
){
  ftp_background_model <- match.arg(ftp_background_model)
  
  stan_input <- validate_construct_stan_input(joint_freq_table,
                                              footprint_prior_diralphas,
                                              ftp_background_model,
                                              background_model_params,
                                              ftp_model_params)
  
  ## get initial values for sampling
  stan_initvals <- init_ftp_abundance_for_inference(dir_alpha = footprint_prior_diralphas,
                                                    nchains = nchains,
                                                    ftp_background_model = ftp_background_model,
                                                    background_model_params = background_model_params,
                                                    ftp_model_params = ftp_model_params
  )
  
 
  ## hmc sampling
  stanfit_out <- rstan::sampling(object = stanmodels[[stan_input$stan_model_name]],
                                 data = stan_input$stan_inputdata,
                                 chains = nchains,
                                 init = stan_initvals,
                                 cores = ncpu,
                                 control = rstan_control,
                                 ...)
  
  return(stanfit_out)
}








