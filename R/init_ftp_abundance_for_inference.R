#' Create initial values for HMC sampling
#'
#' @param dir_alpha vector containing parameters for dirichlet distribution which is used to model footprint abundances.
#' Default vector assumes that around 50\% of data is free of any footprints and all other footprint up-to 200bp have uninformative prior.
#' 
#' @param nchains number of markov chains which will be run.
#' @param ftp_protect_min minimum allowed value for footprint emission probability for 1.
#' @param ftp_protect_max maximum allowed value for footprint emission probability for 1.
#' @param ftp_protect_mean expected mean for footprint emission probability for 1.
#' @param ftp_protect_var expected variance for footprint emission probability for 1.
#'
#' @return list of length `nchains` containing initial values for `ftp_abundances` and `footprint_protect_prob`.
#' @export
#'
#' @examples
#' 
#' ## parameters for dirichlet distribution
#' footprint_prior_diralphas <- c(10,rep(1,14))
#' 
#' #' ## get initial values for sampling for vb
#' stan_initvals <- init_ftp_abundance_for_inference(dir_alpha = footprint_prior_diralphas)
#' 
#' 
init_ftp_abundance_for_inference <- function(dir_alpha = c(200,rep(1,199)),
                                          nchains = 4,
                                          ftp_protect_min = 0.6,
                                          ftp_protect_max = 0.9999,
                                          ftp_protect_mean = 0.95,
                                          ftp_protect_var = 0.01){
  
  ## check footprint_prior_diralphas
  if(length(dir_alpha) == 0){
    stop("dir_alpha seems to be empty")
  }
  if(any(dir_alpha <= 0)){
    stop("Incorrect dir_alpha! non-positive values are prohibited in footprint_prior_diralphas.")
  }
  
  
  init_vals <- lapply(1:nchains,function(x){
    ftp_abundances <- as.vector(extraDistr::rdirichlet(1,dir_alpha))
    
    ## param for ftp beta
    if(ftp_protect_var < ftp_protect_mean * (1 - ftp_protect_mean)){
      ftp_nu <- ftp_protect_mean * (1 - ftp_protect_mean)/ftp_protect_var - 1
      ftp_protect_alpha <- ftp_protect_mean * ftp_nu
      ftp_protect_beta <- (1 - ftp_protect_mean) * ftp_nu
    } else{
      stop("Incorrect values for ftp_protect_mean and ftp_protect_var! Inequality ftp_protect_var < ftp_protect_mean * (1 - ftp_protect_mean) must hold.")
    }
    
    footprint_protect_prob <- truncdist::rtrunc(1,spec="beta",a = ftp_protect_min, b=ftp_protect_max,
                                                shape1 = ftp_protect_alpha, shape2 = ftp_protect_beta)
    return(list("ftp_abundances" = ftp_abundances,
                "footprint_protect_prob" = footprint_protect_prob))
    
  })
  
  return(init_vals)
}
