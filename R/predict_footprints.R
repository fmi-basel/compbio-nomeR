#' @title Function for running probabilistic model to predict binding sites of DNA binding proteins in NOMe-Seq data
#'
#' @param data A matrix containing NOMe-Seq data. Rows represent reads (or fragments), 
#' columns represent positions in amplicon. NOTE that \code{data} must contain only positions of selected GpCs where '0' represent unprotected "C", 
#' '1' represent protected "C" and all other positions must be filled with "NA".
#' 
#' @param footprint_models A list containing footprint models for proteins. 
#' Each element must have 3 slots:
#' PROTECT_PROB - a numeric vector with probability to find protected 'C' within the footprint;
#' COVER_PRIOR - prior coverage probability (numeric) reflecting what fraction of reads you expect to be covered by a footprint;
#' NAME - name (character) of a model, e.g. "Nucleosome".
#' @param bgprotectprob probability to find a protected 'C' in open (uprotected) regions
#' @param bgcoverprior prior probability for percentage of all fragments to be in a free (unprotected, or background) state
#' @param ncpu number of threads to use
#' @param verbose verbose mode for bug fixing
#'
#' @return A list which contains 3 data frames:
#'         \code{START_PROB} - data frame with calculated start probabilities for each NOME-Seq fragment (column \code{seq}),
#'                             each position in the amplicon (column \code{pos}) and each binding model, e.g. Nucleosome, background etc.
#'                             These probabilities reflect how likely it is to find a start in each fragment and at each position of a certain bidning model.
#'          
#'         \code{COVER_PROB} - data frame with calculated coverage probabilities for each NOME-Seq fragment (column \code{seq}),
#'                             each position in the amplicon (column \code{pos}) and each binding model, e.g. Nucleosome, background etc.
#'                             These probabilities reflect how likely it is that a certain position in an amplicon and certain fragment is covered by a certain bidning model.
#' 
#'         \code{SUMMARY} - data frame with summary information for each binding model and for the whole amplicon.
#'                          
#'                          
#'                          Column \code{Statistics} contains the following information:
#'                          
#'                          \code{Prior} - empirical prior probabilities for each model.
#'                          
#'                          
#'                          \code{Expected number of sites} - expected number of sites for each model, i.e. sum of all start probabilities.
#'                          
#'                          
#'                          \code{Coverage} - average coverage of the amiplicon.
#'                          
#'                          
#'                          \code{Sites<0.5} - expected number of sites with start probabilities lower than 0.5.
#'                          
#'                          
#'                          \code{Sites>=0.5} - expected number of sites with start probabilities higher than 0.5.
#'                          
#'                          
#'                          \code{Positions<0.5} - number of positions with start probabilities lower than 0.5.
#'                          
#'                          
#'                          \code{Positionds>=} - number of positions with start probabilities higher than 0.5.
#'                          
#'         
#'                          
#' 
#' @export
#'
#' 
#' 
#' 
#' @examples
#' set.seed(3346)
#' nc <- 50
#' nr <- 50
#' rmatr <- matrix(data = as.integer(rnorm(nc * nr) >= 0.5),
#'                 ncol = nc,nrow=nr)
#' 
#' ## create dummy footprints
#' bg.pr <- 0.5
#' ft.pr <- 1-bg.pr
#' ft.len <- 15
#' 
#' ## creating a list of binding models for nomeR
#' ftp.models <- list(list("PROTECT_PROB" = rep(0.99,ft.len),
#'                         "COVER_PRIOR" = ft.pr,
#'                         "NAME" = "FOOTPRINT"))
#' 
#' nomeR.out <- predict_footprints(data=rmatr,
#'                                 footprint_models = ftp.models,
#'                                 bgprotectprob = 0.05,
#'                                 bgcoverprior = bg.pr)
#' 
#' 
predict_footprints <- function(data,
                          footprint_models,
                          bgprotectprob,
                          bgcoverprior,
                          ncpu = 1L,
                          verbose = F
                          
) {
  
  ## check arguments
  arg.check <- ArgumentCheck::newArgCheck()
  
  ##### CHECKING WHETHER data HAS CORRECT FORMAT #####
  ## check whether data is a matrix with enough rows and columns
  if(!inherits(data,"matrix") | nrow(data) <= 1 | ncol(data) <= 1){
    # stop("Error! data must be a matrix with at least 2 columns and 2 rows!")
    ArgumentCheck::addError(
      msg = "data must be a matrix with at least one column and one row.",
      argcheck = arg.check
    )
  }
  
  
  ##### CHECK THAT binding_models is correct #######
  if(!inherits(footprint_models,"list")){
    # stop("Error! data must be a list!")
    ArgumentCheck::addError(
      msg = "footprint_models must be a list!",
      argcheck = arg.check
    )
  }
  if(length(footprint_models) == 0){
    # stop("Error! data has length 0")
    ArgumentCheck::addError(
      msg = "footprint_models has length 0",
      argcheck = arg.check
    )
  }
  if(any(sapply(footprint_models,function(x){
    !inherits(x,"list")
  },simplify = T))){
    ArgumentCheck::addError(
      msg = "footprint_models must be a list of lists",
      argcheck = arg.check
    )
    # stop("Error! data must be a list of lists")
  }
  ## check whether binding_models has all required slots
  check.res <- sapply(footprint_models,function(x){
    if(any(!(c("PROTECT_PROB","COVER_PRIOR","NAME") %in% names(x)))){
      ArgumentCheck::addError(
        msg = "Each item in footprint_models must contain entries PROTECT_PROB, COVER_PRIOR and NAME",
        argcheck = arg.check
      )
      
    }
    
    ## whether PROTECT_PROB has length > 0
    if(length(x[["PROTECT_PROB"]]) == 0){
      ArgumentCheck::addError(
        msg = "PROTECT_PROB of NAME=",x[["NAME"]]," has length 0",
        argcheck = arg.check
      )
    }
    ## whether NAME is string
    if(!inherits(x[["NAME"]],"character")){
      ArgumentCheck::addError(
        msg = paste0("NAME must be a string"),
        argcheck = arg.check
      )
    }
    if(!inherits(x[["PROTECT_PROB"]],"numeric")){
      ArgumentCheck::addError(
        msg = paste0(x[["NAME"]],": PROTECT_PROB must contain vector of numeric values"),
        argcheck = arg.check
      )
    }
    
    ## whether COVER_PRIOR is between 0 and 1
    if(!inherits(x[["COVER_PRIOR"]],"numeric") | x[["COVER_PRIOR"]] < 0 | x[["COVER_PRIOR"]] > 1){
      ArgumentCheck::addError(
        msg = paste0(x[["NAME"]],": COVER_PRIOR must be a value between 0 and 1."),
        argcheck = arg.check
      )
    }
    
  })
  
  #### check other args
  # bgprotectprob,
  # bgprior,
  # bound_fit_tol = 1e-5,
  # bound_min_fderiv_val = 1e-5,
  # bound_max_steps = 1000,
  # run_priorEM = FALSE,
  # priorEM_fit_tol = 1e-5,
  # priorEM_max_steps = 1000
  
  if(!inherits(bgprotectprob,"numeric") | bgprotectprob < 0 | bgprotectprob > 1){
    ArgumentCheck::addError(
      msg = "bgprotectprob be a numeric value between 0 and 1.",
      argcheck = arg.check
    )
  }
  
  if(!inherits(bgcoverprior,"numeric") | bgcoverprior < 0 | bgcoverprior > 1){
    ArgumentCheck::addError(
      msg = "bgcoverprior be a numeric value between 0 and 1.",
      argcheck = arg.check
    )
  }
  
  
  
  if(ncpu < 0   ){
    ArgumentCheck::addError(
      msg = "ncpu must be non-negative integer",
      argcheck = arg.check
    )
  }
  
  avail_ncpu <- parallel::detectCores()
  
  if(is.na(avail_ncpu)){
    .warning_timestamp("Could not detect number of available cpu. Setting ncpu to 1L.")
    ncpu <- 1L
  } else if(ncpu > avail_ncpu | ncpu == 0) {
    .warning_timestamp("Number of ncpu is 0 or exceeds number of available cpu. Setting ncpu to number of available cpus.")
    ncpu <- avail_ncpu
  }
  
  
  ArgumentCheck::finishArgCheck(arg.check)
  
  ## convert data into list of lists
  if(verbose)
    .message_timestamp("Converting data matrix to data list...")
  data <- .create_data_list(matr = data,
                            ncpu = ncpu)
  
  ## convert cover priors to start priors and add them into binding_models
  if(verbose)
    .message_timestamp("convert cover priors to start priors and add them into binding_models...")
  cover_priors <- c("BG" = bgcoverprior,
                    sapply(footprint_models,function(x){x[["COVER_PRIOR"]]}))
  footpr_lens <- c("BG" = 1,
                   sapply(footprint_models,function(x){length(x[["PROTECT_PROB"]])}))
  start_priors <- .cover_prior2start_prior(cover_prior = cover_priors,
                                           footprint_len = footpr_lens)
  
  footprint_models <- sapply(seq_along(footprint_models),
                             function(i){
                               ft.model <- footprint_models[[i]]
                               c(ft.model,
                                 list("PRIOR" = start_priors[i+1]))
                             },
                             simplify = F,USE.NAMES = T)
  
  if(verbose)
    .message_timestamp("Calling nomeR cpp code...")
  out.list <- run_cpp_nomeR(data,
                            footprint_models,
                            bgprotectprob,
                            start_priors["BG"],
                            # bound_fit_tol, 
                            # bound_min_fderiv_val, 
                            # bound_max_steps, 
                            # run_priorEM, 
                            # priorEM_fit_tol, 
                            # priorEM_max_steps, 
                            ncpu,
                            verbose)
  if(verbose)
    .message_timestamp("convert cpp_nomeR output to data.frame...")
  lapply(out.list,as.data.frame,
         stringsAsFactors = FALSE,
         check.names = FALSE)
  
}
