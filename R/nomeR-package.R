## R functions for running nomeR

#' nomeR: A package that uses a probabilistic model for interpreting NOME-Seq data
#'
#' The package uses a probabilistic model for calculating posterior probabilities that a certain position 
#' in a NOME-Seq read is covered by a certain footprint.
#' 
#' @section nomeR functions:
# \code{create_data_list} - function for preprocessing a NOME-Seq matrix containing information about
# "protected" and "unprotected" positions in each read. Rows must represent reads (fragments) and columns must represent positions in the amplicon.
#' 
#'   
#' \code{run_nomeR} - function which takes as input 1) a matrix with NOME-Seq data where rows represent reads (or fragments)
#' and columns represent positions in the amplicon, 2) models for a set of user-defined footprints,
#' and using a probabilistic model it returns a data.frame with posterior probabilities for each position and fragment to be covered by each footprint.
#' @docType package
#' @name nomeR
NULL


## usethis namespace: start
#' @useDynLib nomeR, .registration = TRUE
## usethis namespace: end
NULL


## usethis namespace: start
#' @importFrom Rcpp sourceCpp
## usethis namespace: end
NULL



#' @title Function for running probabilistic model to predict binding sites of DNA binding proteins in NOMe-Seq data
#'
#' @param data A matrix containing NOMe-Seq data. Rows represent reads (or fragments), 
#' columns represent positions in amplicon. NOTE that \code{data} must contain only positions of selected GpCs where '0' represent unprotected "C", 
#' '1' represent protected "C" and all other positions must be filled with "NA".
#' @param footprint_models A list containing footprint models for proteins. 
#' Each element must have 3 slots:
#' PROTECT_PROB - a numeric vector with probability to find protected 'C' within the footprint
#' PRIOR - prior (numeric) probability reflecting how often we expect to find a footprint
#' NAME - name (character) of a model, e.g. "Nucleosome"
#' @param bgprotectprob probability to find a protected 'C' in open (uprotected) regions
#' @param bgcoverprior prior probability for percentage of all fragments to be in a free (unprotected, or background) state
#' @param bound_fit_tol fitting tolerance for estimating initial values of partition sums
#' @param bound_min_fderiv_val This value used in Haley's numerical method for solving equations and represent a 
#' minimal value of denominator in first derivative of a function
#' @param bound_max_steps Maximum number of iterations for algorithm which estimates boundary values for partition sums
#' @param run_priorEM If TRUE the function runs Expectation Maximization algorithm for fitting prior probabilities of binding objects.
#' Not recommended for regions with transcription factor binding as the EM overestimates shorter footprints
#' @param priorEM_fit_tol Fitting tolerance for running prior EM 
#' @param priorEM_max_steps Maximum number of iterations in prior EM
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
#' nomeR.out <- nomeR_predict(data=rmatr,
#'                            footprint_models = ftp.models,
#'                            bgprotectprob = 0.05,
#'                            bgcoverprior = bg.pr)
#' 
#' 
nomeR_predict <- function(data,
                          footprint_models,
                  bgprotectprob,
                  bgcoverprior,
                  
                  bound_fit_tol = 1e-5,
                  bound_min_fderiv_val = 1e-5,
                  bound_max_steps = 1000,
                  run_priorEM = FALSE,
                  priorEM_fit_tol = 1e-5,
                  priorEM_max_steps = 1000
                  
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
  
  if(bound_fit_tol <=0){
    ArgumentCheck::addError(
      msg = "bound_fit_tol must be a positive number",
      argcheck = arg.check
    )
  }
  if(priorEM_fit_tol <=0){
    ArgumentCheck::addError(
      msg = "priorEM_fit_tol must be a positive number",
      argcheck = arg.check
    )
  }
  ArgumentCheck::finishArgCheck(arg.check)
  
  ## convert data into list of lists
  data <- .create_data_list(matr = data)
  
  ## convert cover priors to start priors and add them into bindin_models
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
  
  
  out.list <- run_cpp_nomeR(data, footprint_models, bgprotectprob, start_priors["BG"], bound_fit_tol, bound_min_fderiv_val, bound_max_steps, run_priorEM, priorEM_fit_tol, priorEM_max_steps)
  
  lapply(out.list,as.data.frame,stringsAsFactors = FALSE)
  
}



.start_prior2cover_prior <- function(start_prior,
                                     footprint_len){
  stopifnot(length(start_prior) == length(footprint_len))
  
  cover_prior <- start_prior * footprint_len
  cover_prior/sum(cover_prior)
}


.cover_prior2start_prior <- function(cover_prior,
                                     footprint_len){
  stopifnot(length(cover_prior) == length(footprint_len))
  stopifnot(all(footprint_len > 0))
  
  
  start_prior <- cover_prior/footprint_len
  start_prior/sum(start_prior)
}


#' Doing data preprocessing steps for NOME-Seq data to amke it suitable for the run_nomeR function
#'
#' @param matr A matrix containing NOMe-Seq data. Rows represent reads (or fragments), 
#' columns represent position in an amplicon. NOTE that \code{data} must contain only positions of selected GpCs where '0' represent unprotected "C", 
#' '1' represent protected "C" and all other positions must be filled with "NA".
#' @return TBA
#' 
#'
#' 
#' 
.create_data_list <- function(matr){
  ## check arguments
  arg.check <- ArgumentCheck::newArgCheck()
  if(!inherits(matr,"matrix")){
    ArgumentCheck::addError(
      msg = "matr must be a matrix",
      argcheck = arg.check
    )
  }
  if(ncol(matr) == 0 | nrow(matr) == 0){
    ArgumentCheck::addError(
      msg = "matr must contain some data",
      argcheck = arg.check
    )
  }
  if(all(is.na(matr))){
    ArgumentCheck::addError(
      msg = "matr contains all NAs",
      argcheck = arg.check
    )
  }
  if(any(!(matr %in% c(0,1,NA)))){
    ArgumentCheck::addError(
      msg = "matr must contain only 0, 1 or NA",
      argcheck = arg.check
    )
  }
  ArgumentCheck::finishArgCheck(arg.check)
  
  
  if(is.null(row.names(matr))){
    warning("matr does not have row names. Fill the row names with indexes.")
    row.names(matr) <- 1:nrow(matr)
  }
  
  
  dat.out <- lapply(row.names(matr),function(nm){
    seq <- matr[nm,]
    seq[is.na(seq)] <- 2
    seq <- paste(seq,collapse = "")
    list("DATA" = seq,
         "NAME" = nm)
  })
  dat.out
}



#' Function for creating contingency tables of frequencing have 0,0; 0,1, 0,NA etc with spacing S between them
#'
#' @param data A list containing preprocessed NOMe-Seq data returned from \code{create_data_list}.
#' @param maxspacing Maximum spcaing between positions.
#' @param maxwmlen A priori maximum length of possible footprint in the data.
#' If \code{maxwmlen} is not 0, then the sequences are flanked by NA filled artificial sequences on the left and on the right.
#'
#' @return data.frame containing frequencies for 0,0; 0,1; etc for each spacing between 0 (no gap, adjacent positions) and \code{maxspacing}
#'
#' 
#' 
count_spacing_frequencies <- function(data, maxspacing, maxwmlen = 0L){
  
  message(paste0("R function for counting spacing frequencies"))
  
  out_data <- count_spacing_freq_cpp(data, maxspacing, maxwmlen)
  as.data.frame(out_data)
}
