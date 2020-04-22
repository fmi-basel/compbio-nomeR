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
#' 
#' 
#' 
#' 
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
#' 
#' @param footprint_models A list containing footprint models for proteins. 
#' Each element must have 3 slots:
#' PROTECT_PROB - a numeric vector with probability to find protected 'C' within the footprint;
#' COVER_PRIOR - prior coverage probability (numeric) reflecting what fraction of reads you expect to be covered by a footprint;
#' NAME - name (character) of a model, e.g. "Nucleosome".
#' @param bgprotectprob probability to find a protected 'C' in open (uprotected) regions
#' @param bgcoverprior prior probability for percentage of all fragments to be in a free (unprotected, or background) state
#' @param ncpu number of threads to use
#' @param bound_fit_tol fitting tolerance for estimating initial values of partition sums
#' @param bound_min_fderiv_val This value used in Haley's numerical method for solving equations and represent a 
#' minimal value of denominator in first derivative of a function
#' @param bound_max_steps Maximum number of iterations for algorithm which estimates boundary values for partition sums
#' @param run_priorEM If TRUE the function runs Expectation Maximization algorithm for fitting prior probabilities of binding objects.
#' Not recommended for regions with transcription factor binding as the EM overestimates shorter footprints.
#' @param priorEM_fit_tol Fitting tolerance for running prior EM 
#' @param priorEM_max_steps Maximum number of iterations in prior EM
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
                          ncpu = 1L,
                          bound_fit_tol = 1e-5,
                          bound_min_fderiv_val = 1e-5,
                          bound_max_steps = 1000,
                          run_priorEM = FALSE,
                          priorEM_fit_tol = 1e-5,
                          priorEM_max_steps = 1000,
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
                            bound_fit_tol, 
                            bound_min_fderiv_val, 
                            bound_max_steps, 
                            run_priorEM, 
                            priorEM_fit_tol, 
                            priorEM_max_steps, 
                            ncpu,
                            verbose)
  if(verbose)
    .message_timestamp("convert cpp_nomeR output to data.frame...")
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
#' @param ncpu number of cores
#' @return TBA
#' 
#'
#' 
#' 
.create_data_list <- function(matr,
                              ncpu = 1L){
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
  
  ## convert NAs to 2
  matr[which(is.na(matr),arr.ind = T)] <- 2
  
  
  
  dat.out <- parallel::mclapply(row.names(matr),
                                function(nm){
                                  seq <- matr[nm,]
                                  #seq[is.na(seq)] <- 2
                                  seq <- paste(seq,collapse = "")
                                  list("DATA" = seq,
                                       "NAME" = nm)
                                }, mc.cores = ncpu)
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
#' @export
#' 
#' 
count_joint_frequencies <- function(data, maxspacing, maxwmlen = 0L,ncpu=1L){
  
  #message(paste0("R function for counting spacing frequencies"))
  
  data <- .create_data_list(matr = data,
                            ncpu = ncpu)
  #browser()
  out_data <- count_spacing_freq_cpp(data, maxspacing, maxwmlen)
  return(as.data.frame(out_data))
}




#' Function for calculating theoretical joint probablities at certain distance
#' 
#' Using assumption that all non-background footprints have the same emission probabilities we calculate 
#' theoretical probabilities 
#' \code{P(0,0|S)}, \code{P(0,1|S)}, \code{P(1,0|S)}, \code{P(1,1|S)} for all distances \code{S} up to \code{max_spacing}.
#' \code{P(0,0|S)} etc, are probabilities to find \code{0} and \code{0} separated by distance \code{S}.
#' @param ftp_cover_priors coverage probabilities for each footprint. 
#' IMPORTANT! This must be a vector where \code{ftp_cover_priors[1]} is coverage probability of background and we assume that indexes in the vector correspond to the footprint lengths.
#' @param bg_protect_prob emission probability for background for \code{1}.
#' @param footprint_protect_prob emission probability for footprints for \code{1}.
#' @param max_spacing maximum distance
#'
#' @return
#' @export
#'
calculate_theor_joint <- function(ftp_cover_priors, # here vector of priors also represent lengths, namely ith element of the vector has length i+1, e.g. ftp_cover_priors[0] is a prior for bg with length 1. Make sure that R function passes correct vector with priors
                                  bg_protect_prob,
                                  footprint_protect_prob,
                                  max_spacing){
  
  # if(any(ftp_cover_priors < 0) | abs(sum(ftp_cover_priors) - 1) > .Machine$double.eps){
  #   stop("ftp_cover_priors must be vector of non-negative numbers which sum up to 1")
  # }
  
  theor_joint_prob <- calculate_theor_joint_prob_cpp(ftp_cover_priors,
                                                     bg_protect_prob,
                                                     footprint_protect_prob,
                                                     max_spacing)
  
  return(as.data.frame(theor_joint_prob))
  
}



#' Generate synthetic NOME-Seq like dataset
#'
#' @param amplicon_len length of a synthetic amplicon
#' @param n_reads number of reads to generate
#' @param footprint_models footprint models
#' @param bgprotectprob background emission probabilitity for 1
#' @param GpC_param if scalar between 0 and 1 treated as percentage of GpC pos. If vector of integers treated as predefined positions of GpC.
#' @param extend_ampl whether to extend amplicon by \code{max(length(footprint_models))} and fill it with \code{NA}.
#'
#' @return
#' @export
#'
#' @examples
#' 
#' \dontrun{
#' 
#' nucl_cover_prior <- 0.3
#' nucl_len <- 10
#' nucl_prot_prob <- 0.99
#' 
#' tf_cover_prior <- 0.05
#' tf_len <- 5
#' tf_prot_prob <- 0.99
#' 
#' 
#' bg_prot_prob <- 0.01
#' 
#' footprint_models <- list(list("PROTECT_PROB" = rep(nucl_prot_prob,nucl_len),
#'                               "COVER_PRIOR" = nucl_cover_prior,
#'                               "NAME" = "Nucleosome"),
#'                          list("PROTECT_PROB" = rep(tf_prot_prob,tf_len),
#'                               "COVER_PRIOR" = tf_cover_prior,
#'                               "NAME" = "TF"))
#' 
#' #### simulate nome data
#' amp.len <- 500
#' n.reads <- 5000
#' 
#' 
#' set.seed(338485)
#' synthetic_nome_data <- generate_synthetic_NOME_data(amplicon_len = amp.len,
#'                                                   n_reads = n.reads,
#'                                                   footprint_models = footprint_models,
#'                                                   bgprotectprob = bg_prot_prob)
#' 
#' }
#' 
generate_synthetic_NOME_data <- function(amplicon_len, # length of the apmlicon
                                       n_reads, # number of reads
                                       footprint_models,
                                       bgprotectprob,
                                       GpC_param = 1, # if scalar between 0 and 1 treated as percentage of GpC pos
                                       # if vector of integers treated as predefined positions of GpC
                                       extend_ampl = FALSE){
  #### CHECK ARGUMENTS ####
  arg.check <- ArgumentCheck::newArgCheck()
  if(!(inherits(amplicon_len,"numeric") | inherits(amplicon_len,"integer")) | amplicon_len <= 0){
    ArgumentCheck::addError(
      msg = "amplicon_len must be a positive number!",
      argcheck = arg.check
    )
    
  }
  
  if(!(inherits(n_reads,"numeric") | inherits(n_reads,"integer")) | n_reads <= 0){
    ArgumentCheck::addError(
      msg = "n_reads must be a positive number!",
      argcheck = arg.check
    )
    
  }
  
  ##### CHECK WHETHER GpC_param is correct #####
  if(length(GpC_param) == 0){
    ArgumentCheck::addError(
      msg = "GpC_param must not be empty!",
      argcheck = arg.check
    )
  } else if(length(GpC_param) == 1){
    if(GpC_param <= 0 | GpC_param > 1){
      ArgumentCheck::addError(
        msg = "GpC_param must be a scalar between 0 and 1 or a vector of integers or numerics.",
        argcheck = arg.check
      )
    }
  } else if(!(is.vector(GpC_param,"integer")  | is.vector(GpC_param,"numeric"))){
    ArgumentCheck::addError(
      msg = "GpC_param must be a scalar between 0 and 1 or a vector of integers or numerics.",
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
  
  ArgumentCheck::finishArgCheck(arg.check)
  
  #### GET GPC POSITIONS ####
  
  if(length(GpC_param) == 1 & GpC_param > 0 & GpC_param <= 1){
    GpC_pos <- sort(sample(1:amplicon_len,size=floor(GpC_param * amplicon_len),
                           replace=F)) # replace = F makes sure that positions will not repeat
    
  } else if(is.vector(GpC_param,"integer") | is.vector(GpC_param,"numeric")){
    GpC_param <- as.integer(GpC_param)
    if(any(!GpC_param %in% 1:amplicon_len)){
      warnings("GpC_param contains positions outsite amplicon length. Removing such positions")
    }
    GpC_pos <- GpC_param[GpC_param %in% 1:amplicon_len]
  } else {
    stop(paste0("GpC_param must be a scalar between 0 and 1 or a vector of integers or numerics."))
  }
  
  
  #### PROCESSING FOOTPRINT MODELS #####
  ## check whether any of the footprint has length 1 or name "BG"
  if(any(sapply(footprint_models,function(x){
    x$NAME
  }) == "BG") | any(sapply(footprint_models,function(x){
    length(x$PROTECT_PROB)
  }) == 1)){
    stop("Footprint with name BG or length 1 are reserved for background. Remove them from footprint models.")
  }
  
  
  ## adding BG
  bgcoverprior <- 1 - sum(sapply(footprint_models,function(x){
    x$COVER_PRIOR
  },simplify = T,USE.NAMES = F))
  
  footprint_models <- c(list(list("NAME" = "BG",
                                  "COVER_PRIOR" = bgcoverprior,
                                  "PROTECT_PROB" = bgprotectprob)),
                        footprint_models)
  #browser()
  model_names <- sapply(footprint_models,function(x){
    x$NAME
  })
  
  if(any(duplicated(model_names))){
    stop("Footprint names must not be duplicated!")
  }
  
  names(footprint_models) <- model_names
  
  model_cover_priors <- sapply(footprint_models,function(x){
    x$COVER_PRIOR
  },simplify = T,USE.NAMES = F)
  
  model_lengths <- sapply(footprint_models,function(x){
    length(x$PROTECT_PROB)
  },simplify = T,USE.NAMES = F)
  
  
  model_start_priors <- .cover_prior2start_prior(model_cover_priors,
                                                 model_lengths)
  
  ## decide whether start of footprints will be outside amplicon
  if(extend_ampl){
    max.model.len <- max(model_lengths)
  } else{
    max.model.len <- 1
  }
  ampl.len.ext <- amplicon_len + max.model.len - 1
  
  
  nome_data <- lapply(1:n_reads,function(i){
    
    dat <- matrix(NA,nrow = 1,ncol = ampl.len.ext)
    true_conf_pos <- vector(mode="integer")
    true_conf_model <- vector(mode = "character")
    pos <- 1
    while(pos <= ampl.len.ext){
      
      #browser()
      
      ## choose which model will fit
      if(extend_ampl){
        curr_model_start_prob <- model_start_priors
        curr_model_names <- model_names
      } else {
        which.fit <- which(pos + model_lengths - 1 <= ampl.len.ext)
        curr_model_names <- model_names[which.fit]
        curr_model_start_prob <- model_start_priors[which.fit]
      }
      
      ## choose which model to emit
      emit_model <- sample(x = curr_model_names,
                           size = 1,
                           prob = curr_model_start_prob)
      
      true_conf_pos <- c(true_conf_pos,pos - max.model.len + 1)
      true_conf_model <- c(true_conf_model,emit_model)
      
      model_len <- model_lengths[emit_model]
      model_prot_prob <- footprint_models[[emit_model]][["PROTECT_PROB"]]
      ## draw random data using binomial distribution
      
      model_synth_data <- rbinom(n = model_len, size = 1, prob = model_prot_prob)
      fill.idx <- pos:(min(pos + model_len - 1,ampl.len.ext))
      dat[1,fill.idx] <- model_synth_data[1:length(fill.idx)]
      pos <- pos + model_len
      names(pos) <- NULL
      
    }
    
    dat.out <- list("DATA" = dat,
                    "TRUE_CONF" = data.frame(seq = i,
                                             pos=true_conf_pos,
                                             model = true_conf_model,
                                             stringsAsFactors = F))
    
  })
  
  nome_sim_data <- do.call(rbind,lapply(nome_data,function(x){
    x[["DATA"]]
  }))
  true_conf <- do.call(rbind,lapply(nome_data,function(x){
    x[["TRUE_CONF"]]
  }))
  ## remove additional flanks
  nome_sim_data <- nome_sim_data[,max.model.len:ampl.len.ext,drop=F]
  ## set non GpC positions to NA
  nome_sim_data[,!(1:ncol(nome_sim_data) %in% GpC_pos)] <- NA
  
  ## return without additional flanking columns
  list("DATA" = nome_sim_data,
       "TRUE_CONF" = true_conf)
}




#' Test function to get prior estimates from empirical joint probs. Not bayesian setting
#'
#' @param count.data 
#' @param max.footprint.len 
#' @param bg.protect.prob 
#' @param footprint.protect.prob 
#'
#' @return
#' @export
#'
infer_cover_priors8 <- function(count.data,
                                max.footprint.len,
                                bg.protect.prob,
                                footprint.protect.prob
){
  #browser()
  P0 <- subset(count.data,S == 1)$P00
  P1 <- subset(count.data,S == 1)$P11
  
  count.data <- count.data[order(count.data$S),]
  
  beta1 <- bg.protect.prob
  beta0 <- 1-beta1
  alpha1 <- footprint.protect.prob
  alpha0 <- 1- alpha1
  
  alpha1_vec <- c(beta1,rep(alpha1,max.footprint.len - 1))
  alpha0_vec <- 1 - alpha1_vec
  
  bg_prior <- (P0 - alpha0)/(beta0 - alpha0)
  
  
  P11 <- count.data$P11 
  names(P11) <- count.data$S
  P00 <- count.data$P00
  names(P00) <- count.data$S
  
  
  sigma_vec <- setNames(vector(mode="numeric",
                               length = max(count.data$S)),
                        1:max(count.data$S))
  sigma_vec[] <- NA
  est_priors <- setNames(vector(mode="numeric",
                                length = max.footprint.len),
                         1:max.footprint.len)
  est_priors[] <- NA
  
  ### calculate initial conditions 
  bg_start_prior <- 1/(beta0 - alpha0) * (alpha1 * P00["2"] + alpha0 * (P11["2"] - alpha1))/(P00["2"] - P11["2"] + alpha1 - P00["1"])
  
  #sigma_vec["1"] <- 1/(beta0 - alpha0) * (alpha1 * P00["2"] + alpha0 * (P11["2"] - alpha1))/(P00["1"] - alpha0)
  sigma_vec["1"] <- bg_start_prior
  C1 <- unname(bg_prior * (alpha1 + (beta1 - alpha1)* bg_start_prior)/bg_start_prior)
  C0 <- unname(bg_prior * (alpha0 + (beta0 - alpha0)* bg_start_prior)/bg_start_prior)
  
  tmp.cnt.dat <- subset(count.data,S >=2)
  Y00 <- setNames((tmp.cnt.dat$P00 - P0 * alpha0)/(beta0 - alpha0),
                  tmp.cnt.dat$S)
  Y01 <- setNames((tmp.cnt.dat$P01 - P0 * alpha1)/(beta1 - alpha1),
                  tmp.cnt.dat$S)
  
  
  Y10 <- setNames((tmp.cnt.dat$P10 - P1 * alpha0)/(beta0 - alpha0),
                  tmp.cnt.dat$S)
  Y11 <- setNames((tmp.cnt.dat$P11 - P1 * alpha1)/(beta1 - alpha1),
                  tmp.cnt.dat$S)
  
  
  sigma_vec["2"] <- 1/C1 * (Y11["3"] + bg_prior * beta1 * sigma_vec["1"])
  est_priors["1"] <- bg_prior
  #browser()
  Y00 <- (Y00 + Y01)/2
  Y11 <- (Y11 + Y10)/2
  for(S in 3:(max.footprint.len + 1)){
    ## calc Z
    Z1 <- sum(est_priors[1:(S-2)]/(1:(S-2)) * sigma_vec[S-1:(S-2)] * alpha1_vec[1:(S-2)])
    Z0 <- sum(est_priors[1:(S-2)]/(1:(S-2)) * sigma_vec[S-1:(S-2)] * alpha0_vec[1:(S-2)])
    
    est_priors[S-1] <- (S-1)/sigma_vec[1] * (C0/(alpha0 * C1 - alpha1)) * (Y11[as.character(S+1)] - C1/C0 * Y00[as.character(S+1)] + Z1 - C1/C0 * Z0)
    
    ## shall we set it to 0 if it is negative?
    # if(est_priors[S-1] < 0){
    #   est_priors[S-1] <- 0
    # }
    
    sigma_vec[S] <- 1/C0 * (Y00[as.character(S+1)] +  Z0 + alpha0 * est_priors[S-1]/(S-1) * sigma_vec[1])
    
  }
  
  return(est_priors)
  
  #browser()
}





#' Add to a \code{data.frame} with empirical joint frequencies columns with empirical joint probabilities as well as mutual information 
#'
#' @param emp_joint_frequencies count table which consists empirical joint frequencies. It must contain the columns with names \code{S, N00, N01, N10, N11}.
#'
#' @return \code{data.frame} with additional columns
#' @export
#'
#' @examples
#' 
#' \dontrun{
#' count_table_spacings <- add_joint_probs_and_MI(emp_joint_frequencies = count_table_spacings)
#' }
#' 
#' 
add_joint_probs_and_MI <- function(emp_joint_frequencies){
  stopifnot(all(c("S","N00","N01","N10","N11") %in% colnames(emp_joint_frequencies)) | 
              all(c("S","P00","P01","P10","P11") %in% colnames(emp_joint_frequencies)))
  
  # tot.counts <- rowSums(wind.count.tbl[,c("N_0_0","N_0_1","N_1_0","N_1_1",
  #                                         "N_0_NA","N_1_NA","N_NA_0","N_NA_1",
  #                                         "N_NA_NA")])
  if(!all(c("P00","P01","P10","P11") %in% colnames(emp_joint_frequencies))){
    tot.counts <- rowSums(emp_joint_frequencies[,c("N00","N01","N10","N11")])
    
    # wind.count.tbl$P_0_0 <- with(data = wind.count.tbl,
    #                              (N_0_0 + 0.5*N_0_NA + 0.5*N_NA_0 + 0.25 * N_NA_NA)/tot.counts)
    # wind.count.tbl$P_0_1 <- with(data = wind.count.tbl,
    #                              (N_0_1 + 0.5*N_0_NA + 0.5*N_NA_1 + 0.25 * N_NA_NA)/tot.counts)
    # wind.count.tbl$P_1_0 <- with(data = wind.count.tbl,
    #                              (N_1_0 + 0.5*N_1_NA + 0.5*N_NA_0 + 0.25 * N_NA_NA)/tot.counts)
    # wind.count.tbl$P_1_1 <- with(data = wind.count.tbl,
    #                              (N_1_1 + 0.5*N_1_NA + 0.5*N_NA_1 + 0.25 * N_NA_NA)/tot.counts)
    emp_joint_frequencies$P00 <- with(data = emp_joint_frequencies,
                                      (N00)/tot.counts)
    emp_joint_frequencies$P01 <- with(data = emp_joint_frequencies,
                                      (N01)/tot.counts)
    emp_joint_frequencies$P10 <- with(data = emp_joint_frequencies,
                                      (N10)/tot.counts)
    emp_joint_frequencies$P11 <- with(data = emp_joint_frequencies,
                                      (N11)/tot.counts)
  }
  indep.prob <- with(emp_joint_frequencies,
                     data.frame("P0any" = P00 + P01,
                                "P1any" = P10 + P11,
                                "Pany0" = P00 + P10,
                                "Pany1" = P01 + P11,
                                stringsAsFactors = F))

  emp_joint_frequencies <- cbind(emp_joint_frequencies,
                                 indep.prob)
  emp_joint_frequencies$MI <- with(emp_joint_frequencies,
                            ifelse(P00 == 0,0,P00 * log(P00/(P0any * Pany0))) + 
                              ifelse(P01 == 0,0,P01 * log(P01/(P0any * Pany1))) +
                              ifelse(P10 == 0,0,P10 * log(P10/(P1any * Pany0))) +
                              ifelse(P11 == 0,0,P11 * log(P11/(P1any * Pany1))))
  
  emp_joint_frequencies$H_ij <- -with(emp_joint_frequencies,
                               ifelse(P00 == 0,0,P00 * log(P00)) + 
                                 ifelse(P01 == 0,0,P01 * log(P01)) +
                                 ifelse(P10 == 0,0,P10 * log(P10)) +
                                 ifelse(P11 == 0,0,P11 * log(P11)))
  
  emp_joint_frequencies$IQR <- with(emp_joint_frequencies,
                                    MI/H_ij)
  emp_joint_frequencies
}





.message_timestamp <- function(msg){
  message(paste0("[",Sys.time(),"]: ",msg))
}

.warning_timestamp <- function(msg){
  warning(paste0("[",Sys.time(),"]: ",msg))
}


.onUnload <- function (libpath) {
  library.dynam.unload("nomeR", libpath)
}

