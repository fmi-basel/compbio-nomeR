


#' Generate synthetic NOMe-Seq like dataset
#'
#' @param amplicon_len length of a synthetic amplicon
#' @param n_reads number of reads to generate
#' @param footprint_models footprint models
#' @param bgprotectprob background emission probabilitity for 1
#' @param GpC_param if scalar between 0 and 1 treated as percentage of GpC pos. If vector of integers treated as predefined positions of GpC.
#' @param extend_ampl whether to extend amplicon by \code{max(length(footprint_models))} and fill it with \code{NA}.
#'
#' @return TODO
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
#' synthetic_nome_data <- generate_synthetic_NOMe_data(amplicon_len = amp.len,
#'                                                   n_reads = n.reads,
#'                                                   footprint_models = footprint_models,
#'                                                   bgprotectprob = bg_prot_prob)
#' 
#' }
#' 
generate_synthetic_NOMe_data <- function(amplicon_len, # length of the apmlicon
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
      
      model_synth_data <- stats::rbinom(n = model_len, size = 1, prob = model_prot_prob)
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






