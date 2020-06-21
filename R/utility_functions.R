


## function for converting vector of start priors to cover priors
.start_prior2cover_prior <- function(start_prior,
                                     footprint_len){
  stopifnot(length(start_prior) == length(footprint_len))
  
  cover_prior <- start_prior * footprint_len
  cover_prior/sum(cover_prior)
}

## function for converting vector of covert priors to start priors
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





.message_timestamp <- function(msg){
  message(paste0("[",Sys.time(),"]: ",msg))
}

.warning_timestamp <- function(msg){
  warning(paste0("[",Sys.time(),"]: ",msg))
}


.onUnload <- function (libpath) {
  library.dynam.unload("nomeR", libpath)
}


