


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


.message_timestamp <- function(msg){
  message(paste0("[",Sys.time(),"]: ",msg))
}

.warning_timestamp <- function(msg){
  warning(paste0("[",Sys.time(),"]: ",msg))
}


.onUnload <- function (libpath) {
  library.dynam.unload("nomeR", libpath)
}


