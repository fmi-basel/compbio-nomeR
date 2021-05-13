#' @title Function for running probabilistic model to predict binding sites of DNA binding proteins in NOMe-Seq data
#'
#' @param data \code{matrix} or \code{list} containing NOMe-Seq data. If data is a \code{matrix} rows represent reads (or fragments), 
#' columns represent positions in amplicon. If data is a \code{list} than each element must contain vector with NOMe-seq data. 
#' NOTE that \code{data} must contain only positions of selected GpCs where '0' represent unprotected "C", 
#' '1' represent protected "C" and all other positions must be filled with "NA".
#' 
#' @param footprint_models A list containing footprint models for proteins. 
#' Each element must have 3 slots:
#' PROTECT_PROB - a numeric vector with probability to find protected 'C' within the footprint;
#' COVER_PRIOR - prior coverage probability (numeric) reflecting what fraction of reads you expect to be covered by a footprint;
#' NAME - name (character) of a model, e.g. "Nucleosome".
#' @param bgprotectprob probability to find a protected 'C' in open (uprotected) regions
#' @param bgcoverprior prior probability for percentage of all fragments to be in a free (unprotected, or background) state
#' @param report_prediction_in_flanks \code{logical} whether to return calculated coverage and start probabilities in flanking regions.
#' In order to take into account partial footprints at edges of fragments the algorithms extends each fragment by maximum footprint lengths from both sides.
#' \code{report_prediction_in_flanks} controls whether calculated probabilities in flanking regions will be returned or not.
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
                          report_prediction_in_flanks = FALSE,
                          ncpu = 1L,
                          verbose = F
                          
) {
  
  ## check arguments
  
  coll = checkmate::makeAssertCollection()
  ### validate data
  data <- validate_data_for_predict(data)
  
  ### validate binding_models
  if(checkmate::check_list(footprint_models,types = "list",any.missing = FALSE, min.len = 1)){
    ftpnames <- sapply(footprint_models,function(x){
      
      checkmate::assert_names(names(x),must.include = c("PROTECT_PROB","COVER_PRIOR","NAME"),add = coll)
      checkmate::assert_numeric(x = x[["PROTECT_PROB"]], lower = 0, upper = 1, any.missing = FALSE, add = coll)
      checkmate::assert_number(x = x[["COVER_PRIOR"]], lower = 0, upper = 1, na.ok = FALSE, add = coll)
      return(x[["NAME"]])
    })
    checkmate::assert_character(x=ftpnames,any.missing = FALSE,unique = TRUE,add = coll)
  }
  
  ### validate bgprotectprob
  checkmate::assert_number(x = bgprotectprob,na.ok = FALSE, lower = 0, upper = 1, add = coll)
  ### validate bgcoverprior
  checkmate::assert_number(x = bgcoverprior,na.ok = FALSE, lower = 0, upper = 1, add = coll)
  ### validate report_prediction_in_flanks
  checkmate::assert_logical(report_prediction_in_flanks,any.missing = FALSE,all.missing = FALSE,len=1,add = coll)
  
  ### validate ncpu
  checkmate::assert_int(x=ncpu,lower = 0,na.ok = TRUE, add = coll)
  avail_ncpu <- parallel::detectCores()
  if(is.na(avail_ncpu)){
    .warning_timestamp("Could not detect number of available cpu. Setting ncpu to 1L.")
    ncpu <- 1L
  } else if(ncpu > avail_ncpu | ncpu == 0) {
    .warning_timestamp("Number of ncpu is 0 or exceeds number of available cpu. Setting ncpu to number of available cpus.")
    ncpu <- avail_ncpu
  }
  
  ## finish argument check
  checkmate::reportAssertions(coll)
  
  
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
    .message_timestamp("Calling run_cpp_nomeR...")
  out.list <- run_cpp_nomeR(data[["data_list"]],
                            data[["fragnames"]],
                            footprint_models,
                            bgprotectprob,
                            start_priors["BG"],
                            report_prediction_in_flanks,
                            ncpu,
                            verbose)
  
  if(all(c(!is.null(out.list[["START_PROB"]]),
           !is.null(out.list[["COVER_PROB"]]),
           !is.null(out.list[["SUMMARY"]])))){
    if(verbose)
      .message_timestamp("convert cpp_nomeR output to data.frame...")
    return(lapply(out.list,as.data.frame,
                  stringsAsFactors = FALSE,
                  check.names = FALSE))
  } else {
    stop("retrieved NULL results from C++ function.")
  }
  
}
