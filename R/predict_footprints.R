#' @title Function for calculating coverage probabilities for footprints in single-molecule-footprinting (SMF) data
#'
#' @param data \code{matrix} or \code{list} containing SMF data. If data is a \code{matrix} rows represent reads (or fragments),
#' columns represent positions in region of interest (ROI). If data is a \code{list} than each element must contain vector with SMF data.
#' NOTE that \code{data} must contain data for informative positions '0' represents accessible (methylated) position,
#' and '1' represents protected (unmethylated) position "C" and all other positions must be filled with "NA".
#'
#' @param footprint_models A list containing footprint models for proteins.
#' Each element must have 3 slots:
#' PROTECT_PROB - a numeric vector with footprint emission probabilities to find protected position within a footprint;
#' COVER_PRIOR - prior coverage probability (abundance) (\code{numeric}) reflecting what fraction of reads you expect to be covered by a footprint;
#' NAME - name (\code{character}) of a model, e.g. "Nucleosome".
#' @param bgprotectprob background emission probability to find a protected position within open (accessible) regions
#' @param bgcoverprior prior probability for percentage of all fragments to be in a free (accessible, or background) state
#' @param report_prediction_in_flanks \code{logical} whether to return calculated start probabilities in left flanking region.
#' In order to take into account partial footprints at left edge of fragments the algorithm extends each fragment by maximum footprint length on the left side.
#' \code{report_prediction_in_flanks} controls whether calculated start probabilities in the left flanking region will be reported in the \code{START_PROB} .
#' @param ncpu number of threads to use
#' @param verbose verbose mode for bug fixing

#'
#' @return A list which contains 3 data frames:
#'         \code{START_PROB} - data frame with calculated start probabilities for each SMF molecule (column \code{seq}),
#'                             each position in ROI (column \code{pos}) and each footprint model, e.g. Nucleosome, background etc.
#'                             These probabilities reflect how likely it is to find a start in each fragment and at each position of a certain footprint model.
#'
#'         \code{COVER_PROB} - data frame with calculated coverage probabilities for each SMF molecule (column \code{seq}),
#'                             each position in ROI (column \code{pos}) and each footprint model, e.g. Nucleosome, background etc.
#'                             These probabilities reflect how likely it is that a certain position in an amplicon and certain fragment is covered by a certain footprint model.
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

  ### validate footprint models
  ftpvalout <- validate_footprint_models(footprint_models,
  																							bgprotectprob,
  																							bgcoverprior,
  																							verbose,
  																							add=coll)
  footprint_models <- ftpvalout[["footprint_models"]]
  start_priors <- ftpvalout[["start_priors"]]
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
           !is.null(out.list[["COVER_PROB"]])))){
    if(verbose)
      .message_timestamp("convert cpp_nomeR output to data.frame...")
    return(lapply(out.list,as.data.frame,
                  stringsAsFactors = FALSE,
                  check.names = FALSE))
  } else {
    stop("retrieved NULL results from C++ function.")
  }

}
