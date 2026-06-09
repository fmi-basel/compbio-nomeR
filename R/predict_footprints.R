#' Calculate posterior probabilities and predict footprints in single-molecule
#' footprinting (SMF) data
#'
#' @description
#' Computes posterior probabilities and predicts footprint configurations
#' in single-molecule footprinting (SMF) datasets.
#'
#' @param data A \code{matrix} or \code{list} containing binarized SMF data.
#'   - If a \code{matrix}, rows correspond to reads (or fragments) and columns
#'   correspond to positions in the region of interest (ROI).
#'   - If a \code{list}, each element must be a vector containing SMF data.
#'   Data must be binarized: '1' represents accessible (methylated) positions,
#'   '0' represents protected (unmethylated) positions, and all other positions
#'   should be \code{NA}.
#' @param footprint_models A list of footprint models for proteins. Each element
#'   must contain:
#'   \describe{
#'     \item{PROTECT_PROB}{Numeric vector of footprint emission probabilities
#'     for protected positions within a footprint.}
#'     \item{COVER_PRIOR}{Numeric prior probability reflecting the expected
#'     fraction of reads covered by the footprint.}
#'     \item{NAME}{Unique name of the model (character), e.g.,
#'     "Nucleosome--149", "Nucleosome--150".}
#'     \item{GROUP}{Optional non-unique group name (character) used for
#'     aggregating probabilities if \code{aggrByGroup = TRUE}. For example,
#'     footprints with names "Nucleosome--149" and "Nucleosome--150" may share
#'     the GROUP "Nucleosome".}
#'   }
#'   NOTE: To account for very short footprints that likely reflect correlated
#'   noise within accessible regions, users can define short footprints (usually 2–5 bp)
#'   and set \code{GROUP} to "background". In this case, if \code{aggrByGroup} is
#'   \code{TRUE}, the posteriors for background will additionally be aggregated across
#'   these short footprints.
#' @param bgprotectprob Background emission probability of a protected position
#'   in open (accessible) regions.
#' @param bgcoverprior Prior probability of a position being in the free
#'   (accessible/background) state.
#' @param aggrByGroup Logical. If \code{TRUE}, probabilities are aggregated by
#'   the \code{GROUP} ID in \code{footprint_models}. If \code{FALSE}, or if
#'   GROUP IDs are missing, probabilities are reported for each individual
#'   footprint \code{NAME}.
#' @param keepStartProb Logical. If \code{TRUE}, posteriors for footprint
#'   start positions are included in the return value.
#' @param ftpConfigMethod Algorithm for decoding footprint configurations:
#'   \describe{
#'     \item{\code{PV}}{Posterior-Viterbi (default): uses footprint coverage
#'     posteriors and a Viterbi-like algorithm to find a valid footprint
#'     configuration maximizing posterior coverage probability (see Fariselli
#'     et al, 2005). The \code{score} for each footprint is the geometric mean
#'     of posterior coverages across all positions covered by the footprint.}
#'     \item{\code{PosteriorDecoding}}{Posterior decoding using calculated
#'     coverage posteriors. The algorithm returns segments where coverage posteriors
#'     for a given footprint are maximum compared to other footprints.
#'     If coverage posteriors are aggregated by \code{GROUP}, i.e. \code{aggrByGroup}
#'     is \code{TRUE}, the \code{ftp_name} and \code{ftp_group} are identical and
#'     correspond to groups defined in \code{footprint_models}. NOTE: The widths of
#'     the segments do not necessarily match the footprint lengths defined by
#'     \code{footprint_models}. The \code{score} for each segment is the geometric
#'     mean of posterior coverage for the corresponding footprint across all positions
#'     within the segment.}
#'     \item{\code{Viterbi}}{Classic Viterbi algorithm to find the most probable
#'     footprint configuration. The \code{score} corresponds to the posterior
#'     start probability for each reported footprint.}
#'   }
#' @param ncpu Number of threads to use.
#' @param verbose Logical. If \code{TRUE}, enables verbose output for debugging.
#'
#' @return A list containing three data frames:
#'   \describe{
#'     \item{\code{COVER_PROB}}{Contains coverage probabilities for each SMF
#'     molecule (\code{seq}), each position in ROI (\code{pos}), and each
#'     footprint model. These probabilities indicate how likely a position is
#'     covered by a given footprint.}
#'     \item{\code{FOOTPRINT_CONF}}{Contains footprint configurations predicted
#'     for each molecule using the selected \code{ftpConfigMethod}. Each row
#'     reports the SMF molecule (\code{seq}), start position (\code{start}),
#'     width (\code{width}), footprint name (\code{ftp_name}), group
#'     (\code{ftp_group}), and confidence score (\code{score}) between 0 and 1.}
#'     \item{\code{START_PROB} (optional, if \code{keepStartProb = TRUE}.)}{
#'     Contains start probabilities for each SMF
#'     molecule (\code{seq}), each position in ROI (\code{pos}), and each
#'     footprint model. These probabilities indicate how likely a footprint
#'     starts at each position in a fragment.}
#'   }
#'
#' @importFrom checkmate makeAssertCollection assert_logical assert_int
#'     reportAssertions
#' @importFrom parallel detectCores
#' @importFrom cli cli_warn
#'
#' @references
#' Fariselli, P., Martelli, P. L., & Casadio, R. (2005).
#' *A new decoding algorithm for Hidden Markov Models improves the prediction of the topology of all-beta membrane proteins.*
#' BMC Bioinformatics, 6(S4), S12. https://doi.org/10.1186/1471-2105-6-S4-S12
#'
#' @export
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
predict_footprints <- function(data,
                               footprint_models,
                               bgprotectprob,
                               bgcoverprior,
                               aggrByGroup = TRUE,
                               ftpConfigMethod = c("PV","PosteriorDecoding", "Viterbi"),
                               keepStartProb = FALSE,
                               ncpu = 1L,
                               verbose = FALSE) {

    cli::cli_warn(c(
        "!" = "{.fn predict_footprints} is a legacy interface and will be deprecated in a future release.",
        "i" = "Input {.arg data} is expected to contain accessibility values: {.val 1} = accessible, {.val 0} = protected.",
        "i" = "Use {.fn predict_footprints_SE} for the current interface."
    ))

    ## check arguments
    coll <- makeAssertCollection()
    ### validate data. The output is a list
    ### ("nonNA_data" = nonNA_data,"fragnames" = fragnames)
    data <- validate_prepare_listOrMat(data)

    ftpConfigMethod <- match.arg(ftpConfigMethod)
    ### validate footprint models
    ftpvalout <- validate_footprint_models(footprint_models,
                                           bgprotectprob,
                                           bgcoverprior,
                                           aggrByGroup,
                                           verbose,
                                           add = coll)
    footprint_models <- ftpvalout[["footprint_models"]]
    start_priors <- ftpvalout[["start_priors"]]

    ### validate ncpu
    assert_int(x = ncpu, lower = 0, na.ok = TRUE, add = coll)
    avail_ncpu <- parallel::detectCores()
    if (is.na(avail_ncpu)) {
        .warning_timestamp(
            "Could not detect number of available cpu. Setting ncpu to 1L.")
        ncpu <- 1L
    } else if (ncpu > avail_ncpu || ncpu == 0) {
        .warning_timestamp(c("Number of ncpu is 0 or exceeds number of ",
                             "available cpu. Setting ncpu to number of ",
                             "available cpus."))
        ncpu <- avail_ncpu
    }

    ## finish argument check
    reportAssertions(coll)

    if (verbose) {
        .message_timestamp("Calling run_cpp_nomeR...")
    }
    ## the C++ needs only fidx_glob, fragpos, mod_prob
    predict_res_list <- calcStartCoverProbs_cpp(
        data[["nonNA_data"]][,"fidx_glob"], ## unique fragment ID or index
        data[["nonNA_data"]][,"fragpos"],   ## position within fragment, 1 - based
        data[["nonNA_data"]][,"mod_prob"],  ## modification probability in [0,1]
        footprint_models,
        bgprotectprob,
        start_priors["BG"],
        ftpConfigMethod,
        aggrByGroup,
        keepStartProb,
        ncpu,
        verbose)

    if (all(c(!is.null(predict_res_list[["START_PROB"]]),
              !is.null(predict_res_list[["COVER_PROB"]]),
              !is.null(predict_res_list[["FOOTPRINT_CONF"]])))) {
        if (verbose) {
            .message_timestamp("convert cpp_nomeR output to data.frame...")
        }

        if(!keepStartProb){
            predict_res_list <- predict_res_list[c("COVER_PROB","FOOTPRINT_CONF")]
        }

        return(lapply(predict_res_list,as.data.frame,
                      stringsAsFactors = FALSE,
                      check.names = FALSE))
    } else {
        stop("retrieved NULL results from C++ function.")
    }
}
