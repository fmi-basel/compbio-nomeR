

#' Get shapes for beta distribution using method of moments
#'
#' @param x vector of modification probabilities
#'
#' @returns \code{list} with \code{n_dat} - number of values used,
#' \code{shape1} and \code{shape2} - shape parameters for beta distribution.
#' @noRd
#' @keywords internal
#' @importFrom stats var
.fit_beta_mom <- function(x) {
    m  <- mean(x)
    v  <- var(x)

    common <- m * (1 - m) / v - 1
    shape1 <- m * common
    shape2 <- (1 - m) * common

    list(n_dat = length(x), shape1 = shape1, shape2 = shape2)
}


#' Get shapes for beta distributions for SummarizedExperiment
#'
#' @param se \code{SummarizedExperiment} object
#' @param assayName name of assay
#' @param control_samplename name of the sample
#' @param min_n_data minimum number of points for sequence context
#'
#' @returns \code{data.table} containing shapes for beta distribution for each sequence context
#' @noRd
#' @keywords internal
#' @importFrom SummarizedExperiment SummarizedExperiment assay rowData colData
#'     colData<-
#' @importFrom SparseArray NaArray nnawhich
#' @importFrom GenomicRanges GPos match seqnames start end strand
#' @importFrom Seqinfo seqinfo
#' @importFrom IRanges subsetByOverlaps IRanges IRangesList
#' @importFrom S4Vectors DataFrame SimpleList metadata metadata<-
#'     make_zero_col_DFrame
#' @importFrom rlang .data
#' @import data.table
.get_beta_shapes <- function(se,assayName,control_samplename,min_n_data = 100){

    cnt_data <- assay(se, assayName)[,control_samplename]

    seqcont <- rowData(se)[,"sequenceContext"]

    ## get shapes for control
    ### get M-indices of non-NAs
    cnt_nonNA_data <- nnawhich(cnt_data, arr.ind = TRUE)
    ## count number of data points per sequence context add add any sequence contexts with less than min_n_data
    ## data points to groups other
    cnt_seq_ndat <- data.frame(table("sequenceContext" = as.character(seqcont[cnt_nonNA_data[,1]]))) %>%
        mutate(new_seq_cont = ifelse(.data$Freq > min_n_data,as.character(.data$sequenceContext),"other"))

    ## add columns with indices for redefined sequence contexts
    cnt_nonNA_data <- cbind(cnt_nonNA_data,
                               match(as.character(seqcont[cnt_nonNA_data[,1]]),cnt_seq_ndat$sequenceContext))
    seqmodprobs <- split(cnt_data[as.matrix(cnt_nonNA_data[,1:2])], ## mod_probs
                         cnt_seq_ndat$new_seq_cont[cnt_nonNA_data[,3]]  ## sequence context
    )
    shapes <- rbindlist(lapply(names(seqmodprobs),
                                  function(scnt){
                                      sh <- .fit_beta_mom(seqmodprobs[[scnt]])
                                      data.table(seqcont = scnt,
                                                 n_dat = sh$n_dat,
                                                 shape1 = sh$shape1,
                                                 shape2 = sh$shape2)
                                  }))
    return(shapes)

}


#' Fetch shapes of Beta distributions for positive and negative controls
#'
#' @param neg_control_sampleName sample name for negative control, i.e. SMF experiment without MTase
#' @param pos_control_sampleName sample name for positive control, i.e. SMF experiment with MTase for naked DNA
#' @param min_n_data minimum number of values. All sequence contexts with fewer data points will be merged into "other"
#' @inheritParams predict_footprints_SE
#'
#' @returns \code{list} with \code{data.table}'s containing shapes of Beta distributions for positive and negative controls
#' @export
#' @importFrom SummarizedExperiment SummarizedExperiment rowData colData
#'     assay assayNames
#' @importFrom SparseArray NaArray nnawhich
#' @importFrom GenomicRanges GPos match seqnames start end
#' @importFrom IRanges subsetByOverlaps
#' @importFrom S4Vectors DataFrame SimpleList metadata
#' @import data.table
get_SeqContext_control_beta_shapes_SE <- function(se,
                                       neg_control_sampleName,
                                       pos_control_sampleName,
                                       assayName = "mod_prob",
                                       min_n_data = 100){

    ## check if sequence context is present
    stopifnot("sequenceContext" %in% colnames(rowData(se)))
    stopifnot(all(c(neg_control_sampleName,pos_control_sampleName) %in% colnames(se)))

    ## get shapes for positive control
    negcnt_shapes <- .get_beta_shapes(se = se,
                                      assayName = assayName,
                                      control_samplename = neg_control_sampleName,
                                      min_n_data = min_n_data)
    poscnt_shapes <- .get_beta_shapes(se = se,
                                      assayName = assayName,
                                      control_samplename = pos_control_sampleName,
                                      min_n_data = min_n_data)

    list("positiveControlShapes" = poscnt_shapes,
         "negativeControlShapes" = negcnt_shapes)
}
