


#' Extract assay data from a SummarizedExperiment as a \code{data.table}
#'
#' @param sampleName Character scalar. Name of the sample to extract data for.
#' @param assayName Character scalar. Name of the assay to extract.
#' @param datColName Character scalar. Name to use for the data column in the
#'   output \code{data.table}. Defaults to \code{"assay_data"} if \code{NULL}.
#' @inheritParams predict_footprints_SE
#'
#' @returns A \code{data.table} with columns:
#' \describe{
#'     \item{\code{seqnames}}{Chromosome name.}
#'     \item{\code{refPos}}{Genomic coordinate.}
#'     \item{\code{fragID}}{Fragment (read) name.}
#'     \item{\code{assay_data}}{Data values (renamed to \code{datColName} if supplied).}
#'   }
#' @export
#' @importFrom SummarizedExperiment SummarizedExperiment rowData colData rowRanges
#'     assay assayNames assay<-
#' @importFrom SparseArray NaArray nnawhich
#' @importFrom GenomicRanges GPos match seqnames start end
#' @importFrom IRanges subsetByOverlaps
#' @importFrom S4Vectors DataFrame SimpleList metadata
#' @importFrom rlang .data
#' @import data.table
extract_assay_datatable <- function(se,
                                  sampleName,
                                  assayName,
                                  datColName = assayName){
    if(!sampleName %in% colnames(se)){
        stop("Can't find ",sampleName," in samples of the input se.")
    }
    if(!assayName %in% assayNames(se)){
        stop("Can't find ",sampleName," in assays of the input se.")
    }

    nNA_idxmat <- nnawhich(assay(se,assayName)[[sampleName]],arr.ind=TRUE)

    nNA_DT <- data.table(seqnames = as.vector(seqnames(rowRanges(se))[nNA_idxmat[,1]]),
                         refPos = start(rowRanges(se))[nNA_idxmat[,1]],
                         fragID = colnames(assay(se,assayName)[[sampleName]])[nNA_idxmat[,2]],
                         assay_data = assay(se,assayName)[[sampleName]][nNA_idxmat])
    if(!is.null(datColName)){
        setnames(nNA_DT,"assay_data",datColName)
    }
    return(nNA_DT)
}
