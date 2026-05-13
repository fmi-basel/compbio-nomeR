


#' Extract data as \code{data.table} object

#' @param sampleName Name of sample for which data needs to be extracted.
#' @param assayName Name of assay for which data needs to be extracted.
#' @param datColName Name of a data containing column in the output \code{data.table}.
#' @inheritParams predict_footprints_SE
#'
#' @returns \code{data.table} object with columns
#' \describe{
#'     \item{\code{seqnames}}{Chromosome name}
#'     \item{\code{refPos}}{Genomic coordinate}
#'     \item{\code{fragID}}{fragment name}
#'     \item{\code{assay_data} or column with the name defined by \code{datColName}}{data containing column}
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
                                  datColName = NULL){
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
