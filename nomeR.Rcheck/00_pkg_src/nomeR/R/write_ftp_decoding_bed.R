
#' Write footprint decoding to a BED12 file
#'
#' @description
#' Note that the molecule coordinates in the output BED12 file span the region
#' from the start of the leftmost to the end of the rightmost footprint and do
#' not correspond to the original genomic coordinates of the sequenced molecules.
#' @param ftp_decode_dt A \code{data.table} with footprint decoding results,
#'   as returned by \code{\link{predict_footprints_SE}} with \code{returnAs="data.table"}.
#' @param file Character scalar. Path to the output BED12 file.
#' @param ... Additional options passed to \code{\link[data.table]{fwrite}}.
#' @import data.table
#' @importFrom checkmate assertDataTable
#' @export
#'
write_ftp_decoding_bed <- function(ftp_decode_dt,
                                   file,
                                   ...){
    chrom <- chromStart <- chromEnd <- name <- fragID <- seqnames <- start <-
        score <- strand <- thickStart <- thickEnd <- itemRgb <- blockCount <- blockSizes <-
        blockStarts <- NULL # due to NSE notes in R CMD check

    assertDataTable(x = ftp_decode_dt,min.cols = 5)
    if(!all(c("seqnames","start","width","strand","fragID") %in% colnames(ftp_decode_dt))){
        stop("ftp_decode_dt must contain columns \"seqnames\",\"start\",\"width\",\"strand\",\"fragID\"")
    }

    setorder(ftp_decode_dt, fragID, seqnames,start)

    bed12_dt <- ftp_decode_dt[, {
        chromStart <- min(start) - 1
        chromEnd   <- max(start + width - 1)

        list(
            chrom       = as.character(seqnames[1]),
            chromStart  = chromStart,
            chromEnd    = chromEnd,
            name        = fragID[1],
            score       = 0L,                       # BED-required, can be 0
            strand      = as.character(strand[1]),
            thickStart  = chromStart,               # usually same as chromStart
            thickEnd    = chromEnd,
            itemRgb     = "0",                      # or "R,G,B"
            blockCount  = .N,
            blockSizes  = paste(width, collapse = ","),
            blockStarts = paste(start - 1L - chromStart, collapse = ",")
        )
    }, by = fragID][,fragID := NULL]

    fwrite(
        bed12_dt,
        file = file,
        sep = "\t",
        col.names = FALSE,
        scipen = 9999,
        ...
    )




}
