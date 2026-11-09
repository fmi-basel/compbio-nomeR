#' Create Adaptive Genomic Chunks from a BAM Index
#'
#' @description
#' Tiles the genome into variable-width chunks so that each chunk contains
#' approximately \code{target_reads} mapped reads, as estimated from the BAM
#' index without reading the BAM file itself.  This avoids the RAM spikes
#' caused by fixed-width chunking when some regions (e.g. blacklisted loci)
#' carry extremely high read density.
#'
#' \strong{Algorithm:}
#' \enumerate{
#'   \item The BAI linear index is parsed directly to obtain the compressed
#'     BGZF-block byte offset at every 16,384 bp interval (the native BAI
#'     resolution).  Differences between consecutive offsets proxy the data
#'     volume — and therefore read density — in each window.
#'   \item \code{\link[Rsamtools]{idxstatsBam}} supplies the true mapped-read
#'     count per chromosome, which calibrates the byte signal into estimated
#'     read counts.  Chromosomes with fewer than \code{min_reads_calib} reads
#'     use the genome-wide bytes-per-read ratio as a fallback.  This
#'     calibration naturally accounts for varying per-read payload sizes
#'     (e.g. the large MM/ML tags in SAMOSA / Fiber-seq BAMs).
#'   \item Windows are greedily merged into chunks: a chunk is closed once its
#'     estimated read count \eqn{\geq} \code{target_reads} \strong{and} its
#'     width \eqn{\geq} \code{min_chunk_width}.
#'   \item Leading and trailing empty (zero-byte) windows are trimmed from each
#'     chunk.  Chromosomes and chunks with no reads are dropped.
#' }
#'
#' @param bamfile Character scalar. Path to the BAM file.
#' @param target_reads Positive number. Desired number of mapped reads per
#'   chunk.
#' @param bam_index Character scalar. Path to the BAM index file (default:
#'   \code{paste0(bamfile, ".bai")}).  Must be a BAI-format index; CSI indexes
#'   are not supported.
#' @param min_chunk_width Positive integer. Minimum chunk width in base pairs
#'   (default: 500,000).  A chunk will not be closed until its width reaches
#'   this value even if \code{target_reads} has already been exceeded.  This
#'   prevents many tiny chunks in very high-density regions.  If a chromosome
#'   is shorter than \code{min_chunk_width}, the whole chromosome forms one
#'   chunk.
#' @param min_reads_calib Positive integer. Minimum number of mapped reads on a
#'   chromosome before a per-chromosome bytes-per-read ratio is computed
#'   (default: 1000).  References below this threshold use the genome-wide
#'   ratio.
#'
#' @return A sorted \code{\link[GenomicRanges]{GRanges}} with metadata columns:
#' \describe{
#'   \item{\code{chunk_id}}{Integer chunk index (1-based, consecutive).}
#'   \item{\code{est_reads}}{Estimated number of mapped reads in the chunk,
#'     derived from the BAI byte signal.}
#' }
#' Names follow the zero-padded format \code{"chunk_id_00001"} matching the
#' width of the total chunk count.  Chromosomes with no mapped reads are absent
#' from the result.
#'
#' @seealso \code{\link[Rsamtools]{idxstatsBam}}
#'
#' @importFrom Rsamtools idxstatsBam
#' @importFrom GenomicRanges GRanges sort
#' @importFrom IRanges IRanges
#' @export
create_adaptive_bam_chunks <- function(bamfile,
                                        target_reads,
                                        bam_index       = paste0(bamfile, ".bai"),
                                        min_chunk_width = 500000L,
                                        min_reads_calib = 1000L) {

    if (!file.exists(bamfile))  stop("BAM file not found: ", bamfile)
    if (!file.exists(bam_index)) stop("BAM index not found: ", bam_index)
    target_reads    <- as.double(target_reads)
    min_chunk_width <- as.integer(min_chunk_width)
    min_reads_calib <- as.integer(min_reads_calib)
    if (target_reads    < 1) stop("target_reads must be >= 1")
    if (min_chunk_width < 1L) stop("min_chunk_width must be >= 1")

    WIN <- 16384L  # BAI linear-index bin size (fixed by spec)

    ## Per-chromosome mapped-read counts from BAI (fast, no BAM scan).
    ## idxstatsBam includes a "*" row for unmapped reads; the BAI linear index
    ## has no such entry, so drop it before any length comparison.
    idx <- Rsamtools::idxstatsBam(bamfile, index = bam_index)
    idx <- droplevels(idx[idx$seqnames != "*", , drop = FALSE])

    ## Parse BAI linear index: one compressed-block byte offset per 16,384 bp
    lin <- .read_bai_linear(bam_index)

    if (length(lin) != nrow(idx))
        stop("BAI/header mismatch: index has ", length(lin),
             " references, BAM header has ", nrow(idx))

    ## Byte differences between consecutive compressed-block offsets proxy
    ## read density.  cummax() enforces monotonicity (coordinate-sorted BAM).
    n_win   <- lengths(lin)
    all_off <- cummax(c(unlist(lin, use.names = FALSE), file.size(bamfile)))
    all_dif <- diff(all_off)

    ## Split per chromosome. Use explicit factor levels so chromosomes with
    ## n_win == 0 (no linear-index entries) are preserved as empty elements,
    ## which prevents off-by-one indexing when some chromosomes are skipped.
    grp   <- factor(rep.int(seq_along(lin), n_win), levels = seq_along(lin))
    bytes <- split(all_dif, grp)

    ## Bytes-per-read calibration: per-chromosome where data allow, else genome-wide
    tot_bytes  <- vapply(bytes, sum, numeric(1L))
    bpr_global <- sum(tot_bytes) / max(sum(idx$mapped), 1L)
    bpr <- ifelse(idx$mapped >= min_reads_calib & tot_bytes > 0,
                  tot_bytes / pmax(idx$mapped, 1L),
                  bpr_global)

    min_win <- max(1L, as.integer(ceiling(min_chunk_width / WIN)))

    ## Build chunks per chromosome
    out_list <- lapply(seq_len(nrow(idx)), function(k) {
        nb <- bytes[[k]]
        if (!length(nb) || sum(nb) == 0) return(NULL)

        est <- nb / bpr[k]   # estimated reads per 16,384-bp window
        n   <- length(est)
        cs  <- cumsum(est)
        nz  <- nb > 0

        ## Greedy forward scan using binary search for efficiency
        s_win <- integer(n); e_win <- integer(n)
        base  <- 0; i <- 1L; m <- 0L
        while (i <= n) {
            ## First window where cumulative reads from `base` hit target_reads
            j <- findInterval(base + target_reads, cs, left.open = TRUE) + 1L
            ## Enforce minimum span; clamp to last window
            j    <- min(max(j, i + min_win - 1L), n)
            m    <- m + 1L
            s_win[m] <- i; e_win[m] <- j
            base <- cs[j]; i <- j + 1L
        }
        if (m == 0L) return(NULL)
        s_win <- s_win[seq_len(m)]; e_win <- e_win[seq_len(m)]

        ## Trim leading/trailing zero-byte windows per chunk; pad trimmed span
        ## back to min_win within the greedy allocation (s_win[p], e_win[p]).
        ## The terminal chunk of a chromosome may end up narrower than min_win
        ## when the BAI has no entries beyond the last read.
        keep <- logical(m)
        for (p in seq_len(m)) {
            w <- which(nz[s_win[p]:e_win[p]])
            if (!length(w)) next
            lo  <- s_win[p] + w[1L]        - 1L
            hi  <- s_win[p] + w[length(w)] - 1L
            gap <- min_win - (hi - lo + 1L)
            if (gap > 0L) {
                lo <- max(s_win[p], lo - ceiling(gap / 2L))
                hi <- min(e_win[p], hi + (min_win - (hi - lo + 1L)))
                lo <- max(s_win[p], hi - min_win + 1L)
            }
            s_win[p] <- lo; e_win[p] <- hi; keep[p] <- TRUE
        }
        s_win <- s_win[keep]; e_win <- e_win[keep]
        if (!length(s_win)) return(NULL)

        chr_len   <- as.integer(idx$seqlength[k])
        est_reads <- vapply(seq_along(s_win),
                            function(p) sum(est[s_win[p]:e_win[p]]),
                            numeric(1L))

        GenomicRanges::GRanges(
            seqnames  = idx$seqnames[k],
            ranges    = IRanges::IRanges(
                start = (s_win - 1L) * WIN + 1L,
                end   = pmin(e_win  *  WIN,  chr_len)
            ),
            est_reads = est_reads
        )
    })

    ## Combine, sort, assign zero-padded chunk IDs
    result <- Filter(Negate(is.null), out_list)
    if (length(result) == 0L)
        return(GenomicRanges::GRanges())

    chunks_gr <- do.call(c, result)
    chunks_gr <- GenomicRanges::sort(chunks_gr, ignore.strand = TRUE)

    n_chunks <- length(chunks_gr)
    chunks_gr$chunk_id <- seq_len(n_chunks)
    id_fmt <- paste0("chunk_id_%0", nchar(n_chunks), "d")
    names(chunks_gr) <- sprintf(id_fmt, seq_len(n_chunks))

    chunks_gr
}


## Internal: parse the BAI linear index.
## Returns a list of length n_ref; each element is a numeric vector of
## compressed BGZF-block byte offsets (one per 16,384 bp window), or
## numeric(0) for references with no linear-index entries.
.read_bai_linear <- function(index) {
    con <- file(index, "rb")
    on.exit(close(con))

    ## BAI magic: 0x42 0x41 0x49 0x01 ("BAI\1")
    magic <- readBin(con, "raw", 4L)
    expected <- as.raw(c(0x42, 0x41, 0x49, 0x01))
    if (!identical(magic, expected))
        stop("Not a BAI index file (CSI indexes are not supported): ", index)

    n_ref <- readBin(con, "integer", 1L, size = 4L)

    lapply(seq_len(n_ref), function(i) {
        n_bin <- readBin(con, "integer", 1L, size = 4L)
        ## Skip all bin records: each is 4 (bin id) + 4 (n_chunk) + 16*n_chunk bytes
        for (b in seq_len(n_bin)) {
            readBin(con, "integer", 1L, size = 4L)           # bin id (uint32)
            n_chunk <- readBin(con, "integer", 1L, size = 4L) # n_chunk (uint32)
            seek(con, 16L * n_chunk, origin = "current")      # skip chunk offsets
        }
        n_intv <- readBin(con, "integer", 1L, size = 4L)
        if (n_intv == 0L) return(numeric(0L))

        ## Each entry: uint64 little-endian BGZF virtual offset.
        ## Upper 48 bits (bytes 3–8, 1-indexed) = compressed block file offset.
        raw_bytes <- readBin(con, "raw", 8L * n_intv)
        r <- matrix(as.double(raw_bytes), nrow = 8L)
        colSums(r[3:8, , drop = FALSE] * (256 ^ (0:5)))
    })
}
