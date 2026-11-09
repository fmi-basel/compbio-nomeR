## Tests for create_adaptive_bam_chunks() and .read_bai_linear()
##
## Requires the mESC SAMOSA test BAM included in the package repository.
## The BAM has 2000 reads mapped to chr1 (mm10, 195,471,971 bp).

## Absolute path: this BAM lives in the dev source tree, not in the installed
## package. skip_if_no_bam() skips gracefully in R CMD check or any other
## environment where the file is absent.
BAM <- "/tachyon/groups/scratch/gpeters/ozonevge/Programming/RPackages/compbio-nomeR/TEST_MODEL/19_test_clinomeR_scripts/mESC-SAMOSA-test.bam"
BAI <- paste0(BAM, ".bai")

## Skip all tests gracefully when the test BAM is absent (e.g. installed package
## without the TEST_MODEL directory).
skip_if_no_bam <- function() {
    testthat::skip_if(!file.exists(BAM), "Test BAM not found; skipping adaptive-chunking tests.")
}


## ── .read_bai_linear ──────────────────────────────────────────────────────────

test_that(".read_bai_linear returns a list of the correct length", {
    skip_if_no_bam()
    lin <- nomeR:::.read_bai_linear(BAI)
    idx <- Rsamtools::idxstatsBam(BAM)
    ## BAI has one entry per reference; the '*' pseudo-reference is NOT in the
    ## linear index, so the list length equals the number of real references.
    n_real <- sum(idx$seqnames != "*")
    expect_type(lin, "list")
    expect_length(lin, n_real)
})

test_that(".read_bai_linear entries are numeric vectors", {
    skip_if_no_bam()
    lin <- nomeR:::.read_bai_linear(BAI)
    for (i in seq_along(lin)) {
        expect_true(is.numeric(lin[[i]]),
                    info = paste("reference", i, "is not numeric"))
    }
})

test_that(".read_bai_linear aborts on a non-BAI file", {
    skip_if_no_bam()
    expect_error(nomeR:::.read_bai_linear(BAM), "Not a BAI index file")
})


## ── create_adaptive_bam_chunks: error handling ────────────────────────────────

test_that("create_adaptive_bam_chunks errors on missing BAM", {
    expect_error(
        create_adaptive_bam_chunks("/nonexistent/file.bam", target_reads = 500),
        "BAM file not found"
    )
})

test_that("create_adaptive_bam_chunks errors on missing BAI", {
    skip_if_no_bam()
    expect_error(
        create_adaptive_bam_chunks(BAM,
                                   target_reads = 500,
                                   bam_index    = "/nonexistent.bai"),
        "BAM index not found"
    )
})

test_that("create_adaptive_bam_chunks errors on target_reads < 1", {
    skip_if_no_bam()
    expect_error(
        create_adaptive_bam_chunks(BAM, target_reads = 0, bam_index = BAI),
        "target_reads must be >= 1"
    )
})

test_that("create_adaptive_bam_chunks errors on min_chunk_width < 1", {
    skip_if_no_bam()
    expect_error(
        create_adaptive_bam_chunks(BAM, target_reads = 500,
                                   bam_index = BAI, min_chunk_width = 0L),
        "min_chunk_width must be >= 1"
    )
})


## ── create_adaptive_bam_chunks: output structure ──────────────────────────────

test_that("create_adaptive_bam_chunks returns a GRanges", {
    skip_if_no_bam()
    res <- create_adaptive_bam_chunks(BAM, target_reads = 500, bam_index = BAI)
    expect_s4_class(res, "GRanges")
})

test_that("result has chunk_id (integer) and est_reads (numeric) metadata", {
    skip_if_no_bam()
    res <- create_adaptive_bam_chunks(BAM, target_reads = 500, bam_index = BAI)
    expect_true("chunk_id"  %in% names(GenomicRanges::mcols(res)))
    expect_true("est_reads" %in% names(GenomicRanges::mcols(res)))
    expect_type(res$chunk_id,  "integer")
    expect_type(res$est_reads, "double")
})

test_that("chunk_id values are 1:n consecutive integers", {
    skip_if_no_bam()
    res <- create_adaptive_bam_chunks(BAM, target_reads = 500, bam_index = BAI)
    n   <- length(res)
    expect_identical(res$chunk_id, seq_len(n))
})

test_that("names are zero-padded chunk_id strings", {
    skip_if_no_bam()
    res <- create_adaptive_bam_chunks(BAM, target_reads = 500, bam_index = BAI)
    n   <- length(res)
    expected_names <- sprintf(paste0("chunk_id_%0", nchar(n), "d"), seq_len(n))
    expect_identical(names(res), expected_names)
})

test_that("seqinfo contains only real chromosomes (no '*' level)", {
    skip_if_no_bam()
    res <- create_adaptive_bam_chunks(BAM, target_reads = 500, bam_index = BAI)
    expect_false("*" %in% GenomeInfoDb::seqlevels(res))
})


## ── create_adaptive_bam_chunks: genomic correctness ──────────────────────────

test_that("chunks are sorted by genomic coordinate", {
    skip_if_no_bam()
    res    <- create_adaptive_bam_chunks(BAM, target_reads = 500, bam_index = BAI)
    sorted <- GenomicRanges::sort(res, ignore.strand = TRUE)
    expect_identical(as.character(GenomicRanges::seqnames(res)),
                     as.character(GenomicRanges::seqnames(sorted)))
    expect_identical(BiocGenerics::start(res), BiocGenerics::start(sorted))
})

test_that("chunks are non-overlapping within a chromosome", {
    skip_if_no_bam()
    res <- create_adaptive_bam_chunks(BAM, target_reads = 500, bam_index = BAI)
    ol  <- GenomicRanges::findOverlaps(res, res, ignore.strand = TRUE)
    ## only self-hits should remain
    self_hits <- S4Vectors::queryHits(ol) == S4Vectors::subjectHits(ol)
    expect_true(all(self_hits),
                info = "Chunks overlap each other — they should be non-overlapping.")
})

test_that("chunk starts are >= 1 and ends are within chromosome length", {
    skip_if_no_bam()
    res    <- create_adaptive_bam_chunks(BAM, target_reads = 500, bam_index = BAI)
    idx    <- Rsamtools::idxstatsBam(BAM)
    idx    <- idx[idx$seqnames != "*", , drop = FALSE]
    len_by_chr <- setNames(idx$seqlength, as.character(idx$seqnames))
    for (chr in unique(as.character(GenomicRanges::seqnames(res)))) {
        sub <- res[as.character(GenomicRanges::seqnames(res)) == chr]
        expect_true(all(BiocGenerics::start(sub) >= 1L))
        expect_true(all(BiocGenerics::end(sub) <= len_by_chr[[chr]]),
                    info = paste("chunk end exceeds seqlength on", chr))
    }
})

test_that("non-terminal chunks per chromosome are at least min_chunk_width bp wide", {
    skip_if_no_bam()
    ## The terminal chunk on each chromosome receives whatever windows remain
    ## after the greedy split; it may be narrower than min_chunk_width when the
    ## BAI has no entries beyond the last mapped read.  Only intermediate chunks
    ## are guaranteed to satisfy the minimum width constraint.
    min_w <- 100000L
    res   <- create_adaptive_bam_chunks(BAM, target_reads = 500,
                                        bam_index = BAI, min_chunk_width = min_w)
    chrs  <- as.character(GenomicRanges::seqnames(res))
    for (chr in unique(chrs)) {
        idx_chr <- which(chrs == chr)
        if (length(idx_chr) < 2L) next  # single chunk: no intermediate chunks
        non_terminal <- res[idx_chr[-length(idx_chr)]]
        widths <- BiocGenerics::width(non_terminal)
        expect_true(all(widths >= min_w),
                    info = sprintf("%s: non-terminal chunk width %d < min_chunk_width %d",
                                   chr, min(widths), min_w))
    }
})

test_that("est_reads are positive for all returned chunks", {
    skip_if_no_bam()
    res <- create_adaptive_bam_chunks(BAM, target_reads = 500, bam_index = BAI)
    expect_true(all(res$est_reads > 0))
})

test_that("total est_reads approximately equals total mapped reads", {
    skip_if_no_bam()
    res       <- create_adaptive_bam_chunks(BAM, target_reads = 500, bam_index = BAI)
    total_est <- sum(res$est_reads)
    idx       <- Rsamtools::idxstatsBam(BAM)
    total_map <- sum(idx$mapped)
    ## Calibration should bring the sum within 5% of the true count
    expect_lt(abs(total_est - total_map) / total_map, 0.05,
              label = sprintf("sum(est_reads)=%.0f vs mapped=%d (>5%% off)",
                              total_est, total_map))
})


## ── target_reads controls chunk count ────────────────────────────────────────

test_that("large target_reads produces fewer chunks than small target_reads", {
    skip_if_no_bam()
    res_small <- create_adaptive_bam_chunks(BAM, target_reads = 200,  bam_index = BAI)
    res_large <- create_adaptive_bam_chunks(BAM, target_reads = 2000, bam_index = BAI)
    expect_lt(length(res_large), length(res_small))
})

test_that("target_reads >= total reads gives exactly one chunk", {
    skip_if_no_bam()
    idx    <- Rsamtools::idxstatsBam(BAM)
    n_map  <- sum(idx$mapped)
    res    <- create_adaptive_bam_chunks(BAM, target_reads = n_map * 2L,
                                         bam_index = BAI)
    expect_length(res, 1L)
})
