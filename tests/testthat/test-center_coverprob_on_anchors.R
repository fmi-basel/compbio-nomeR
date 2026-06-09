## ---- helpers -----------------------------------------------------------------

make_cover_dt <- function(seqnames = "chr1",
                          start    = 100L:120L,
                          strand   = "+",
                          sample   = "sA",
                          seed     = 1L) {
    set.seed(seed)
    data.table::data.table(
        seqnames = seqnames,
        start    = as.integer(start),
        strand   = strand,
        sample   = sample,
        fragID   = paste0("frag", seq_along(start)),
        score    = runif(length(start))
    )
}

make_anchors <- function(seqnames = "chr1",
                         start    = 110L,
                         end      = 110L,
                         strand   = "+",
                         ...) {
    GenomicRanges::GRanges(
        seqnames = seqnames,
        ranges   = IRanges::IRanges(start = as.integer(start),
                                    end   = as.integer(end)),
        strand   = strand,
        ...
    )
}


## ---- basic correctness -------------------------------------------------------

test_that("output has expected columns and no i.* artefacts", {
    dt  <- make_cover_dt()
    anc <- make_anchors()

    res <- center_coverprob_on_anchors(dt, anc, window = 5L)

    expect_s3_class(res, "data.table")
    # All original columns preserved
    expect_true(all(names(dt) %in% names(res)))
    # New anchor columns added
    expect_true(all(c("anchor_pos", "anchor_idx", "anchor_strand",
                      "dist_to_anchor", "fragID_anchorID") %in% names(res)))
    # No i.* columns from foverlaps internals
    expect_false(any(grepl("^i\\.", names(res))))
})

test_that("dist_to_anchor is start - anchor_pos for plus-strand anchor", {
    dt  <- make_cover_dt()            # positions 100:120
    anc <- make_anchors(start = 110L, end = 110L, strand = "+")

    res <- center_coverprob_on_anchors(dt, anc, window = 5L)
    data.table::setorder(res, start)

    # anchor_pos = 110, window = 5 => positions 105..115 retained
    expect_equal(res$start, 105L:115L)
    expect_equal(res$dist_to_anchor, -5L:5L)
})

test_that("anchor_pos = 'start' uses range start as reference", {
    dt  <- make_cover_dt()
    anc <- GenomicRanges::GRanges("chr1",
                                   IRanges::IRanges(start = 108L, end = 114L),
                                   strand = "+")

    res_start  <- center_coverprob_on_anchors(dt, anc, window = 3L,
                                              anchor_pos = "start")
    res_end    <- center_coverprob_on_anchors(dt, anc, window = 3L,
                                              anchor_pos = "end")
    res_center <- center_coverprob_on_anchors(dt, anc, window = 3L,
                                              anchor_pos = "center")

    expect_true(all(res_start$anchor_pos  == 108L))
    expect_true(all(res_end$anchor_pos    == 114L))
    expect_true(all(res_center$anchor_pos == 111L))   # (108+114) %/% 2 = 111
})

test_that("minus-strand anchor flips sign of dist_to_anchor", {
    dt  <- make_cover_dt()
    anc       <- make_anchors(strand = "+")  # plus-strand anchor at 110
    anc_minus <- make_anchors(strand = "-")  # minus-strand anchor at same position

    res_plus  <- center_coverprob_on_anchors(dt, anc,       window = 3L)
    res_minus <- center_coverprob_on_anchors(dt, anc_minus, window = 3L)

    # Sort both by start for comparison
    data.table::setorder(res_plus,  start)
    data.table::setorder(res_minus, start)

    # Same positions selected
    expect_equal(res_plus$start, res_minus$start)
    # Distances are negated
    expect_equal(res_minus$dist_to_anchor, -res_plus$dist_to_anchor)
})

test_that("ignore_strand = TRUE disables sign flip for minus-strand anchors", {
    dt  <- make_cover_dt()
    anc <- make_anchors(strand = "-")

    res_flip   <- center_coverprob_on_anchors(dt, anc, window = 3L,
                                              ignore_strand = FALSE)
    res_noflip <- center_coverprob_on_anchors(dt, anc, window = 3L,
                                              ignore_strand = TRUE)

    data.table::setorder(res_flip,   start)
    data.table::setorder(res_noflip, start)

    expect_equal(res_noflip$dist_to_anchor, -res_flip$dist_to_anchor)
})

test_that("positions outside window are excluded", {
    dt  <- make_cover_dt()   # 100:120
    anc <- make_anchors()    # anchor at 110

    res <- center_coverprob_on_anchors(dt, anc, window = 2L)
    expect_equal(sort(res$start), 108L:112L)
})

test_that("window = 0 keeps only exact-position matches", {
    dt  <- make_cover_dt()
    anc <- make_anchors(start = 110L, end = 110L)

    res <- center_coverprob_on_anchors(dt, anc, window = 0L)
    expect_equal(nrow(res), 1L)
    expect_equal(res$start, 110L)
    expect_equal(res$dist_to_anchor, 0L)
})


## ---- mult parameter ----------------------------------------------------------

test_that("mult = 'all' duplicates positions inside multiple anchor windows", {
    dt <- make_cover_dt(start = 100L:120L)
    # Two anchors whose windows overlap around position 110
    anc <- GenomicRanges::GRanges(
        seqnames = c("chr1", "chr1"),
        ranges   = IRanges::IRanges(start = c(108L, 112L), end = c(108L, 112L)),
        strand   = "+"
    )
    res_all   <- center_coverprob_on_anchors(dt, anc, window = 5L, mult = "all")
    res_first <- center_coverprob_on_anchors(dt, anc, window = 5L, mult = "first")

    # Positions 107-113 fall in both windows (overlap region 107-113 for anchor1,
    # 107-117 for anchor2, overlap at 107-113 = 7 positions)
    expect_true(nrow(res_all) > nrow(res_first))
    # mult = "first": each position appears at most once
    expect_equal(nrow(res_first), length(unique(res_first$start)))
})


## ---- anchor metadata --------------------------------------------------------

test_that("anchor mcols are propagated to output", {
    dt  <- make_cover_dt()
    anc <- make_anchors(score = 0.99, motif_name = "CTCF")

    res <- center_coverprob_on_anchors(dt, anc, window = 2L)

    expect_true("score" %in% names(res) || "motif_name" %in% names(res))
    expect_true(all(res$motif_name == "CTCF"))
})

test_that("fragID_anchorID uses anchor 'ID' column when present", {
    dt  <- make_cover_dt()
    anc <- make_anchors(ID = "CTCF_site_001")

    res <- center_coverprob_on_anchors(dt, anc, window = 2L)

    expect_true("fragID_anchorID" %in% names(res))
    expect_true(all(grepl(":CTCF_site_001$", res$fragID_anchorID)))
})

test_that("fragID_anchorID falls back to anchor_idx when no 'ID' column", {
    dt  <- make_cover_dt()
    anc <- make_anchors()   # no ID mcol

    res <- center_coverprob_on_anchors(dt, anc, window = 2L)

    expect_true(all(grepl(":1$", res$fragID_anchorID)))
})

test_that("anchor_idx correctly indexes into anchors when multiple anchors given", {
    dt  <- make_cover_dt(start = 100L:130L)
    anc <- GenomicRanges::GRanges(
        seqnames = c("chr1", "chr1"),
        ranges   = IRanges::IRanges(start = c(105L, 120L), end = c(105L, 120L)),
        strand   = "+"
    )
    res <- center_coverprob_on_anchors(dt, anc, window = 2L)

    expect_setequal(unique(res$anchor_idx), c(1L, 2L))
    expect_equal(res[anchor_idx == 1L, unique(anchor_pos)], 105L)
    expect_equal(res[anchor_idx == 2L, unique(anchor_pos)], 120L)
})


## ---- no-match behaviour -----------------------------------------------------

test_that("no matches returns zero-row result with a warning", {
    dt  <- make_cover_dt(seqnames = "chr1", start = 100L:120L)
    anc <- make_anchors(seqnames = "chr2")  # different chromosome

    expect_warning(
        res <- center_coverprob_on_anchors(dt, anc, window = 5L),
        regexp = "No positions fell"
    )
    expect_equal(nrow(res), 0L)
    expect_true(all(names(dt) %in% names(res)))
})


## ---- input validation --------------------------------------------------------

test_that("non-data.table cover_prob_dt triggers cli_abort", {
    anc <- make_anchors()
    expect_error(
        center_coverprob_on_anchors(as.data.frame(make_cover_dt()), anc),
        class = "rlang_error"
    )
})

test_that("empty cover_prob_dt triggers cli_abort", {
    dt  <- make_cover_dt()[integer(0L)]
    anc <- make_anchors()
    expect_error(
        center_coverprob_on_anchors(dt, anc),
        class = "rlang_error"
    )
})

test_that("missing required columns triggers cli_abort", {
    dt  <- make_cover_dt()[, .(seqnames, start)]  # fragID, strand, sample missing
    anc <- make_anchors()
    expect_error(
        center_coverprob_on_anchors(dt, anc),
        class = "rlang_error"
    )
})

test_that("non-GRanges anchors triggers cli_abort", {
    dt <- make_cover_dt()
    expect_error(
        center_coverprob_on_anchors(dt, data.frame(seqnames = "chr1",
                                                    start = 110L, end = 110L)),
        class = "rlang_error"
    )
})

test_that("empty GRanges triggers cli_abort", {
    dt  <- make_cover_dt()
    anc <- GenomicRanges::GRanges()
    expect_error(
        center_coverprob_on_anchors(dt, anc),
        class = "rlang_error"
    )
})

test_that("negative window triggers cli_abort", {
    dt  <- make_cover_dt()
    anc <- make_anchors()
    expect_error(
        center_coverprob_on_anchors(dt, anc, window = -1L),
        class = "rlang_error"
    )
})

test_that("non-finite window triggers cli_abort", {
    dt  <- make_cover_dt()
    anc <- make_anchors()
    expect_error(
        center_coverprob_on_anchors(dt, anc, window = Inf),
        class = "rlang_error"
    )
    expect_error(
        center_coverprob_on_anchors(dt, anc, window = NA_real_),
        class = "rlang_error"
    )
})

test_that("non-logical ignore_strand triggers cli_abort", {
    dt  <- make_cover_dt()
    anc <- make_anchors()
    expect_error(
        center_coverprob_on_anchors(dt, anc, ignore_strand = "yes"),
        class = "rlang_error"
    )
})
