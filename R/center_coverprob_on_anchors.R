## Silence R CMD check 'no visible binding' NOTEs for data.table column
## references inside := expressions.
utils::globalVariables(c("dist_to_anchor", "fragID", "fragID_anchorID"))

#' Center coverage probability profiles on genomic anchor positions
#'
#' @description
#' Given a coverage probability data.table produced by
#' \code{\link{predict_footprints_SE}} (the \code{COVER_PROB} element of the
#' \code{returnAs = "data.table"} output), and a set of genomic anchor
#' positions, this function centers each per-molecule coverage profile relative
#' to the nearest (or overlapping) anchor.  Positions are re-expressed as
#' signed distances from the anchor, enabling aggregation and meta-profile
#' analysis across many loci.
#'
#' @param cover_prob_dt A \code{data.table} as returned in the
#'   \code{COVER_PROB} slot of \code{\link{predict_footprints_SE}} with
#'   \code{returnAs = "data.table"}.  Required columns:
#'   \code{seqnames} (chromosome), \code{start} (1-based genomic position),
#'   \code{strand}, \code{sample}, \code{fragID}, and one or more footprint
#'   probability columns (e.g. \code{Nucl}, \code{TF}, \code{background}).
#' @param anchors A \code{\link[GenomicRanges]{GRanges}} object defining the
#'   anchor loci.  Each range is reduced to its centre (mid-point) unless
#'   \code{anchor_pos} is \code{"start"} or \code{"end"}.  The object may
#'   carry arbitrary metadata columns (e.g. motif name, score); these are
#'   propagated to the output.
#' @param window Single finite non-negative integer scalar.  Half-width of the
#'   window (in base pairs) around each anchor to retain.  Positions with
#'   \code{|start - anchor| > window} are dropped.  Default: \code{1000L}.
#' @param anchor_pos Character scalar specifying which part of each anchor
#'   range to use as the reference point.  One of \code{"center"} (default),
#'   \code{"start"}, or \code{"end"}.
#' @param ignore_strand Logical.  If \code{FALSE} (default), positions on
#'   minus-strand anchors are sign-flipped so that upstream is always
#'   negative and downstream is always positive.  Set to \code{TRUE} to
#'   ignore strand and always report raw distance.
#' @param mult Character scalar passed to \code{\link[data.table]{foverlaps}}.
#'   \code{"all"} (default) returns one row per (position, anchor) pair, so a
#'   position inside two anchor windows appears twice.  \code{"first"} keeps
#'   only the first overlapping anchor per position, avoiding row duplication
#'   at the cost of discarding additional matches.
#'
#' @return A \code{data.table} with the same columns as \code{cover_prob_dt}
#'   plus:
#'   \describe{
#'     \item{\code{anchor_pos}}{Integer. Genomic coordinate of the anchor
#'       reference point.}
#'     \item{\code{anchor_idx}}{Integer. Index of the matched anchor in
#'       \code{anchors} (1-based).}
#'     \item{\code{dist_to_anchor}}{Integer. Signed distance from the
#'       position to the anchor (\code{start - anchor_pos}, flipped for
#'       minus-strand anchors when \code{ignore_strand = FALSE}).}
#'     \item{\code{fragID_anchorID}}{Character. Composite key formed as
#'       \code{fragID:anchor_idx} (or \code{fragID:ID} when the \code{anchors}
#'       GRanges carries an \code{ID} metadata column).  Useful for uniquely
#'       identifying each (molecule, anchor) pair.}
#'   }
#'   Rows that do not fall within \code{window} bp of any anchor are
#'   excluded.  If a position is within range of multiple anchors it is
#'   duplicated once per anchor (unless \code{mult = "first"}).
#'
#' @seealso \code{\link{predict_footprints_SE}}
#'
#' @importFrom GenomicRanges GRanges seqnames start end strand
#' @importFrom S4Vectors mcols
#' @importFrom methods is
#' @import data.table
#'
#' @export
center_coverprob_on_anchors <- function(cover_prob_dt,
                                        anchors,
                                        window        = 1000L,
                                        anchor_pos    = c("center", "start", "end"),
                                        ignore_strand = FALSE,
                                        mult          = c("all", "first")) {

    anchor_pos <- match.arg(anchor_pos)
    mult       <- match.arg(mult)

    # ---- input validation -----------------------------------------------

    if (!is.data.table(cover_prob_dt))
        cli::cli_abort("{.arg cover_prob_dt} must be a {.cls data.table}.")
    if (nrow(cover_prob_dt) == 0L)
        cli::cli_abort("{.arg cover_prob_dt} has no rows.")
    req_cols <- c("seqnames", "start", "strand", "sample", "fragID")
    missing_cols <- setdiff(req_cols, names(cover_prob_dt))
    if (length(missing_cols) > 0L)
        cli::cli_abort(
            "{.arg cover_prob_dt} is missing required column(s): {.val {missing_cols}}.")

    if (!methods::is(anchors, "GRanges"))
        cli::cli_abort("{.arg anchors} must be a {.cls GRanges} object.")
    if (length(anchors) == 0L)
        cli::cli_abort("{.arg anchors} must contain at least one range.")

    if (!is.numeric(window) || length(window) != 1L ||
            !is.finite(window) || window < 0L)
        cli::cli_abort(
            "{.arg window} must be a single finite non-negative integer scalar.")
    window <- as.integer(window)

    if (!is.logical(ignore_strand) || length(ignore_strand) != 1L ||
            is.na(ignore_strand))
        cli::cli_abort(
            "{.arg ignore_strand} must be a single non-missing logical value.")

    # ---- compute anchor reference coordinates ---------------------------

    anc_coord <- switch(anchor_pos,
        center = as.integer(
            (GenomicRanges::start(anchors) + GenomicRanges::end(anchors)) %/% 2L),
        start  = as.integer(GenomicRanges::start(anchors)),
        end    = as.integer(GenomicRanges::end(anchors))
    )

    # ---- convert anchors to data.table with window intervals ------------

    anchor_dt <- data.table(
        seqnames      = as.character(GenomicRanges::seqnames(anchors)),
        start         = anc_coord - window,
        end           = anc_coord + window,
        anchor_pos    = anc_coord,
        anchor_idx    = seq_along(anchors),
        anchor_strand = as.character(GenomicRanges::strand(anchors))
    )
    ## propagate metadata columns carried by anchors
    anc_meta <- as.data.frame(S4Vectors::mcols(anchors))
    if (ncol(anc_meta) > 0L)
        anchor_dt <- cbind(anchor_dt, as.data.table(anc_meta))

    setkey(anchor_dt, seqnames, start, end)

    # ---- build a minimal index table for foverlaps ----------------------
    ## Avoids copying the large cover_prob_dt. We extract only the three key
    ## columns plus a row-index (.row_idx), which is cheap. foverlaps runs on
    ## this small table; results are then used to subscript back into the
    ## original cover_prob_dt by row index, leaving the input untouched.
    ## Each position is a single base (point interval), so end = start.

    pos_idx_dt <- cover_prob_dt[, .(
        seqnames  = as.character(seqnames),
        start     = as.integer(start),
        end       = as.integer(start),
        .row_idx  = .I
    )]
    setkey(pos_idx_dt, seqnames, start, end)

    # ---- overlap positions with anchor windows --------------------------

    overlaps_idx <- foverlaps(pos_idx_dt, anchor_dt,
                              by.x    = c("seqnames", "start", "end"),
                              by.y    = c("seqnames", "start", "end"),
                              type    = "within",
                              mult    = mult,
                              nomatch = NULL)

    if (nrow(overlaps_idx) == 0L) {
        cli::cli_warn("No positions fell within {window} bp of any anchor.")
        return(cover_prob_dt[integer(0L)])
    }

    # ---- assemble result ---------------------------------------------------
    ## Subscript matched rows from the original (large) data.table by row index.
    ## This is cheaper than a full copy: only matched rows are materialised.

    result_dt <- cover_prob_dt[overlaps_idx$.row_idx]

    # ---- compute signed distance to anchor --------------------------------

    result_dt[, dist_to_anchor := as.integer(start) - overlaps_idx$anchor_pos]

    if (!ignore_strand) {
        minus_anchor <- overlaps_idx$anchor_strand == "-"
        if (any(minus_anchor))
            result_dt[minus_anchor,
                      dist_to_anchor := -dist_to_anchor]
    }

    # ---- attach anchor metadata -------------------------------------------
    ## foverlaps prefixes x-side key columns with "i." in the output
    ## (i.seqnames, i.start, i.end).  Exclude those together with the y-side
    ## key columns (window boundaries) and the internal row index.

    anc_cols_keep <- setdiff(
        names(overlaps_idx),
        c("seqnames", "start", "end", ".row_idx",
          grep("^i\\.", names(overlaps_idx), value = TRUE))
    )
    result_dt[, (anc_cols_keep) := overlaps_idx[, .SD, .SDcols = anc_cols_keep]]

    # ---- composite (molecule, anchor) key ---------------------------------
    ## Use anchors' own ID column when available, otherwise fall back to
    ## anchor_idx.  Mirrors the fragID:ID pattern used in mot_dt-based
    ## foverlaps workflows.

    anchor_id_col <- if ("ID" %in% names(result_dt)) "ID" else "anchor_idx"
    result_dt[, fragID_anchorID := paste0(fragID, ":", get(anchor_id_col))]

    result_dt
}
