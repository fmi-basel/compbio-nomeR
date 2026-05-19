#' Locate and Describe the nomeR Command-Line Scripts
#'
#' @description
#' Returns the full paths to the two Rscript-based command-line workflows
#' bundled with nomeR and prints a short usage guide.
#'
#' \describe{
#'   \item{\code{footprint_spectral_analysis-SAMOSA.R}}{
#'     Estimates footprint length spectra and emission probabilities from a
#'     BAM file.  Run this first to generate the FSA YAML required by the
#'     prediction script.}
#'   \item{\code{predict_footprints-SAMOSA.R}}{
#'     Genome-wide footprint prediction and enrichment scoring.  Takes the
#'     FSA YAML produced above and a BAM file, and writes BED12 footprints
#'     and BigWig enrichment tracks.}
#' }
#'
#' @param print_usage Logical. If \code{TRUE} (default) a short usage guide
#'   is printed to the console.
#'
#' @return A named character vector with elements \code{footprint_spectral_analysis}
#'   and \code{predict_footprints}, each containing the absolute path to the
#'   corresponding installed script.  Returns \code{NA} for any script that
#'   cannot be found in the current installation.
#'
#' @examples
#' paths <- nomeR_scripts(print_usage = FALSE)
#' paths
#'
#' @export
nomeR_scripts <- function(print_usage = TRUE) {

    fsa_script  <- system.file("exec", "footprint_spectral_analysis-SAMOSA.R",
                               package = "nomeR")
    pred_script <- system.file("exec", "predict_footprints-SAMOSA.R",
                               package = "nomeR")

    paths <- c(
        footprint_spectral_analysis = if (nzchar(fsa_script))  fsa_script  else NA_character_,
        predict_footprints          = if (nzchar(pred_script)) pred_script else NA_character_
    )

    if (print_usage) {
        cli::cli_h1("nomeR command-line scripts")

        cli::cli_h2("Step 1 — Footprint Spectral Analysis")
        cli::cli_text("Estimates the footprint length spectrum and emission probabilities
        from a sample of reads. Run once per experiment / condition.")
        cli::cli_text("")
        cli::cli_code(paste(
            "Rscript", shQuote(paths[["footprint_spectral_analysis"]]),
            "--bamfile    <input.bam>",
            "--outfsayaml <output_fsa.yaml>",
            "--nfragfsa   2000",
            "--threads    4"
        ))
        cli::cli_text("")
        cli::cli_text("{.strong Key options:}")
        cli::cli_dl(c(
            "--bamfile"       = "Input BAM file with 6mA modification probabilities (SAMOSA / Fiber-seq).",
            "--outfsayaml"    = "Output YAML with FSA results (required by predict script).",
            "--outpdf"        = "Optional PDF with the footprint length spectrum plot.",
            "--correctseqbias"= "'no_correction' (default) or 'BC_KMF' for k-mer bias correction.",
            "--nfragfsa"      = "Number of randomly sampled fragments for FSA (default: 2000).",
            "--threads"       = "Number of parallel threads (default: 1).",
            "--libdir"        = "Colon-separated extra R library paths (e.g. for mcprogress)."
        ))

        cli::cli_h2("Step 2 — Genome-wide Footprint Prediction")
        cli::cli_text("Uses the FSA YAML from Step 1 to predict footprints genome-wide
        and compute per-tile enrichment BigWig tracks.")
        cli::cli_text("")
        cli::cli_code(paste(
            "Rscript", shQuote(paths[["predict_footprints"]]),
            "--bamfile    <input.bam>",
            "--fsayaml    <output_fsa.yaml>",
            "--outputdir  <output_directory>",
            "--threads    4"
        ))
        cli::cli_text("")
        cli::cli_text("{.strong Key options:}")
        cli::cli_dl(c(
            "--bamfile"       = "Input BAM file (same as Step 1).",
            "--fsayaml"       = "FSA YAML produced in Step 1.",
            "--outputdir"     = "Output directory (created automatically; default: nomeR_output/).",
            "--overwrite"     = "Overwrite existing output directory.",
            "--ftpmodeltype"  = "'fast' (default), 'medium', or 'slow' — controls footprint model resolution.",
            "--ftpdecoding"   = "'PosteriorDecoding' (default), 'PV', or 'Viterbi'.",
            "--chunksize"     = "Genome is processed in chunks of this size in bp (default: 5,000,000).",
            "--tilewidthstep" = "Sliding-window width and step for enrichment tracks, e.g. '500,250' (default).",
            "--correctseqbias"= "Sequence bias correction method (must match the setting used in Step 1).",
            "--noftps"        = "Skip BED12 footprint output.",
            "--noenrich"      = "Skip BigWig enrichment output.",
            "--threads"       = "Number of parallel threads (default: 1).",
            "--libdir"        = "Colon-separated extra R library paths."
        ))

        cli::cli_h2("Script paths")
        cli::cli_dl(setNames(paths, names(paths)))
    }

    invisible(paths)
}
