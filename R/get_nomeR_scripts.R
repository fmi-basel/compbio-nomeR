#' Locate and Describe the nomeR Command-Line Scripts
#'
#' @description
#' Returns the full paths to the two independent Rscript-based workflows
#' bundled with nomeR and prints a short usage guide.
#'
#' The two scripts implement separate statistical models and can be used
#' independently or in combination:
#'
#' \describe{
#'   \item{\code{footprint_spectral_analysis-SAMOSA.R}}{
#'     \strong{Footprint Spectral Analysis (FSA).}
#'     Estimates the footprint length spectrum and emission probabilities
#'     from a BAM file.  The results are saved as a YAML file that can
#'     optionally be passed to the prediction script to derive footprint
#'     models.  Can also be used standalone to characterise the footprint
#'     size distribution of an experiment.}
#'   \item{\code{predict_footprints-SAMOSA.R}}{
#'     \strong{Genome-wide Footprint Prediction.}
#'     Predicts per-molecule footprints genome-wide and computes per-tile
#'     enrichment BigWig tracks.  Footprint models can be supplied either
#'     as an FSA YAML (from the script above) or as a pre-built model YAML
#'     (\code{--ftpmodelyaml}), making the script usable without running
#'     FSA first.}
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
#' paths <- get_nomeR_scripts(print_usage = FALSE)
#' paths
#'
#' @export
get_nomeR_scripts <- function(print_usage = TRUE) {

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

        cli::cli_h2("Footprint Spectral Analysis")
        cli::cli_text(
            "Bayesian model that estimates the footprint length spectrum and
            emission probabilities from a sample of single molecules.
            Results are saved as a YAML file and can be used to derive
            footprint models for the prediction script, or analysed on their own."
        )
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
            "--bamfile"        = "Input BAM file with 6mA modification probabilities (SAMOSA / Fiber-seq).",
            "--outfsayaml"     = "Output YAML with FSA results.",
            "--outpdf"         = "Optional PDF with the footprint length spectrum plot.",
            "--outrds"         = "Optional RDS with the full FSA result data frame.",
            "--correctseqbias" = "'no_correction' (default) or 'BC_KMF' for k-mer bias correction.",
            "--nfragfsa"       = "Number of randomly sampled fragments for FSA (default: 2000).",
            "--threads"        = "Number of parallel threads (default: 1).",
            "--libdir"         = "Colon-separated extra R library paths (e.g. for mcprogress)."
        ))

        cli::cli_h2("Genome-wide Footprint Prediction")
        cli::cli_text(
            "HMM-based statistical model for per-molecule footprint prediction
            and per-tile enrichment scoring across the whole genome.
            Footprint models can be supplied as an FSA YAML (--fsayaml) or
            as a pre-built model YAML (--ftpmodelyaml); neither is mandatory
            if the other is provided."
        )
        cli::cli_text("")
        cli::cli_code(paste(
            "Rscript", shQuote(paths[["predict_footprints"]]),
            "--bamfile    <input.bam>",
            "--fsayaml    <fsa.yaml>",
            "--outputdir  <output_directory>",
            "--threads    4"
        ))
        cli::cli_text("")
        cli::cli_text("{.strong Key options:}")
        cli::cli_dl(c(
            "--bamfile"        = "Input BAM file with 6mA modification probabilities (SAMOSA / Fiber-seq).",
            "--fsayaml"        = "FSA YAML used to derive footprint models (alternative to --ftpmodelyaml).",
            "--ftpmodelyaml"   = "Pre-built footprint model YAML (alternative to --fsayaml; takes priority).",
            "--outputdir"      = "Output directory (created automatically; default: nomeR_output/).",
            "--overwrite"      = "Overwrite existing output directory.",
            "--ftpdecoding"    = "'PosteriorDecoding' (default), 'PV', or 'Viterbi'.",
            "--chunksize"      = "Genome processed in chunks of this size in bp (default: 5,000,000).",
            "--tilewidthstep"  = "Sliding-window width and step for enrichment tracks, e.g. '500,250' (default).",
            "--correctseqbias" = "Sequence bias correction method.",
            "--noftps"         = "Skip BED12 footprint output.",
            "--noenrich"       = "Skip BigWig enrichment output.",
            "--threads"        = "Number of parallel threads (default: 1).",
            "--libdir"         = "Colon-separated extra R library paths."
        ))

        cli::cli_h2("Script paths")
        cli::cli_dl(base::setNames(paths, names(paths)))
    }

    invisible(paths)
}
