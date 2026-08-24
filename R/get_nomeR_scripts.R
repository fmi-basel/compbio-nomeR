#' Locate and Describe the nomeR Command-Line Scripts and Bundled Data
#'
#' @description
#' Returns the full paths to the Rscript-based workflows bundled with nomeR
#' and to three pre-computed reference data files, and prints a short usage
#' guide.
#'
#' \strong{Scripts:}
#' \describe{
#'   \item{\code{footprint_spectral_analysis-SAMOSA.R}}{
#'     \strong{Footprint Spectral Analysis (FSA).}
#'     Estimates the footprint length spectrum and emission probabilities from
#'     a SAMOSA / Fiber-seq BAM file.  Results are saved as a YAML that can be
#'     passed to the SAMOSA prediction script via \code{--fsayaml}.}
#'   \item{\code{predict_footprints-SAMOSA.R}}{
#'     \strong{Footprint prediction — SAMOSA / Fiber-seq.}
#'     Per-molecule footprint prediction and per-tile enrichment
#'     BigWig tracks for 6mA SAMOSA / Fiber-seq data.  Supports Bayesian
#'     k-mer sequence-bias correction (\code{--correctseqbias BetaCorrect} or
#'     \code{BetaUniform}).  Footprint models supplied via \code{--fsayaml} or
#'     \code{--ftpmodelyaml} (takes priority).}
#'   \item{\code{predict_footprints-NOMe-Seq.R}}{
#'     \strong{Footprint prediction — NOMe-seq.}
#'     Per-molecule footprint prediction and per-tile enrichment
#'     BigWig tracks for bisulfite GCH NOMe-seq data from QuasR- or
#'     Bismark-aligned BAM files.  Footprint models supplied via
#'     \code{--fsayaml} or \code{--ftpmodelyaml}.  No sequence-bias
#'     correction is currently available.}
#' }
#'
#' \strong{Bundled reference data (\code{inst/extdata}):}
#' \describe{
#'   \item{\code{mESC_Abdulhai_etal_SAMOSA_ftp_models_10_200_by_10_flat_prior.yaml}}{
#'     Pre-built SAMOSA / Fiber-seq footprint models (Abdulhay et al. 2023,
#'     mESC): footprints of 10, 20, \ldots, 200 bp with equal prior
#'     probabilities.  Pass to \code{--ftpmodelyaml}.}
#'   \item{\code{mESC_Guo_etal_2017_NOMe-Seq_ftp_models_10_200_by10_flat_prior.yaml}}{
#'     Pre-built NOMe-seq footprint models (Guo et al. 2017, mESC): footprints
#'     of 10, 20, \ldots, 200 bp with equal prior probabilities.  Pass to
#'     \code{--ftpmodelyaml}.}
#'   \item{\code{mESC_Abdulhay_etal_2023_SAMOSA_BetaCorrect_fit_kmer_7.txt}}{
#'     7-mer Beta-distribution parameters for Bayesian sequence-bias correction
#'     of SAMOSA data (Abdulhay et al. 2023, mESC).  Pass to
#'     \code{--betaparams} when using \code{--correctseqbias BetaCorrect}.}
#' }
#'
#' @param print_usage Logical. If \code{TRUE} (default) a short usage guide is
#'   printed to the console.
#'
#' @return A named list with two elements:
#' \describe{
#'   \item{\code{scripts}}{Named character vector of installed script paths:
#'     \code{footprint_spectral_analysis}, \code{predict_footprints_SAMOSA},
#'     \code{predict_footprints_NOMe_Seq}.}
#'   \item{\code{extdata}}{Named character vector of bundled data-file paths:
#'     \code{SAMOSA_footprint_models}, \code{NOMe_Seq_footprint_models},
#'     \code{SAMOSA_BetaCorrect_params}.}
#' }
#' Any file not found in the current installation is returned as
#' \code{NA_character_}.
#'
#' @examples
#' paths <- get_nomeR_scripts(print_usage = FALSE)
#' paths$scripts
#' paths$extdata
#'
#' @export
get_nomeR_scripts <- function(print_usage = TRUE) {

    ## scripts
    fsa_script  <- system.file("exec", "footprint_spectral_analysis-SAMOSA.R", package = "nomeR")
    pred_samosa <- system.file("exec", "predict_footprints-SAMOSA.R",           package = "nomeR")
    pred_nome   <- system.file("exec", "predict_footprints-NOMe-Seq.R",         package = "nomeR")

    ## bundled reference data
    samosa_models <- system.file("extdata",
                                  "mESC_Abdulhai_etal_SAMOSA_ftp_models_10_200_by_10_flat_prior.yaml",
                                  package = "nomeR")
    nome_models   <- system.file("extdata",
                                  "mESC_Guo_etal_2017_NOMe-Seq_ftp_models_10_200_by10_flat_prior.yaml",
                                  package = "nomeR")
    beta_params   <- system.file("extdata",
                                  "mESC_Abdulhay_etal_2023_SAMOSA_BetaCorrect_fit_kmer_7.txt",
                                  package = "nomeR")

    scripts <- c(
        footprint_spectral_analysis = if (nzchar(fsa_script))    fsa_script    else NA_character_,
        predict_footprints_SAMOSA   = if (nzchar(pred_samosa))   pred_samosa   else NA_character_,
        predict_footprints_NOMe_Seq = if (nzchar(pred_nome))     pred_nome     else NA_character_
    )

    extdata <- c(
        SAMOSA_footprint_models   = if (nzchar(samosa_models)) samosa_models else NA_character_,
        NOMe_Seq_footprint_models = if (nzchar(nome_models))   nome_models   else NA_character_,
        SAMOSA_BetaCorrect_params = if (nzchar(beta_params))   beta_params   else NA_character_
    )

    paths <- list(scripts = scripts, extdata = extdata)

    if (print_usage) {
        cli::cli_h1("nomeR command-line scripts")

        ## ── 1. Footprint Spectral Analysis ────────────────────────────────────
        cli::cli_h2("1. Footprint Spectral Analysis — SAMOSA / Fiber-seq")
        cli::cli_text(
            "Estimates the footprint length spectrum and emission probabilities
            from a random sample of single molecules. Results are saved as a
            YAML file that can be passed to the SAMOSA prediction script via
            {.code --fsayaml} or analysed independently."
        )
        cli::cli_text("")
        cli::cli_code(paste(
            "Rscript", shQuote(scripts[["footprint_spectral_analysis"]]),
            "--bamfile    <input.bam>",
            "--outfsayaml <output_fsa.yaml>",
            "--nfragfsa   2000",
            "--threads    4"
        ))
        cli::cli_text("")
        cli::cli_text("{.strong Key options:}")
        cli::cli_dl(c(
            "--bamfile"    = "Input BAM file with 6mA modification probabilities (SAMOSA / Fiber-seq).",
            "--outfsayaml" = "Output YAML with FSA results.",
            "--outpdf"     = "Optional PDF with the footprint length spectrum plot.",
            "--outrds"     = "Optional RDS with the full FSA result data frame.",
            "--nfragfsa"   = "Number of randomly sampled fragments for FSA (default: 2000).",
            "--threads"    = "Number of parallel threads (default: 1).",
            "--libdir"     = "Colon-separated extra R library paths."
        ))

        ## ── 2. SAMOSA Prediction ──────────────────────────────────────────────
        cli::cli_h2("2. Footprint Prediction — SAMOSA / Fiber-seq")
        cli::cli_text(
            "Per-molecule footprint prediction and per-tile enrichment
            BigWig tracks for 6mA SAMOSA / Fiber-seq data. Footprint models are
            supplied via {.code --fsayaml} or {.code --ftpmodelyaml} (takes
            priority). A bundled model YAML is available at
            {.code paths$extdata['SAMOSA_footprint_models']}. Optional Bayesian
            k-mer sequence-bias correction is supported."
        )
        cli::cli_text("")
        cli::cli_text("{.emph Basic usage (bundled models, no bias correction):}")
        cli::cli_code(paste(
            "Rscript", shQuote(scripts[["predict_footprints_SAMOSA"]]),
            "--bamfile      <input.bam>",
            "--ftpmodelyaml <footprint_models.yaml>",
            "--outputdir    <output_directory>",
            "--threads      4"
        ))
        cli::cli_text("")
        cli::cli_text("{.emph With Bayesian k-mer sequence-bias correction (BetaCorrect):}")
        cli::cli_code(paste(
            "Rscript", shQuote(scripts[["predict_footprints_SAMOSA"]]),
            "--bamfile        <input.bam>",
            "--ftpmodelyaml   <footprint_models.yaml>",
            "--correctseqbias BetaCorrect",
            "--betaparams     <beta_params.txt>",
            "--refseq         <reference.fasta>",
            "--outputdir      <output_directory>",
            "--threads        4"
        ))
        cli::cli_text("")
        cli::cli_text("{.strong Key options:}")
        cli::cli_dl(c(
            "--bamfile"            = "Input BAM file with 6mA modification probabilities.",
            "--fsayaml"            = "FSA YAML used to derive footprint models (alternative to --ftpmodelyaml).",
            "--ftpmodelyaml"       = "Pre-built footprint model YAML (takes priority over --fsayaml).",
            "--bgnoiselengths"     = "Length range for background noise footprints, e.g. '2:5:1' (default).",
            "--tflengths"          = "Length range for TF footprints, e.g. '20:50:10' (default).",
            "--nucllengths"        = "Length range for nucleosome footprints, e.g. '100:200:10' (default).",
            "--outputdir"          = "Output directory (created automatically; default: nomeR_output/).",
            "--overwrite"          = "Overwrite existing output directory.",
            "--ftpdecoding"        = "'PosteriorDecoding' (default), 'PV', or 'Viterbi'.",
            "--chunksize"          = "Genome processed in chunks of this size in bp (default: 5,000,000).",
            "--tilewidthstep"      = "Sliding-window width,step for enrichment tracks, e.g. '500,250' (default).",
            "--coverpsc"           = "Pseudo-count for posterior coverage log-transformation (default: 0.1).",
            "--correctseqbias"     = "Bias correction: 'no_correction' (default), 'BetaCorrect', or 'BetaUniform'.",
            "--betaparams"         = "[BetaCorrect] TXT file with Beta distribution parameters (see SAMOSA_BetaCorrect_params).",
            "--betaunifparams"     = "[BetaUniform] TXT file with Beta-Uniform mixture parameters.",
            "--refseq"             = "FASTA file with reference sequence (required for k-mer bias correction).",
            "--kmer"               = "Width of k-mer sequence context for bias correction (default: 7).",
            "--kmerblacklist"      = "TXT file with k-mers to exclude (default: bundled mESC blacklist).",
            "--filterkmerblacklist"= "Filter positions with blacklisted k-mer contexts (requires --kmerblacklist).",
            "--noftps"             = "Skip BED12 footprint output.",
            "--noenrich"           = "Skip BigWig enrichment output.",
            "--threads"            = "Number of parallel threads (default: 1).",
            "--libdir"             = "Colon-separated extra R library paths."
        ))

        ## ── 3. NOMe-Seq Prediction ────────────────────────────────────────────
        cli::cli_h2("3. Footprint Prediction — NOMe-seq")
        cli::cli_text(
            "Per-molecule footprint prediction and per-tile enrichment
            BigWig tracks for bisulfite GCH NOMe-seq data from QuasR- or
            Bismark-aligned BAM files. Footprint models are supplied via
            {.code --fsayaml} or {.code --ftpmodelyaml} (takes priority). A
            bundled model YAML is available at
            {.code paths$extdata['NOMe_Seq_footprint_models']}. Sequence-bias
            correction is not performed for NOMe-seq data."
        )
        cli::cli_text("")
        cli::cli_code(paste(
            "Rscript", shQuote(scripts[["predict_footprints_NOMe_Seq"]]),
            "--bamfile      <input.bam>",
            "--bamformat    QuasR",
            "--refgenome    <reference.fasta>",
            "--ftpmodelyaml <footprint_models.yaml>",
            "--outputdir    <output_directory>",
            "--threads      4"
        ))
        cli::cli_text("")
        cli::cli_text("{.strong Key options:}")
        cli::cli_dl(c(
            "--bamfile"        = "Input BAM file with bisulfite GCH methylation calls.",
            "--bamformat"      = "Aligner that produced the BAM: 'QuasR' or 'Bismark' (required).",
            "--refgenome"      = "FASTA file with the reference genome used for alignment (required).",
            "--fsayaml"        = "FSA YAML used to derive footprint models (alternative to --ftpmodelyaml).",
            "--ftpmodelyaml"   = "Pre-built footprint model YAML (takes priority over --fsayaml).",
            "--bgnoiselengths" = "Length range for background noise footprints, e.g. '2:5:1' (default).",
            "--tflengths"      = "Length range for TF footprints, e.g. '20:50:10' (default).",
            "--nucllengths"    = "Length range for nucleosome footprints, e.g. '100:200:10' (default).",
            "--outputdir"      = "Output directory (created automatically; default: nomeR_output/).",
            "--overwrite"      = "Overwrite existing output directory.",
            "--ftpdecoding"    = "'PosteriorDecoding' (default), 'PV', or 'Viterbi'.",
            "--chunksize"      = "Genome processed in chunks of this size in bp (default: 1,000,000).",
            "--tilewidthstep"  = "Sliding-window width,step for enrichment tracks, e.g. '500,250' (default).",
            "--coverpsc"       = "Pseudo-count for posterior coverage log-transformation (default: 0.1).",
            "--noftps"         = "Skip BED12 footprint output.",
            "--noenrich"       = "Skip BigWig enrichment output.",
            "--threads"        = "Number of parallel threads (default: 1).",
            "--libdir"         = "Colon-separated extra R library paths."
        ))

        ## ── Bundled reference data ────────────────────────────────────────────
        cli::cli_h2("Bundled reference data")
        cli::cli_text(
            "Three pre-computed reference files are shipped with nomeR and can
            be used as starting points without running a footprint spectral
            analysis first. Retrieve their paths from
            {.code paths$extdata}."
        )
        cli::cli_text("")
        cli::cli_dl(c(
            "SAMOSA_footprint_models"   =
                "SAMOSA/Fiber-seq footprint models, 10-200 bp flat prior
                (Abdulhay et al. 2023, mESC). Pass to --ftpmodelyaml.",
            "NOMe_Seq_footprint_models" =
                "NOMe-seq footprint models, 10-200 bp flat prior
                (Guo et al. 2017, mESC). Pass to --ftpmodelyaml.",
            "SAMOSA_BetaCorrect_params" =
                "7-mer Beta-distribution parameters for SAMOSA sequence-bias
                correction (Abdulhay et al. 2023, mESC). Pass to --betaparams
                with --correctseqbias BetaCorrect."
        ))
        cli::cli_text("")
        cli::cli_text("{.strong Paths:}")
        cli::cli_dl(extdata)

        ## ── Script paths ──────────────────────────────────────────────────────
        cli::cli_h2("Script paths")
        cli::cli_dl(scripts)
    }

    invisible(paths)
}
