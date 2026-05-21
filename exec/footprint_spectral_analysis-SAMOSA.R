#!/usr/bin/env Rscript
## R script (v6) for footprint spectral analysis in SAMOSA/FiberSeq generated BAM file

## Resolve the lib directory this script was installed into so that the nomeR
## version loaded is always the one bundled with this script.
## Installed layout: <lib>/nomeR/exec/<this-script>.R
## Three dirname() calls climb: script -> exec/ -> nomeR/ -> <lib>/
local({
    argv <- commandArgs(trailingOnly = FALSE)
    f    <- sub("--file=", "", grep("--file=", argv, value = TRUE))
    if (length(f) == 1L) {
        pkg_lib <- dirname(dirname(dirname(normalizePath(f, mustWork = TRUE))))
        .libPaths(c(pkg_lib, .libPaths()))
    }
})


library(optparse)

option_list <- list(
    ### input options
    make_option(c("-b", "--bamfile"),
                type="character",
                help="BAM file containing 6mA modification probabilities for SAMOSA/FiberSeq data"),
    make_option(c("-m", "--thresholdmod"),
                type="numeric",
                default = 0.5,
                help="Threshold for modification probability to binarize the data into accessible and protected positions.
                Probabilities equal or higher will be considered accessible, otherwise protected. [default %default]"),

    ### options for correction of sequence bias
    make_option(c("--correctseqbias"),
                type="character",
                default = "no_correction",
                help="Method for correcting sequence biases. Can be 'no_correction', 'BC_KMF' - bayesian correction of modification probabilities followed by filtering of non-informative k-mers. [default %default]"),
    make_option(c("--negbetas"),
                type="character",
                default = system.file("extdata",
                                      "SAMOSA_mESC_negativeControl_betaShapes_kmer_7.txt",
                                      package = "nomeR"),
                help="path to TXT file containing shapes for beta distribution inferred from negative controls.
							Used for correction of modification probabilities. [default: bundled SAMOSA mESC (Abdulhay et al, 2023) shapes in %default]"),
    make_option(c("--posbetas"),
                type="character",
                default = system.file("extdata",
                                      "SAMOSA_mESC_positiveControl_betaShapes_kmer_7.txt",
                                      package = "nomeR"),
                help="path to TXT file containing shapes for beta distribution inferred from positive controls.
							Used for correction of modification probabilities. [default: bundled SAMOSA mESC (Abdulhay et al, 2023) shapes in %default]"),
    make_option(c("--refseq"),
                type="character",
                default = NULL,
                help="path to fasta file containing reference sequence"),
    make_option(c("--kmer"),
                type="integer",
                default = 7,
                help="width of sequence context for correction. [default %default]"),
    make_option(c("--kmerblacklist"),
                type="character",
                default = system.file("extdata",
                                      "SAMOSA_mESC_blacklist_kmer_7_cutoff_0.2.txt",
                                      package = "nomeR"),
                help="path to TXT file containing k-mers to ignore due to their strong sequence biases. [default: bundled SAMOSA mESC (Abdulhay et al, 2023) blacklist in %default]"),
    make_option(c("--quantnorm"),
                type="logical",
                action="store_true",
                default=FALSE,
                help="Perform quantile normalization of modification probabilities to match distribution of uncorrected probabilities. [default: FALSE]"),

    ### output options
    make_option(c("-s", "--outfsayaml"),
                type="character",
                help="[OPTIONAL] Path to an YAML file to save results of footprint spectrum analysis."),
    make_option(c("-p", "--outpdf"),
                type="character",
                help="[OPTIONAL] Path to an PDF file to save footprint spectrum plot."),
    make_option(c("-o", "--outrds"),
                type="character",
                help="[OPTIONAL] Path to an RDS file to save results of footprint spectrum analysis."),

    ### running options
    make_option(c("-n","--nfragfsa"),
                type = "integer",
                default = 2000,
                help="Number of random fragments to sample for footprint spectral analysis. [default %default]"
    ),
    make_option(c("-t", "--threads"),
                type="integer",
                default=1,
                help="Number of threads [default %default]"),
    make_option(c("-v", "--verbose"),
                action="store_true",
                default=FALSE,
                help="Print extra output [default]"),

    ### library options
    make_option(c("--libdir"),
                type="character",
                default=NULL,
                help="Colon-separated path(s) to additional R library directories to search for packages
                (e.g. '/path/to/libs' or '/path/a:/path/b'). Added before the default library search path. [default: none]")
)


### parse and validate options
opt <- parse_args(OptionParser(option_list = option_list))

if (!is.null(opt$libdir)) {
    extra_libs <- trimws(strsplit(opt$libdir, ":", fixed = TRUE)[[1]])
    extra_libs <- extra_libs[nchar(extra_libs) > 0 & dir.exists(extra_libs)]
    if (length(extra_libs) > 0)
        .libPaths(c(extra_libs, .libPaths()))
}

#### LOAD LIBRARIES ####
if(opt$verbose)
    cli::cli_h1("Loading R libraries")

suppressPackageStartupMessages({
    library(data.table)
    library(dplyr)
    library(stringr)
    library(SingleMoleculeGenomicsIO)
    library(SummarizedExperiment)
    library(Seqinfo)
    library(GenomicRanges)
    library(SparseArray)
    library(BiocParallel)
    library(yaml)
    library(Rsamtools)
    library(ggplot2)
    library(parallel)
    library(nomeR)    
    library(Biostrings)
})


#### VALIDATE INPUT PARAMETERS ####
##### check input bam file #####
if (is.null(opt$bamfile)) {
    cli::cli_abort("Input BAM file is required\n")
} else{
    if(!file.exists(opt$bamfile)){
        cli::cli_abort("Couldn't find input BAM file: ",opt$bamfile,"\n")
    }
}

##### check input model options #####

if(opt$thresholdmod <= 0 || opt$thresholdmod >= 1){
    cli::cli_abort("--thresholdmod must be within (0,1). Provided --thresholdmod {.val opt$thresholdmod}")
}


##### check parameters for sequence bias correction #####
if(!opt$correctseqbias %in% c("no_correction","BC_KMF")){
    cli::cli_abort("--correctseqbias allowed to be only 'no_correction' or 'BC_KMF'.")
}

if(opt$correctseqbias != "no_correction"){
    if(!is.null(opt$negbetas) && !file.exists(opt$negbetas))
        cli::cli_abort("Couldn't find file specified by --negbetas {opt$negbetas}")
    if(!is.null(opt$posbetas) && !file.exists(opt$posbetas))
        cli::cli_abort("Couldn't find file specified by --posbetas {opt$posbetas}")
    if(!is.null(opt$refseq) && !file.exists(opt$refseq))
        cli::cli_abort("Couldn't find file specified by --refseq {opt$refseq}")
    if(!is.null(opt$kmerblacklist) && !file.exists(opt$kmerblacklist))
        cli::cli_abort("Couldn't find file specified by --kmerblacklist {opt$kmerblacklist}")
    if(opt$kmer <= 0)
        cli::cli_abort("Incorrect parameter --kmer: {opt$kmer}")
}

##### check output options #####
if(is.null(opt$outfsayaml) && is.null(opt$outpdf) && is.null(opt$outrds)){
    cli::cli_abort("No output files specified.
                   At least one of the output files set by --outfsayaml, --outpdf, --outrds must be specified. ")
}


if(opt$threads < 1){
    cli::cli_abort("--threads can't be lower than 1. Provided --threads {opt$threads}")
}

if(opt$nfragfsa < 1000){
    cli::cli_abort("--nfragfsa can't be lower than 1000. Provided --nfragfsa {opt$nfragfsa}")
}


cli::cli_h1("Start analysis with the following options:")
cli::cli_dl(opt)
cli::cli_h1("")

#### load data for correction of sequence bias ###
if(opt$correctseqbias != "no_correction"){
    if(!is.null(opt$negbetas) && !is.null(opt$posbetas)){
        negcontrol_shapes <- data.table(read.table(opt$negbetas,
                                                   header = F,
                                                   col.names = c("seqcont",
                                                                 "n_dat",
                                                                 "shape1",
                                                                 "shape2")))
        poscontrol_shapes <- data.table(read.table(opt$posbetas,
                                                   header = F,
                                                   col.names = c("seqcont",
                                                                 "n_dat",
                                                                 "shape1",
                                                                 "shape2")))
    } else{
        negcontrol_shapes <- NULL
        poscontrol_shapes <- NULL
    }
} else{
    negcontrol_shapes <- NULL
    poscontrol_shapes <- NULL
}
## filter k-mers if a blacklist is provided
if(!is.null(opt$kmerblacklist)){
    ## load blacklist
    kmer_blacklist <- read.table(opt$kmerblacklist,header=F,col.names = c("kmer"))
} else{
    kmer_blacklist <- NULL
}

## by default assayName is mod_prob. but if sequence bias correction was done it should be changed to mod_prob_corrected
assayName <- "mod_prob"
if(opt$correctseqbias != "no_correction" && !is.null(negcontrol_shapes) && !is.null(poscontrol_shapes))
    assayName <- "mod_prob_corrected"


#### PERFORM FOOTPRINT SPECTRUM ANALYSIS ####
##### construct seqinfo object from bam header ######
bam_header <- Rsamtools::scanBamHeader(files = opt$bamfile,what = c("targets"))
seqinfo_bam <- Seqinfo::Seqinfo(seqnames = names(bam_header[[1]]$targets),
                                bam_header[[1]]$targets)



longest_chrom <- Seqinfo::seqnames(seqinfo_bam)[which.max(Seqinfo::seqlengths(seqinfo_bam))]
if(opt$correctseqbias == "no_correction"){
    cli::cli_progress_step("Loading modification probabilities for random fragments from BAM file")
    opt$kmer <- 0L
} else{
    cli::cli_progress_step("Loading modification probabilities and sequence contexts for random fragments from BAM file")
}
se <- SingleMoleculeGenomicsIO::readModBam(bamfiles = opt$bamfile,
                                           nAlnsToSample = opt$nfragfsa,
                                           seqnamesToSampleFrom = longest_chrom,
                                           modbase = "a",
                                           sequenceContextWidth = opt$kmer,
                                           sequenceReference = opt$refseq,
                                           BPPARAM = BiocParallel::MulticoreParam(workers = opt$threads),
                                           verbose = F)

rownames(se) <- seq_along(se)

# ----- 2.1 (optional) correction of sequence bias ------
if(opt$correctseqbias != "no_correction"){
    cli::cli_progress_step("Correction of sequence biases")
    if(!is.null(negcontrol_shapes) && !is.null(poscontrol_shapes)){
        cli::cli_inform("Bayesian correction of sequence biases")
        se <- nomeR::correct_modprob_SE(se,
                                 neg_control_shapes = negcontrol_shapes,
                                 pos_control_shapes = poscontrol_shapes,
                                 qnorm_to_raw = opt$quantnorm)

    }
    if(!is.null(kmer_blacklist)){
        cli::cli_inform("Filtering blacklisted kmers")
        kmer_blacklist <- Biostrings::DNAStringSet(kmer_blacklist[["kmer"]])
        ## remove k-mers
        keep_rows <- which(!(rowData(se)[,"sequenceContext"] %in% kmer_blacklist))
        se <- se[keep_rows, ]
    }
}

cli::cli_progress_step("Footprint Spectral Analysis (FSA)")

## parameters for footprint spectral analysis
ftp_bg_model <- "informative_prior"
bg_model_params <- list(bg_protect_prob_fixed = 0.05, bg_protect_min = 0.01,
                        bg_protect_max = 0.49, bg_protect_mean = 0.1, bg_protect_totcount = 1000)
ftp_model_params <- list(ftp_protect_prob_fixed = 0.95, ftp_protect_min = 0.51,
                         ftp_protect_max = 0.99, ftp_protect_mean = 0.9, ftp_protect_totcount = 1000)


fsa_data <- nomeR::ftp_spectral_analysis_SE(se = se,
                                            assayName = assayName,
                                            threshMod = opt$thresholdmod,
                                            ftp_lengths = 2:200,
                                            ftp_bg_model = ftp_bg_model,
                                            bg_model_params = bg_model_params,
                                            ftp_model_params = ftp_model_params,
                                            ncpu = opt$threads)
if(!is.null(opt$outrds)){
    cli::cli_progress_step("Saving inference results in {.file {opt$outrds}}")
    saveRDS(fsa_data,
            file = opt$outrds,
            compress = TRUE)
}


if(!is.null(opt$outpdf)){
    ## save plot for ftp spectrum
    cli::cli_progress_step("Saving footprint spectrum in {.file {opt$outpdf}}")
    ftp_spec_plot <- nomeR::plot_ftp_spectra_DF(fsa_data)

    ggplot2::ggsave(filename = opt$outpdf,
                    plot = ftp_spec_plot,
                    width=6,height=5.5
    )
}

if(!is.null(opt$outfsayaml)){
    cli::cli_progress_step("Saving FSA results into YAML file  {.file {opt$outfsayaml}}")
    fsa2exp <- list("ftp_spectrum" = fsa_data$ftp_spectrum[[1]] %>% 
                                        dplyr::select(ftp_length,mean) %>% 
                                        dplyr::filter(ftp_length > 1),
                    "bgcoverprior" = fsa_data$bg_coverage_mean[[1]],
                    "bgprotectprob" = fsa_data$bg_emis_mean[[1]],
                    "ftpprotectprob" = fsa_data$ftp_emis_mean[[1]])

    yaml::write_yaml(fsa2exp,
                     file = opt$outfsayaml)

}

cli::cli_progress_done()


cli::cli_alert_success("Footprint spectral analysis prediction has been finished!")
cli::cli_h1("The following output files have been stored:")
cli::cli_ol()
cli::cli_li("{.field Footprint spectral analysis results:}")
ulid <- cli::cli_ul()
if(!is.null(opt$outrds) && file.exists(opt$outrds))
    cli::cli_li("RDS file with estimates: {.file {opt$outrds}}")
if(!is.null(opt$outpdf) && file.exists(opt$outpdf))
    cli::cli_li("PDF file with plot for footprint spectrum: {.file {opt$outpdf}}")
if(!is.null(opt$outfsayaml) && file.exists(opt$outfsayaml))
    cli::cli_li("YAML file with results of Footprint Spectral Analysis: {.file {opt$outfsayaml}}")
cli::cli_end(ulid)

closeAllConnections()
