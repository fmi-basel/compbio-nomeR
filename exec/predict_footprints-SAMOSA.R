#!/usr/bin/env Rscript
## R script (v7) for predicting footprints in SAMOSA/FiberSeq generated BAM file
## BAM file is expected to have MM and ML tags

## Enzyme: EcoGII or Hia5
## Informative position: Adenines
## Accessibility decoding: Methylated - accessible, Unmethylated  - protected


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
    #### input options ####
    make_option(c("-b", "--bamfile"),
                type="character",
                help="BAM file containing 6mA modification probabilities for SAMOSA/FiberSeq data"),
    make_option(c("-s", "--fsayaml"),
                type="character",
                help="[OPTIONAL] Path to an YAML file containing results of footprint spectrum analysis.
                If neither of --fsayaml --ftpmodelyaml are not provided, the script will abort with an error.
                If --ftpmodelyaml is provided the --fsayaml will be ignored.
                If --ftpmodelyaml is absent and --fsayaml is provided the footprint models will be constructed for
                TF and Nucl based on footprint spectrum from the --fsayaml file."),
    make_option(c("-f", "--ftpmodelyaml"),
                type="character",
                help="[OPTIONAL] Path to an YAML file containing footprint models.
                If not provided, the footprint models will be constructed using footprint spectrum in --fsayaml."),
    make_option(c("-r", "--regionbed"),
                type="character",
                help="[OPTIONAL] BED file with region of interest. If not provided the analysis will be performed genome-wide."),

    #### options specifying model ####
    make_option(c("--bgnoiselengths"),
                type="character",
                default="2:5:1",
                help="Length ranges and steps for background noise footprints (format: \"min:max:step\", e.g. \"2:5:1\" for lengths 2,3,4,5). [default %default]
                Set to \"none\" to skip background noise model.
                If --ftpmodelyaml is set the --bgnoiselengths has no effect."),
    make_option(c("--tflengths"),
                type="character",
                default="20:50:10",
                help="Length ranges and steps for transcription factor footprints (format: \"min:max:step\", e.g. \"20:50:10\" for lengths 20,30,40,50). [default %default]
                Set to \"none\" to skip transcription factor model.
                If --ftpmodelyaml is set the --bgnoiselengths has no effect."),
    make_option(c("--nucllengths"),
                type="character",
                default="100:200:10",
                help="Length ranges and steps for nucleosome footprints (format: \"min:max:step\", e.g. \"100:200:10\" for lengths 100,110,120,...,200). [default %default]
                Set to \"none\" to skip nucleosome model.
                If --ftpmodelyaml is set the --bgnoiselengths has no effect."),
    make_option(c("-d", "--ftpdecoding"),
                type="character",
                default="PosteriorDecoding",
                help="Method for footprint decoding. Can be 'PV' - Posterior Viterbi; 'PosteriorDecoding' or 'Viterbi' [default %default]"),
    make_option(c("-a", "--tilewidthstep"),
                type="character",
                default = "500,250",
                help="Tile width and step (format: \"width,step\") for aggregating and calculating p-values for TF and BG. [default %default]"),
    make_option(c("--coverpsc"),
                type="numeric",
                default = "0.1",
                help="Pseudo-count for stabilization of coverage probabilities [default %default]"),


    #### options for correction of sequence bias ####
    make_option(c("--correctseqbias"),
                type="character",
                default = "no_correction",
                help="Method for correcting sequence biases. Can be 'no_correction',
                'BetaCorrect' - Bayesian correction using Beta distributions (requires --betaparams),
                'BetaUniform' - [!!!EXPERIMENTAL!!!] Bayesian correction using Beta-Uniform mixtures (requires --betaunifparams; isotonic regression applied by default, disable with --noisotonic). [default %default]"),
    make_option(c("--betaparams"),
                type="character",
                default = NULL,
                help="[BetaCorrect] Path to TXT file with Beta distribution parameters for positive and negative controls, as produced by get_SeqContext_control_beta_shapes_SE() (columns: seqcont n_pos alpha_pos beta_pos n_neg alpha_neg beta_neg)."),
    make_option(c("--betaunifparams"),
                type="character",
                default = NULL,
                help="[BetaUniform] Path to TXT file with Beta-Uniform mixture parameters produced by fit_SeqContext_BetaUnif_params_SE() (columns: seqcont alpha_pos beta_pos eps_pos alpha_neg beta_neg eps_neg n_pos n_neg converged mlrp_ok)."),
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
                default = NULL,
                help="path to TXT file containing k-mers to ignore due to their strong sequence biases. [default: bundled SAMOSA mESC (Abdulhay et al, 2023) blacklist in %default]"),
    make_option(c("--filterkmerblacklist"),
                action="store_true",
                default=FALSE,
                help="Filter positions with blacklisted k-mer sequence contexts. Requires --kmerblacklist. Can be used independently of --correctseqbias. [default %default]"),
    make_option(c("--quantnorm"),
                type="logical",
                action="store_true",
                default=TRUE,
                help="Perform quantile normalization of modification probabilities to match distribution of uncorrected probabilities. [default: TRUE]"),
    make_option(c("--noisotonic"),
                type="logical",
                action="store_true",
                default=FALSE,
                help="[BetaUniform] Disable pool-adjacent-violators (PAV) isotonic regression that enforces monotonicity between raw and corrected probabilities. By default isotonic regression is applied. [default: FALSE]"),
    #### output options ####
    make_option(c("-o", "--outputdir"),
                default = "nomeR_output/",
                type="character",
                help="Folder to store output files. [default %default]"),
    make_option(c("-z", "--tempdir"),
                #default = "nomeR_output/tempfolder",
                type="character",
                help="Folder to store temporary files. [default <output_dir>/tempfolder]"),
    make_option(c("-w", "--noenrich"),
                action="store_true",
                default=FALSE,
                help="Skip writing BIGWIG files containing log-enrichments neg log10(Pval) for TF-covered and accessible positions in genomic tiles? [default %default]"),
    make_option(c("-p", "--noftps"),
                action="store_true",
                default=FALSE,
                help="Skip writing output BED12 files containing footprints for Nucl, TF and Linkers? [default %default]"),

    #### running options ####
    make_option(c("-c","--chunksize"),
                type = "integer",
                default = 5000000,
                help="For genome-wide analyses, the footprint prediction is carried out consecutively for genomic regions of size --chunksize. [default %default]"
    ),
    make_option(c("-t", "--threads"),
                type="integer",
                default=1,
                help="Number of threads [default %default]"),
    make_option(c("-v", "--verbose"),
                action="store_true",
                default=FALSE,
                help="Print extra output [default]"),
    make_option(c("-e", "--overwrite"),
                action="store_true",
                default=FALSE,
                help="Overwrite existing output folder if already exists?
                WARNING: If the folder defined by --outputdir exists and -e is set all existing files will be deleted. [default %default]"),

    #### library options ####
    make_option(c("--libdir"),
                type="character",
                default=NULL,
                help="Colon-separated path(s) to additional R library directories to search for packages
                (e.g. '/path/to/libs' or '/path/a:/path/b'). Added before the default library search path. [default: none]")
)


### parse and validate options
opt <- parse_args(OptionParser(option_list = option_list))

if(is.null(opt$tempdir))
    opt$tempdir <- file.path(opt$outputdir,"tempfolder")

if (!is.null(opt$libdir)) {
    extra_libs <- trimws(strsplit(opt$libdir, ":", fixed = TRUE)[[1]])
    extra_libs <- extra_libs[nchar(extra_libs) > 0 & dir.exists(extra_libs)]
    if (length(extra_libs) > 0)
        .libPaths(c(extra_libs,.libPaths()))
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
    library(rtracklayer)
    library(Rsamtools)
    library(ggplot2)
    library(parallel)
    library(nomeR)
    library(mcprogress)
    library(Biostrings)
})

#### VALIDATE INPUT PARAMETERS ####

##### check input bam file #####
if (is.null(opt$bamfile)) {
    cli::cli_abort("Input BAM file is required\n")
} else{
    if(!file.exists(opt$bamfile)){
        cli::cli_abort("Couldn't find input BAM file: {.file {opt$bamfile}}")
    }
}

##### check input FSA/footprintmodels files #####
if(is.null(opt$fsayaml) && is.null(opt$ftpmodelyaml)){
    cli::cli_abort("Footprint models are defined by --fsayaml or --ftpmodelyaml. Either of these files must be provided.")
}

if(!is.null(opt$fsayaml) && !is.null(opt$ftpmodelyaml)){
    cli::cli_warn("Both --fsayaml or --ftpmodelyaml are provided. --ftpmodelyaml has priority and user-provided footprint models are used.")
}

if(!is.null(opt$fsayaml) && !file.exists(opt$fsayaml)){
    cli::cli_abort("Couldn't find file specified by --fsayaml {opt$fsayaml}")
}
if(!is.null(opt$ftpmodelyaml) && !file.exists(opt$ftpmodelyaml)){
    cli::cli_abort("Couldn't find file specified by --ftpmodelyaml {opt$ftpmodelyaml}")
}

if(!opt$ftpdecoding %in% c("PV","PosteriorDecoding","Viterbi")){
    cli::cli_abort("ftpdecoding must be one of the following: \"PV\",\"PosteriorDecoding\",\"Viterbi\"")
}

##### check input tile width options #####
if(!grepl(pattern = "\\d+,\\d+",opt$tilewidthstep)){
    cli::cli_abort("--tilewidthstep must be of format \"int,int\", e.g. \"500,250\". Provided: --tilewidthstep {opt$tilewidthstep}")
} else{
    strprs <- str_split_fixed(opt$tilewidthstep,pattern = ",",n=2)
    tile_width <- as.numeric(strprs[1])
    tile_step <- as.numeric(strprs[2])
}

##### check parameters for sequence bias correction #####
if(!opt$correctseqbias %in% c("no_correction","BetaCorrect","BetaUniform")){
    cli::cli_abort("--correctseqbias must be one of 'no_correction', 'BetaCorrect', 'BetaUniform'.")
}

if(opt$correctseqbias == "BetaCorrect"){
    if(is.null(opt$betaparams))
        cli::cli_abort("--correctseqbias BetaCorrect requires --betaparams.")
    if(!file.exists(opt$betaparams))
        cli::cli_abort("Couldn't find file specified by --betaparams {opt$betaparams}")
    if(!is.null(opt$refseq) && !file.exists(opt$refseq))
        cli::cli_abort("Couldn't find file specified by --refseq {opt$refseq}")
    if(opt$kmer <= 0)
        cli::cli_abort("Incorrect parameter --kmer: {opt$kmer}")
}

if(opt$correctseqbias == "BetaUniform"){
    if(is.null(opt$betaunifparams))
        cli::cli_abort("--correctseqbias BetaUniform requires --betaunifparams.")
    if(!file.exists(opt$betaunifparams))
        cli::cli_abort("Couldn't find file specified by --betaunifparams {opt$betaunifparams}")
    if(!is.null(opt$refseq) && !file.exists(opt$refseq))
        cli::cli_abort("Couldn't find file specified by --refseq {opt$refseq}")
    if(opt$kmer <= 0)
        cli::cli_abort("Incorrect parameter --kmer: {opt$kmer}")
}

if(opt$filterkmerblacklist){
    if(is.null(opt$kmerblacklist))
        cli::cli_abort("--filterkmerblacklist requires --kmerblacklist to be specified.")
    if(!is.null(opt$kmerblacklist) && !file.exists(opt$kmerblacklist))
        cli::cli_abort("Couldn't find file specified by --kmerblacklist {opt$kmerblacklist}")
    if(opt$kmer <= 0)
        cli::cli_abort("Incorrect parameter --kmer: {opt$kmer}")
}

##### check output options #####
if(!dir.exists(opt$outputdir)){
    dir.create(opt$outputdir)
} else if(dir.exists(opt$outputdir) && !opt$overwrite){
    cli::cli_abort("Output directory {opt$outputdir} already exists. Use --overwrite (-e) to overwrite existing files.")
} else if(dir.exists(opt$outputdir) && opt$overwrite){
    cli::cli_alert_warning("Output directory {opt$outputdir} already exists. All existing files in this directory will be deleted.")
    unlink(opt$outputdir,recursive=TRUE)
    dir.create(opt$outputdir)
}

if(opt$noenrich && opt$noftps){
    cli::cli_abort("Nothing to output as all options, i.e. --noenrich (-w) --noftps (-p) are set.")
}


##### check chunk size option #####
if(opt$chunksize <= 100000){
    cli::cli_abort("--chunksize can't be lower than 100,000. Provided --chunksize {opt$chunksize}")
}
if(opt$threads < 1){
    cli::cli_abort("--threads can't be lower than 1. Provided --threads {opt$threads}")
}


cli::cli_h1("Start analysis with the following options:")
cli::cli_dl(opt)
cli::cli_h1("")

#####  DEFINE FOOTPRINT MODELS #####
## if --ftpmodelyaml is defined we load user-provided footprint models from yaml file
## if not, we load the footprint spectrum and construct TF and Nucl footprints from it
## we have to define variables ftp_models and ftp_group_names
if(!is.null(opt$ftpmodelyaml)){
    ## load footprint model from a YAML file
    ftp_models <- yaml::read_yaml(opt$ftpmodelyaml)
    ## check if all fields are present
    if(is.null(ftp_models$ftp_models) ||
       is.null(ftp_models$bgprotectprob) ||
       is.null(ftp_models$bgcoverprior)){
        cli::cli_abort("Can't find one of \"ftp_models\" \"bgprotectprob\" \"bgcoverprior\". All three must be present in the input YAML: {opt$ftpmodelyaml}")
    }

    ## retrieve names of groups for which we will calculate enrichments etc.
    ftp_group_names <- unique(sapply(ftp_models$ftp_models,function(x){x$GROUP}))
} else if(!is.null(opt$fsayaml)){

    ## load footprint spectrum from a YAML file
    fsa_data <- yaml::read_yaml(opt$fsayaml)

    ## by default we only use TF and Nucl footprint groups
    ## define footprint lengths for background noise models
    if(!is.null(opt$bgnoiselengths) && opt$bgnoiselengths != "none"){
        bg_len_strs <- str_split_fixed(opt$bgnoiselengths,pattern = ":",n=3)
        bg_len_mat <- matrix(as.numeric(bg_len_strs),nrow=1)
        if(any(is.na(bg_len_mat))){
            cli::cli_abort("Couldn't parse --bgnoiselengths. --bgnoiselengths must be of format \"min:max:step\", e.g. \"2:5:1\" for lengths 2,3,4,5. Provided --bgnoiselengths {opt$bgnoiselengths}")
        }
        colnames(bg_len_mat) <- c("min_ftp_len","max_ftp_len","by")
        row.names(bg_len_mat) <- "background"
    } else{
        bg_len_mat <- NULL
    }
    ## define footprint lengths for TF models
    if(!is.null(opt$tflengths) && opt$tflengths != "none"){
        tf_len_strs <- str_split_fixed(opt$tflengths,pattern = ":",n=3)
        tf_len_mat <- matrix(as.numeric(tf_len_strs),nrow=1)
        if(any(is.na(tf_len_mat))){
            cli::cli_abort("Couldn't parse --tflengths. --tflengths must be of format \"min:max:step\", e.g. \"20:50:10\" for lengths 20,30,40,50. Provided --tflengths {opt$tflengths}")
        }
        colnames(tf_len_mat) <- c("min_ftp_len","max_ftp_len","by")
        row.names(tf_len_mat) <- "TF"
    } else{
        tf_len_mat <- NULL
    }
    ## define footprint lengths for Nucl models
    if(!is.null(opt$nucllengths) && opt$nucllengths != "none"){
        nucl_len_strs <- str_split_fixed(opt$nucllengths,pattern = ":",n=3)
        nucl_len_mat <- matrix(as.numeric(nucl_len_strs),nrow=1)
        if(any(is.na(nucl_len_mat))){
            cli::cli_abort("Couldn't parse --nucllengths. --nucllengths must be of format \"min:max:step\", e.g. \"100:200:10\" for lengths 100,110,120,...,200. Provided --nucllengths {opt$nucllengths}")
        }
        colnames(nucl_len_mat) <- c("min_ftp_len","max_ftp_len","by")
        row.names(nucl_len_mat) <- "Nucl"
    } else{
        nucl_len_mat <- NULL
    }
    if(all(is.null(c(bg_len_mat,tf_len_mat,nucl_len_mat)) && is.null(opt$ftpmodelyaml))){
        cli::cli_abort("At least one of the models (background noise, TF or nucleosome) must be defined. Please provide at least one of the following parameters: --bgnoiselengths, --tflengths, --nucllengths or --ftpmodelyaml.")
    }

    ftp_len_mat <- rbind(bg_len_mat,tf_len_mat,nucl_len_mat)

    ## check if minimum and maximum ftp lengths in the spectrum and in the ftp_len_mat are consistent
    ftplen_range <- range(fsa_data$ftp_spectrum$ftp_length)
    ftp_len_mat[which(ftp_len_mat[,1] < ftplen_range[1]),1] <- ftplen_range[1]
    ftp_len_mat[which(ftp_len_mat[,2] > ftplen_range[2]),2] <- ftplen_range[2]
    ftp_models <- nomeR::get_ftp_PWMs_for_prediction(ftp_spectrum = as.data.frame(fsa_data$ftp_spectrum),
                                                     ftp_len_mat = ftp_len_mat,
                                                     bg_cover = fsa_data$bgcoverprior,
                                                     ftp_protect_prob = fsa_data$ftpprotectprob)
    ftp_models <- list("ftp_models" = ftp_models,
                       "bgprotectprob" = fsa_data$bgprotectprob,
                       "bgcoverprior" = fsa_data$bgcoverprior)
    ## retrieve names of groups for which we will calculate enrichments etc.
    ftp_group_names <- unique(sapply(ftp_models$ftp_models,function(x){x$GROUP}))
}


#### DEFINE REGIONS AND GENOMIC BINNING FOR PREDICTION ####
##### construct seqinfo object from bam header ######
bam_header <- Rsamtools::scanBamHeader(files = opt$bamfile,what = c("targets"))
seqinfo_bam <- Seqinfo::Seqinfo(seqnames = names(bam_header[[1]]$targets),
                                bam_header[[1]]$targets)

if(!is.null(opt$regionbed)){
    ## if bed file is provided load the region
    sel_regions <- rtracklayer::import.bed(opt$regionbed)
} else{
    cli::cli_alert_info(text = "Creating genomic tiles of length {opt$chunksize}bp for genome-wide analysis.")
    ## if no bed file provided, load bam header and extract names and seqlengths
    sel_regions <- GRanges(seqinfo_bam)
}
## tile the regions into opt$chunksize chunks
sel_regions <- unlist(tile(sel_regions,width = opt$chunksize))
#names(sel_regions) <- as.character(sel_regions)
sel_regions <- GenomicRanges::sort(sel_regions,ignore.strand=TRUE)

## define chunk_id for parquet partitioning and other
sel_regions$chunk_id <- 1:length(sel_regions)
names(sel_regions) <- paste0("chunk_id_",sel_regions$chunk_id)
## export bed with genomic chunks
rtracklayer::export.bed(sel_regions,
                        con = file.path(opt$outputdir,
                                        "regions.bed"))


#### function for writing data.table into wig file ####

write_dt_wig <- function(dt,
                         tilewid,
                         scoreColumn,
                         file,
                         ...){
    ## delete file if exists
    unlink(file,force=TRUE)

    ## keep seqnames, start and scoreColumn only
    sel_cols = c("seqnames", "start", scoreColumn)
    dt2 <- dt[, ..sel_cols]
    ## remove data with NAs in scoreColumn

    dt2 <- dt2[!is.na(get(scoreColumn))]
    if(nrow(dt2) == 0){
        cli::cli_abort("All values in column {scoreColumn} are NAs.")
    }

    ## split dt2 by chr
    dt2 <- split(dt2,f=dt2[["seqnames"]])

    ## write each list item into file
    scol <- c("start", scoreColumn)
    fout <- lapply(seq_along(dt2),
                   function(i){

                       ## write a header
                       chrom <- dt2[[i]][["seqnames"]][1]
                       data.table::fwrite(
                           x = list(paste("variableStep",
                                          paste0("chrom=",chrom),
                                          paste0("span=",tilewid),
                                          sep=" ")),
                           file = file,
                           append=TRUE,
                           sep = "\t",
                           col.names = FALSE,
                           scipen = 9999
                       )
                       ## write data
                       data.table::fwrite(
                           dt2[[i]][, ..scol],
                           file = file,
                           append=TRUE,
                           sep = "\t",
                           col.names = FALSE,
                           scipen = 9999,
                           ...
                       )
                   })
}

#### function for convert character to safe filename ####
safe_filename <- function(x, replacement = "_", max_length = 255) {
    x <- as.character(x)
    # Remove reserved characters (Windows/Unix): < > : " / \ | ? * and control chars
    x <- gsub("[<>:\"/\\\\|?*[:cntrl:]]", replacement, x)
    # Collapse whitespace to replacement
    x <- gsub("\\s+", replacement, x)
    # Collapse repeated replacements
    x <- gsub(paste0("(", replacement, ")+"), replacement, x, fixed = FALSE)
    # Strip leading/trailing replacements, spaces, and dots (trailing dots/spaces illegal on Windows)
    x <- gsub(paste0("^[", replacement, ". ]+|[", replacement, ". ]+$"), "", x)
    # Guard against Windows reserved device names
    reserved <- "^(CON|PRN|AUX|NUL|COM[1-9]|LPT[1-9])$"
    x <- ifelse(grepl(reserved, x, ignore.case = TRUE), paste0(x, replacement), x)
    # Truncate
    x <- substr(x, 1, max_length)
    # Fallback for empty result
    x <- ifelse(nchar(x) == 0, "unnamed", x)
    x
}


### load data for correction of sequence bias ###
control_params <- NULL
if(opt$correctseqbias == "BetaCorrect"){
    control_params <- data.table::data.table(read.table(opt$betaparams, header = TRUE))

    ## check whether the control parameters are valid
    if(nrow(control_params) == 0){
        stop("No control parameters found for sequence bias correction.")
    }
    ## check whether the control parameters contain all required columns
    req_cols <- c("seqcont","alpha_pos", "beta_pos", "alpha_neg", "beta_neg")
    missing  <- setdiff(req_cols, names(control_params))
    if (length(missing))
        cli::cli_abort("File {.file {opt$betaparams}} is missing columns for method 'BetaCorrect': {paste(missing, collapse=', ')}")
    ## check if sequence contexts in control parameters are consistent with the k-mer width and reference sequence (if provided)
    seqconts <- grep("OTHER", control_params[["seqcont"]], invert = TRUE, value = TRUE)
    if(any(nchar(seqconts) != opt$kmer)){
        cli::cli_abort("Sequence contexts in {.file {opt$betaparams}} don't match the specified k-mer width of {opt$kmer}. Please check the input file and parameters.")
    }
} else if(opt$correctseqbias == "BetaUniform"){
    control_params <- data.table::data.table(read.table(opt$betaunifparams, header = TRUE))
    ## check whether the control parameters contain all required columns
    req_cols <- c("seqcont","alpha_pos", "beta_pos", "eps_pos", "alpha_neg", "beta_neg", "eps_neg")
    missing  <- setdiff(req_cols, names(control_params))
    if (length(missing))
        cli::cli_abort("File {.file {opt$betaunifparams}} is missing columns for method 'BetaUniform': {paste(missing, collapse=', ')}")
    ## check if sequence contexts in control parameters are consistent with the k-mer width and reference sequence (if provided)
    seqconts <- grep("OTHER", control_params[["seqcont"]], invert = TRUE, value = TRUE)
    if(any(nchar(seqconts) != opt$kmer)){
        cli::cli_abort("Sequence contexts in {.file {opt$betaunifparams}} don't match the specified k-mer width of {opt$kmer}. Please check the input file and parameters.")
    }


} else{
    if(!opt$filterkmerblacklist){
        opt$kmer <- 0L
        opt$refseq <- NULL
    }
}
## filter k-mers if a blacklist is provided
if(!is.null(opt$kmerblacklist)){
    ## load blacklist
    kmer_blacklist <- read.table(opt$kmerblacklist,header=F,col.names = c("kmer"))
} else{
    kmer_blacklist <- NULL
}




#### FOOTPRINT PREDICTION ####
options(cli.width = 500)

cli::cli_progress_step("Footprint prediction using mode {opt$ftpdecoding}. Output will be stored in {opt$outputdir}")
##### define temporary filenames and files for final output #####
## all temporary files for each genomic chunk will be stored into a separate folder
dir.create(opt$tempdir,recursive = TRUE,showWarnings = FALSE)

get_outfile_structure <- function(ftp_grnames,
                                  loc_folder = ".",
                                  noftps = FALSE,
                                  noenrich = FALSE){

    ### structure for the enrichments folder (in outputdir or tempfolder/chunk_id the structure is similar)
    # enrichments
    #     |
    #     |- frag_data_cover.bw (file with data coverage)
    #     |- mod_prob_mean.bw (file with average mod_prob)
    #     |
    #     |- Z-statistics
    #     |         |
    #     |         |- ftp1.bw
    #     |         |- ftp2.bw
    #     |         ....
    #     |
    #     |- mean_scores
    #     |         |- ftp1.bw
    #     |         |...
    #     |- neglPVal
    #     |         |- ftp1.bw
    #     |         |...
    #     |- neglFDR
    #     |         |- ftp1.bw
    #     |         |...
    # footprints
    #     |
    #     |- ftp1.bed
    #     |- ftp2.bed
    #     ...

    ## get BED files
    if(!noftps){
        ## folder BED files with footprints
        ftp_dir <- file.path(loc_folder,"footprints")

        ## for every ftp_group from ftp_grnames define to which file they should be stored
        ## create safe filenames. Special for background: convert name to linkers
        bedfnames <- safe_filename(paste0(gsub("background","linkers",ftp_grnames),".bed"))
        ## bed_outfile has relative path and the final path should be preceded either by opt$outputdir or opt$tempfolder/chunk_id
        ftp_bed_outfiles_df <- data.frame(ftp_group = ftp_grnames,
                                          filetype = "footprints",
                                          fileext = "bed",
                                          file = file.path(ftp_dir,bedfnames))
    } else{
        ftp_bed_outfiles_df <- NULL
    }

    if(!noenrich){
        enrich_dir <- file.path(loc_folder,"enrichments")
        wig_raw_data_files <- data.frame(ftp_group = NA,
                                         filetype = c("frag_data_cover",
                                                      "mod_prob_mean"),
                                         fileext = "wig",
                                         scoreColumn = c("frag_cover","mod_prob_mean"),
                                         file = file.path(enrich_dir,c("frag_data_cover.wig",
                                                                       "mod_prob_mean.wig")))

        ## create data.frame describing where wig should be stored
        ## mean scores
        scores_wig_files <- data.frame(ftp_group = ftp_grnames,
                                       filetype = "mean_score",
                                       fileext = "wig",
                                       scoreColumn = paste0(ftp_grnames,"_score_mean"),
                                       file = file.path(enrich_dir,"mean_scores",safe_filename(paste0(ftp_grnames,"_mean_score.wig"))))
        ## z-statistics
        zstats_wig_files <- data.frame(ftp_group = ftp_grnames,
                                       filetype = "Z-statistics",
                                       fileext = "wig",
                                       scoreColumn = paste0(ftp_grnames,"_score_mean_Zstat"),
                                       file = file.path(enrich_dir,"Z-statistics",safe_filename(paste0(ftp_grnames,"_Zstat.wig"))))

        ## negative log10 p.val
        neglPval_wig_files <- data.frame(ftp_group = ftp_grnames,
                                         filetype = "neglPval",
                                         fileext = "wig",
                                         scoreColumn = paste0(ftp_grnames,"_neglogPval"),
                                         file = file.path(enrich_dir,"neglPval",safe_filename(paste0(ftp_grnames,"_neglPval.wig"))))

        ## negative log10 FDR
        neglFDR_wig_files <- data.frame(ftp_group = ftp_grnames,
                                        filetype = "neglFDR",
                                        fileext = "wig",
                                        scoreColumn = paste0(ftp_grnames,"_neglogFDR"),
                                        file = file.path(enrich_dir,"neglFDR",safe_filename(paste0(ftp_grnames,"_neglFDR.wig"))))
        wig_outfiles_df <- bind_rows(wig_raw_data_files,
                                     scores_wig_files,
                                     zstats_wig_files,
                                     neglPval_wig_files,
                                     neglFDR_wig_files)
    } else{
        wig_outfiles_df <- NULL
    }

    file_struct <- bind_rows(ftp_bed_outfiles_df,
                             wig_outfiles_df)
    return(file_struct)

}


### define threads share between mclapply and predict_footprints_SE
#mclapply_threads_share <- 0.5
#threads_mclap <- max(floor(opt$threads * mclapply_threads_share),1)
#threads_predict <- max(floor(opt$threads/threads_mclap),1)
threads_mclap <- opt$threads
threads_predict <- 1


format_index <- function(i, max_i) {
    sprintf(paste0("%0", nchar(max_i), "d"), i)
}

##### LOOP FOR PREDICTING FOOTPRINTS IN GENOMIC CHUNKS ######
n_chunks <- length(sel_regions)

## by default assayName is mod_prob. but if sequence bias correction was done it should be changed to mod_prob_corrected
assayName <- "mod_prob"
if(opt$correctseqbias != "no_correction" && !is.null(control_params))
    assayName <- "mod_prob_corrected"
if(!is.null(kmer_blacklist))
    kmer_blacklist <- Biostrings::DNAStringSet(kmer_blacklist[["kmer"]])

files2concat <- bind_rows(mcprogress::pmclapply(
    mc.cores = threads_mclap,
    #lapply(
    X = seq_along(sel_regions),
    FUN = function(i){
        reg <- sel_regions[i]
        # ---- 1. create a folder for storing temporary files ----
        #chunk_id <- sprintf("%05d", i)
        chunk_id <- format_index(i,max_i = n_chunks)
        chunk_dir <- file.path(opt$tempdir, paste0("chunk_", chunk_id))

        # Safe: only this worker touches this directory
        dir.create(chunk_dir, showWarnings = FALSE, recursive = TRUE)

        ## get file structure for the current chunk
        chunk_file_struct <- get_outfile_structure(ftp_grnames = ftp_group_names,
                                                   loc_folder = chunk_dir,
                                                   noftps = opt$noftps,
                                                   noenrich = opt$noenrich)

        # ---- 2. run  existing predictions -----
        ### 05.03.2026: do not save the fst files. they take way too much space for genome-wide analysis

        ## load data
        #cli::cli_progress_step("Fetching data from BAM file for region {as.character(reg)} [{i} out of {n_chunks}]")
        se <- SingleMoleculeGenomicsIO::readModBam(bamfiles = opt$bamfile,
                                                   regions = reg,
                                                   modbase = "a",
                                                   sequenceContextWidth = opt$kmer,
                                                   sequenceReference = opt$refseq,
                                                   trim = TRUE,
                                                   BPPARAM = BiocParallel::MulticoreParam(workers = threads_predict),
                                                   verbose = F)

        if(se$n_reads == 0){
            return(NULL)
        }

        rownames(se) <- seq_along(se)

        # ----- 2.1 (optional) correction of sequence bias ------
        if(opt$correctseqbias != "no_correction" && !is.null(control_params)){
            cli::cli_inform("Correction of sequence biases ({opt$correctseqbias})")
            se <- nomeR::correct_modprob_SE(se,
                                            control_params = control_params,
                                            method         = opt$correctseqbias,
                                            isotonic       = !opt$noisotonic,
                                            qnorm_to_raw   = opt$quantnorm)
        }
        if(opt$filterkmerblacklist && !is.null(kmer_blacklist)){
            cli::cli_inform("Filtering blacklisted k-mer sequence contexts")
            keep_rows <- which(!(rowData(se)[,"sequenceContext"] %in% kmer_blacklist))
            se <- se[keep_rows, ]
        }
        #cli::cli_progress_step("Running nomeR::predict_footprints_SE for {se$n_reads} fragments in the region {as.character(reg)} [{i} out of {n_chunks}]")
        cli::cli_inform("Running prediction on assay {assayName} for {se$n_reads} fragments in the region {as.character(reg)} [{i} out of {n_chunks}]")
        pred_list_dt <- nomeR::predict_footprints_SE(se = se,
                                                     assayName = assayName,
                                                     footprint_models = ftp_models$ftp_models,
                                                     bgprotectprob = ftp_models$bgprotectprob,
                                                     bgcoverprior = ftp_models$bgcoverprior,
                                                     ftpConfigMethod = opt$ftpdecoding,
                                                     aggrByGroup = TRUE,
                                                     returnAs = "data.table",
                                                     ncpu = threads_predict,
                                                     #ncpu = opt$threads,
                                                     min_frag_data_len = 0L,
                                                     min_frag_data_dens = 0.0,
                                                     verbose = F)

        ## remove sample prefix from molecule names
        pred_list_dt$FOOTPRINT_CONF <- pred_list_dt$FOOTPRINT_CONF[,fragID := str_remove(fragID,"^s1-")]
        pred_list_dt$COVER_PROB <- pred_list_dt$COVER_PROB[,fragID := str_remove(fragID,"^s1-")]

        # ---- 4. write footprints into temporary file for the current genomic chunk ----
        if(!opt$noftps){
            cli::cli_inform("Writing predicted footprints into temporary BED files for region {as.character(reg)} [{i} out of {n_chunks}]")

            ## temporary files
            ftp_bed_files <- chunk_file_struct %>%
                filter(filetype == "footprints")

            for(i in 1:nrow(ftp_bed_files)){
                dir.create(dirname(ftp_bed_files$file[i]),
                           recursive = TRUE,showWarnings = FALSE)
                nomeR::write_ftp_decoding_bed(ftp_decode_dt = pred_list_dt$FOOTPRINT_CONF[ftp_group == ftp_bed_files$ftp_group[i]],
                                              file = ftp_bed_files$file[i],
                                              append=FALSE)
            }
        }

        # ---- 5. calculate and write enrichments for tiles ----
        if(!opt$noenrich){
            cli::cli_inform("Calculating enrichments and writing BIGWIG files for region {as.character(reg)} [{i} out of {n_chunks}]")

            enr_dt <- nomeR::calculate_tile_footprint_enrichments(cover_dt = pred_list_dt$COVER_PROB,
                                                                  tile_width = tile_width,
                                                                  tile_step = tile_step,
                                                                  ftp_names = ftp_group_names,
                                                                  threads = threads_predict,
                                                                  psc = opt$coverpsc
            )

            ## adjust end so that width is tile_step to remove tile overlapping
            enr_dt <- enr_dt[,end := start + tile_step - 1]
            data.table::setorder(enr_dt, seqnames, start)

            ## add frag coverage, -log10(p-val) and -log10(FDR)
            pfdr_cols <- c(paste0(ftp_group_names,"_score_mean_Ztest_pval"),
                           paste0(ftp_group_names,"_score_mean_FDR"))
            neglognames <-c(paste0(ftp_group_names,"_neglogPval"),
                            paste0(ftp_group_names,"_neglogFDR"))
            pval_psc <- 1e-100
            enr_dt <- enr_dt[, (c("frag_cover", neglognames)) := {
                res <- lapply(.SD, function(x) { -log10(x + pval_psc) })
                c(list(n_data_points / n_inf_pos), res)
            }, .SDcols = pfdr_cols]

            ## write wig files
            cur_wig_files <- chunk_file_struct %>%
                filter(fileext == "wig")
            for(i in 1:nrow(cur_wig_files)){
                dir.create(dirname(cur_wig_files$file[i]),recursive = TRUE,showWarnings = FALSE)
                write_dt_wig(dt = enr_dt,
                             tilewid = tile_step,
                             scoreColumn = cur_wig_files$scoreColumn[i],
                             file = cur_wig_files$file[i])


            }
        }

        ### define data.frame for concatenating wig files into one wig file
        curwigs2cat <- chunk_file_struct %>% mutate(chunk_id = chunk_id,scoreColumn=NULL)
        return(curwigs2cat)


    }
))


### create data.frame for concatenating bed files.
merged_bed_files_df <- get_outfile_structure(ftp_grnames = ftp_group_names,
                                             loc_folder = opt$outputdir,
                                             noftps = opt$noftps,
                                             noenrich = TRUE) %>%
    select(ftp_group,filetype,fileext,file)
### create data.frame for concatenating wig files. The reason why we store them in tempdir is bacause we need to convert them into bw files later and place in outputdir
merged_wig_files_df <- get_outfile_structure(ftp_grnames = ftp_group_names,
                                             loc_folder = opt$tempdir,
                                             noftps = TRUE,
                                             noenrich = opt$noenrich)%>%
    select(ftp_group,filetype,fileext,file)

merged_files_df <- bind_rows(merged_bed_files_df,merged_wig_files_df) %>%
    rename(file="merged_file")

## check if any files are missing
missing_files <- files2concat$file[!file.exists(files2concat$file)]
if(length(missing_files) > 0){
    missing_files <- paste(missing_files,collapse = "\n")
    cli::cli_warn("Some files are missing:\n {missing_files}")
}

files2concat <- files2concat %>%
    filter(file.exists(file)) %>%
    left_join(merged_files_df, by = join_by(ftp_group,filetype,fileext))

#### functions for concatenating temporary files ####
concat_files <- function(files,
                         output,
                         buffer_size = 8 * 1024 * 1024) {
    stopifnot(length(files) > 0)

    out <- file(output, open = "wb")
    for (f in files) {
        inp <- file(f, open = "rb")
        repeat {
            buf <- readBin(inp, what = "raw", n = buffer_size)
            if (!length(buf)) break
            writeBin(buf, out)
        }
        close(inp)
    }
    close(out)
    invisible(output)
}

cli::cli_progress_step("Merging temporary files")

merged_out_files_df <- do.call(rbind,(by(data = files2concat, INDICES = files2concat$merged_file,
                                         FUN = function(df){

                                             mrg_file <- unique(df$merged_file)
                                             if(length(mrg_file) != 1){
                                                 cli::cli_warn("More than one file to store merged data. Using only the first one.")
                                             }
                                             dir.create(dirname(mrg_file),recursive = TRUE,showWarnings = F)
                                             mrg_file <- concat_files(files = df$file,
                                                                      output = mrg_file[1])

                                             outd <- df[1,] %>% select(ftp_group,filetype,fileext) %>%
                                                 mutate(merged_file = mrg_file)

                                             return(outd)

                                         },simplify=FALSE)))
row.names(merged_out_files_df) <- NULL

## convert concatenated wig files into bigwig files
merged_wg_files <- merged_out_files_df %>%
    filter(fileext == "wig")
merged_bw_files_df <- merged_wg_files %>%
    left_join(get_outfile_structure(ftp_grnames = ftp_group_names,
                                    loc_folder = opt$outputdir,
                                    noftps = TRUE,
                                    noenrich = opt$noenrich) %>%
                  mutate(fileext = "bw",
                         file = gsub("\\.wig$",".bw",file)) %>%
                  rename(file="bw_file"),
              by = c("ftp_group","filetype"))
cli::cli_progress_step("Converting Wig to BigWig")
if(!opt$noenrich){
    fo <- lapply(1:nrow(merged_bw_files_df),
                 function(i){
                     if(!is.na(merged_bw_files_df$bw_file[i]) && file.exists(merged_bw_files_df$merged_file[i])){
                         dir.create(dirname(merged_bw_files_df$bw_file[i]),recursive = TRUE,showWarnings = F)
                         rtracklayer::wigToBigWig(x = merged_bw_files_df$merged_file[i],
                                                  seqinfo = seqinfo_bam,
                                                  dest = merged_bw_files_df$bw_file[i])
                     }
                 })
}
cli::cli_progress_done()

## final set of create files
final_file_set <- bind_rows(merged_out_files_df %>%
                                filter(fileext == "bed") %>%
                                rename(merged_file = "file"),
                            merged_bw_files_df %>% select(ftp_group,filetype,fileext.y,bw_file) %>%
                                rename(fileext.y = 'fileext',bw_file = "file"))


cli::cli_alert_success("Footprint prediction has been finished!")
cli::cli_h1("The following output files have been stored:")
cli::cli_ol()
final_file_set <- final_file_set %>%
    filter(file.exists(file))
if(!opt$noftps){
    cli::cli_li("{.field BED12 files with footprints are located in the folder {.file footprints}}")
}
if(!opt$noenrich){
    cli::cli_li("{.field Genomic BigWig files in sliding windows of width {tile_width} with step {tile_step}:}")
    ulid <- cli::cli_ul()
    cov_fls <- final_file_set %>% filter(filetype == "frag_data_cover")
    if(file.exists(cov_fls$file))
        cli::cli_li("Data coverage: {.file enrichments/frag_data_cover.bw}")
    modpr_fls <- final_file_set %>% filter(filetype == "mod_prob_mean")
    if(file.exists(cov_fls$file))
        cli::cli_li("Average modification probabilities: {.file enrichments/mod_prob_mean.bw}")

    cli::cli_li("Average predicted log10-transformed posterior coverage probabilities for each footprint model are located in the folder: {.file enrichments/mean_scores}")
    cli::cli_li("Z-statistics for enrichment of posterior coverage probabilities for each footprint model are located in the folder: {.file enrichments/Z-statistics}")
    cli::cli_li("Corresponding -log10(p-values) are located in the folder: {.file enrichments/neglPval}")
    cli::cli_li("Similarly -log10(FDR) for enrichments are located in the folder: {.file enrichments/neglFDR}")
    cli::cli_end(ulid)
}




cli::cli_end()
closeAllConnections()
## delete temprary folder
unlink(x = opt$tempdir,recursive = TRUE)
