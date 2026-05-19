#!/usr/bin/env Rscript
## R script (v6) for predicting footprints in SAMOSA/FiberSeq generated BAM file

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
    make_option(c("-s", "--fsayaml"),
                type="character",
                help="[OPTIONAL] Path to an YAML file containing results of footprint spectrum analysis.
                If --fsayaml --ftpmodelyaml are not provided, the footprint spectral analysis will be done for input BAM file.
                If --ftpmodelyaml is provided the --fsayaml will be ignored."),
    make_option(c("-f", "--ftpmodelyaml"),
                type="character",
                help="[OPTIONAL] Path to an YAML file containing footprint models.
                If not provided, the footprint models will be constructed using footprint spectrum in --fsayaml or
                based on footprint spectral analysis of the input BAM file."),
    make_option(c("-r", "--regionbed"),
                type="character",
                help="[OPTIONAL] BED file with region of interest. If not provided the analysis will be performed genome-wide."),

    ### options specifying model
    make_option(c("-l", "--ftpmodeltype"),
                type="character",
                default="fast",
                help="Type of footprint models to use for prediction {fast,medium,slow} [default %default].
                If --ftpmodelyaml is set the --ftpmodeltype has no effect."),
    make_option(c("-d", "--ftpdecoding"),
                type="character",
                default="PosteriorDecoding",
                help="Method for footprint decoding. Can be 'PV' - Posterior Viterbi; 'PosteriorDecoding' or 'Viterbi' [default %default]"),
    make_option(c("-m", "--thresholdmod"),
                type="numeric",
                default = 0.5,
                help="Threshold for modification probability to binarize the data into accessible and protected positions.
                Probabilities equal or higher will be considered accessible, otherwise protected. [default %default]"),
    make_option(c("-a", "--tilewidthstep"),
                type="character",
                default = "500,250",
                help="Tile width and step (format: \"width,step\") for aggregating and calculating p-values for TF and BG. [default %default]"),

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
    ### output options
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

    ### running options
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

    ### library options
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

##### check input model options #####
if(!opt$ftpmodeltype %in% c("slow","medium","fast")){
    cli::cli_abort("--ftpmodeltype must be one of the following: \"slow\",\"medium\",\"fast\"")
}

if(!opt$ftpdecoding %in% c("PV","PosteriorDecoding","Viterbi")){
    cli::cli_abort("ftpdecoding must be one of the following: \"PV\",\"PosteriorDecoding\",\"Viterbi\"")
}

if(opt$thresholdmod <= 0 || opt$thresholdmod >= 1){
    cli::cli_abort("--thresholdmod must be within (0,1). Provided --thresholdmod {.val opt$thresholdmod}")
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



#### DEFINE FOOTPRINT MODELS ####
if(opt$ftpmodeltype == "slow"){
    ## create footprint models
    ftp_len_mat <- rbind(c(2,5,1),
                         c(20,50,1),
                         c(100,150,1))

} else if(opt$ftpmodeltype == "medium"){
    ftp_len_mat <- rbind(c(2,5,1),
                         c(20,50,5),
                         c(100,150,5))
} else{
    ftp_len_mat <- rbind(c(2,5,1),
                         c(20,50,10),
                         c(100,150,10))

}

colnames(ftp_len_mat) <- c("min_ftp_len","max_ftp_len","by")
row.names(ftp_len_mat) <- c("background","TF","Nucl")


#### LOAD FOOTPRINT MODELS OR RESULTS OF FOOTPRINT SPECTRUM ANALYSIS ####
##### construct seqinfo object from bam header ######
bam_header <- Rsamtools::scanBamHeader(files = opt$bamfile,what = c("targets"))
seqinfo_bam <- Seqinfo::Seqinfo(seqnames = names(bam_header[[1]]$targets),
                                bam_header[[1]]$targets)


if(!is.null(opt$ftpmodelyaml)){
    ## load footprint model from a YAML file
    ftp_models <- yaml::read_yaml(opt$ftpmodelyaml)
    ## check if all fields are present
    if(is.null(ftp_models$ftp_models) ||
       is.null(ftp_models$bgprotectprob) ||
       is.null(ftp_models$bgcoverprior)){
        cli::cli_abort("Can't find one of \"ftp_models\" \"bgprotectprob\" \"bgcoverprior\". All three must be present in the input YAML: {opt$ftpmodelyaml}")
    }


} else if(!is.null(opt$fsayaml)){
    ## load footprint spectrum from a YAML file
    fsa_data <- yaml::read_yaml(opt$fsayaml)
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
}


#### DEFINE REGIONS AND GENOMIC BINNING FOR PREDICTION ####
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

#### PREDICT FOOTPRINTS ####
options(cli.width = 500)

cli::cli_progress_step("Footprint prediction using mode {opt$ftpmodeltype}. Output will be stored in {opt$outputdir}")
##### define temporary filenames and files for final output #####
## all temporary files for each genomic chunk will be stored into a separate folder
dir.create(opt$tempdir,recursive = TRUE,showWarnings = FALSE)


if(!opt$noftps){
    ## folder with temporary files
    ftp_dir <- file.path(opt$outputdir,"footprints")
    ## delete existing directory and create a new one
    unlink(ftp_dir,recursive=TRUE)
    dir.create(ftp_dir)

    nucl_bed <- file.path(ftp_dir,"nucleosomes_ftp.bed")
    tf_bed <- file.path(ftp_dir,"tf_ftp.bed")
    linker_bed <- file.path(ftp_dir,"linkers_ftp.bed")

}
if(!opt$noenrich){

    bw_dir <- file.path(opt$outputdir,"enrichments")
    ## delete existing directory and create a new one
    unlink(bw_dir,recursive=TRUE)
    dir.create(bw_dir)

    ## nDataPoints - bw file with tile coverage
    n_frag_cover_wig <- file.path(opt$tempdir,"frag_data_cover.wig")
    n_frag_cover_bw <- file.path(bw_dir,"frag_data_cover.bw")

    ## files for BG
    #### average BGscore - mean of BGscores across informative positions overlapping each tile
    bg_score_wig <- file.path(opt$tempdir,"bg_score_mean.wig")
    bg_score_bw <- file.path(bw_dir,"bg_score_mean.bw")
    #### BG score Z-statistics
    bg_score_Zstat_wig <- file.path(opt$tempdir,"bg_score_mean_Zstat.wig")
    bg_score_Zstat_bw <- file.path(bw_dir,"bg_score_mean_Zstat.bw")
    #### -log10(p-vals)
    bg_pval_wig <- file.path(opt$tempdir,"bg_neglogPval.wig")
    bg_pval_bw <- file.path(bw_dir,"bg_neglogPval.bw")
    #### -log10(FDR)
    bg_fdr_wig <- file.path(opt$tempdir,"bg_neglogFDR.wig")
    bg_fdr_bw <- file.path(bw_dir,"bg_neglogFDR.bw")


    ## files for TF
    #### average TFscore - mean of TFscores across informative positions overlapping each tile
    tf_score_wig <- file.path(opt$tempdir,"tf_score_mean.wig")
    tf_score_bw <- file.path(bw_dir,"tf_score_mean.bw")
    #### TF score Z-statistics
    tf_score_Zstat_wig <- file.path(opt$tempdir,"tf_score_mean_Zstat.wig")
    tf_score_Zstat_bw <- file.path(bw_dir,"tf_score_mean_Zstat.bw")
    #### -log10(p-vals)
    tf_pval_wig <- file.path(opt$tempdir,"tf_neglogPval.wig")
    tf_pval_bw <- file.path(bw_dir,"tf_neglogPval.bw")
    #### -log10(FDR)
    tf_fdr_wig <- file.path(opt$tempdir,"tf_neglogFDR.wig")
    tf_fdr_bw <- file.path(bw_dir,"tf_neglogFDR.bw")

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


### load data for correction of sequence bias ###
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
        cli::cli_abort("For sequence bias correction both --negbetas and --posbetas files must be provided. Please provide both files or set --correctseqbias to 'no_correction'.")
    }
} else{
    opt$kmer <- 0L
    opt$refseq <- NULL
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



##### LOOP FOR PREDICTING FOOTPRINTS IN GENOMIC CHUNKS ######
n_chunks <- length(sel_regions)

## by default assayName is mod_prob. but if sequence bias correction was done it should be changed to mod_prob_corrected
assayName <- "mod_prob"
if(opt$correctseqbias != "no_correction" && !is.null(negcontrol_shapes) && !is.null(poscontrol_shapes))
    assayName <- "mod_prob_corrected"
if(!is.null(kmer_blacklist))
    kmer_blacklist <- Biostrings::DNAStringSet(kmer_blacklist[["kmer"]])

pred_out <- mcprogress::pmclapply(
    mc.cores = threads_mclap,
    X = seq_along(sel_regions),
    FUN = function(i){

        reg <- sel_regions[i]
        # ---- 1. create a folder for storing temporary files ----
        #chunk_id <- sprintf("%05d", i)
        chunk_id <- format_index(i,max_i = n_chunks)
        chunk_dir <- file.path(opt$tempdir, paste0("chunk_", chunk_id))

        # Safe: only this worker touches this directory
        dir.create(chunk_dir, showWarnings = FALSE, recursive = TRUE)

        # ---- 2. run  existing predictions -----
        ### 05.03.2026: do not save the fst files. they take way too much space for genome-wide analysis
        #if(!file.exists(cov_fst) || !file.exists(ftp_fst)){

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
        if(opt$correctseqbias != "no_correction"){
            if(!is.null(negcontrol_shapes) && !is.null(poscontrol_shapes)){
                cli::cli_inform("Bayesian beta correction of sequence biases")
                se <- nomeR::correct_modprob_SE(se,
                                         neg_control_shapes = negcontrol_shapes,
                                         pos_control_shapes = poscontrol_shapes,
                                         qnorm_to_raw = TRUE)

            }
            if(!is.null(kmer_blacklist)){
                cli::cli_inform("Filtering blacklisted kmers")                                
                ## remove k-mers
                keep_rows <- which(!(rowData(se)[,"sequenceContext"] %in% kmer_blacklist))
                se <- se[keep_rows, ]
            }
        }
        #cli::cli_progress_step("Running nomeR::predict_footprints_SE for {se$n_reads} fragments in the region {as.character(reg)} [{i} out of {n_chunks}]")
        cli::cli_inform("Running prediction on assay {assayName} for {se$n_reads} fragments in the region {as.character(reg)} [{i} out of {n_chunks}]")
        pred_list_dt <- nomeR::predict_footprints_SE(se = se,
                                                     assayName = assayName,
                                                     threshMod = opt$thresholdmod,
                                                     threshUnmod = opt$thresholdmod,
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
            nucl_temp_bed <- file.path(chunk_dir,"nucl_ftps.bed")
            tf_temp_bed <- file.path(chunk_dir,"tf_ftps.bed")
            linkers_temp_bed <- file.path(chunk_dir,"linkers_ftps.bed")

            ## write nucleosome ftp BED12
            nomeR::write_ftp_decoding_bed(ftp_decode_dt = pred_list_dt$FOOTPRINT_CONF[ftp_group == "Nucl"],
                                          file = nucl_temp_bed,
                                          append=FALSE)
            ## write TF ftp BED12
            nomeR::write_ftp_decoding_bed(ftp_decode_dt = pred_list_dt$FOOTPRINT_CONF[ftp_group == "TF"],
                                          file = tf_temp_bed,
                                          append=FALSE)
            ## reduce background footprints to get linkers
            linkers <- pred_list_dt$FOOTPRINT_CONF[ftp_group == "background"]
            linkers <- GenomicRanges::GRanges(seqnames = linkers$seqnames,
                               IRanges(start = linkers$start,
                                       width=linkers$width),
                               strand = linkers$strand,
                               fragID = linkers$fragID)
            linkers <- split(linkers,linkers$fragID)
            linkers <- GenomicRanges::reduce(linkers)
            linkers <- unlist(linkers)
            linkers$fragID <- names(linkers)
            linkers <- data.table::as.data.table(linkers)
            ## write linkers BED12
            nomeR::write_ftp_decoding_bed(ftp_decode_dt = linkers,
                                          file = linkers_temp_bed,
                                          append=FALSE)
        }

        # ---- 5. calculate and write enrichments for tiles ----
        if(!opt$noenrich){
            cli::cli_inform("Calculating enrichments and writing BIGWIG files for region {as.character(reg)} [{i} out of {n_chunks}]")

            ## temporary files
            n_frag_cover_temp_wig <- file.path(chunk_dir,"frag_data_cover.wig")
            ### BG
            bg_score_temp_wig <- file.path(chunk_dir,"bg_score_mean.wig")
            bg_score_Zstat_temp_wig <- file.path(chunk_dir,"bg_score_mean_Zstat.wig")
            bg_pval_temp_wig <- file.path(chunk_dir,"bg_pval.wig")
            bg_fdr_temp_wig <- file.path(chunk_dir,"bg_fdr.wig")

            ### TF
            tf_score_temp_wig <- file.path(chunk_dir,"tf_score_mean.wig")
            tf_score_Zstat_temp_wig <- file.path(chunk_dir,"tf_score_mean_Zstat.wig")
            tf_pval_temp_wig <- file.path(chunk_dir,"tf_pval.wig")
            tf_fdr_temp_wig <- file.path(chunk_dir,"tf_fdr.wig")            
            enr_dt <- nomeR::calculate_tile_BG_TF_enrichments(cover_dt = pred_list_dt$COVER_PROB,
                                                        tile_width = tile_width,
                                                        tile_step = tile_step,
                                                        threads = threads_predict
            )


            ## adjust end so that width is tile_step to remove tile overlapping
            enr_dt <- enr_dt[,end := start + tile_step - 1]
            data.table::setorder(enr_dt, seqnames, start)

            ## add frag coverage, -log10(p-val) and -log10(FDR)
            pval_psc <- 1e-100
            enr_dt <- enr_dt[,':='(frag_cover = n_data_points/n_inf_pos,
                                   bg_neglogPval = -log10(bg_score_mean_Ztest_pval + pval_psc),
                                   bg_neglogFDR  = -log10(bg_score_mean_FDR        + pval_psc),
                                   tf_neglogPval = -log10(tf_score_mean_Ztest_pval + pval_psc),
                                   tf_neglogFDR  = -log10(tf_score_mean_FDR        + pval_psc))]

            ## write mean_frag_cover - average fragment coverage per informative position
            write_dt_wig(dt = enr_dt,
                         tilewid = tile_step,
                         scoreColumn = "frag_cover",
                         file = n_frag_cover_temp_wig)

            #write bg wig files
            write_dt_wig(dt = enr_dt,
                         tilewid = tile_step,
                         scoreColumn = "bg_score_mean",
                         file = bg_score_temp_wig)
            write_dt_wig(dt = enr_dt,
                         tilewid = tile_step,
                         scoreColumn = "bg_score_mean_Zstat",
                         file = bg_score_Zstat_temp_wig)
            write_dt_wig(dt = enr_dt,
                         tilewid = tile_step,
                         scoreColumn = "bg_neglogPval",
                         file = bg_pval_temp_wig)
            write_dt_wig(dt = enr_dt,
                         tilewid = tile_step,
                         scoreColumn = "bg_neglogFDR",
                         file = bg_fdr_temp_wig)

            #write tf wig files
            write_dt_wig(dt = enr_dt,
                         tilewid = tile_step,
                         scoreColumn = "tf_score_mean",
                         file = tf_score_temp_wig)
            write_dt_wig(dt = enr_dt,
                         tilewid = tile_step,
                         scoreColumn = "tf_score_mean_Zstat",
                         file = tf_score_Zstat_temp_wig)
            write_dt_wig(dt = enr_dt,
                         tilewid = tile_step,
                         scoreColumn = "tf_neglogPval",
                         file = tf_pval_temp_wig)
            write_dt_wig(dt = enr_dt,
                         tilewid = tile_step,
                         scoreColumn = "tf_neglogFDR",
                         file = tf_fdr_temp_wig)


        }
    }

)



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

merge_temp_files <- function(filepattern,
                             tempdir,
                             chunk_id_vec,
                             outputfile){
    ## construct temprorary files names for merging
    #chunk_id <- sprintf("%05d", chunk_id_vec)
    chunk_id <- format_index(chunk_id_vec,
                             max_i = max(chunk_id_vec))
    chunk_dirs <- file.path(tempdir, paste0("chunk_", chunk_id))
    files <- file.path(chunk_dirs,
                       filepattern)
    ## filter files that do not exist
    files <- files[file.exists(files)]
    if(length(files) == 0)
        return(NULL)
    invisible(concat_files(files = files,
                           output = outputfile))
}


cli::cli_progress_step("Merging temporary files")
merge_files_df <- data.frame(fpat = character(0), fout = character(0), fbw = character(0))

if (!opt$noenrich) {
    merge_files_df <- rbind(merge_files_df, data.frame(
        fpat = c("frag_data_cover.wig",
                 "bg_score_mean.wig",
                 "bg_score_mean_Zstat.wig",
                 "bg_pval.wig",
                 "bg_fdr.wig",
                 "tf_score_mean.wig",
                 "tf_score_mean_Zstat.wig",
                 "tf_pval.wig",
                 "tf_fdr.wig"),
        fout = c(n_frag_cover_wig,
                 bg_score_wig,
                 bg_score_Zstat_wig,
                 bg_pval_wig,
                 bg_fdr_wig,
                 tf_score_wig,
                 tf_score_Zstat_wig,
                 tf_pval_wig,
                 tf_fdr_wig),
        fbw  = c(n_frag_cover_bw,
                 bg_score_bw,
                 bg_score_Zstat_bw,
                 bg_pval_bw,
                 bg_fdr_bw,
                 tf_score_bw,
                 tf_score_Zstat_bw,
                 tf_pval_bw,
                 tf_fdr_bw)
    ))
}

if (!opt$noftps) {
    merge_files_df <- rbind(merge_files_df, data.frame(
        fpat = c("nucl_ftps.bed", "tf_ftps.bed", "linkers_ftps.bed"),
        fout = c(nucl_bed, tf_bed, linker_bed),
        fbw  = c(NA_character_, NA_character_, NA_character_)
    ))
}

fo <- lapply(1:nrow(merge_files_df),
             function(i){
                 merge_temp_files(filepattern = merge_files_df$fpat[i],
                                  tempdir = opt$tempdir,
                                  chunk_id_vec = seq_along(sel_regions),
                                  outputfile = merge_files_df$fout[i])
             })


### convert enrichments wig to bigwigs
cli::cli_progress_step("Converting Wig to BigWig")
if(!opt$noenrich){
    fo <- lapply(1:nrow(merge_files_df),
                 function(i){
                     if(!is.na(merge_files_df$fbw[i]) && file.exists(merge_files_df$fout[i])){
                         rtracklayer::wigToBigWig(x = merge_files_df$fout[i],
                                                  seqinfo = seqinfo_bam,
                                                  dest = merge_files_df$fbw[i])
                     }
                 })
}
cli::cli_progress_done()


cli::cli_alert_success("Footprint prediction has been finished!")
cli::cli_h1("The following output files have been stored:")
cli::cli_ol()

if(!opt$noftps){
    cli::cli_li("{.field BED12 files with footprints:}")
    ulid <- cli::cli_ul()
    if(file.exists(nucl_bed))
        cli::cli_li("Nucleosomes: {.file {nucl_bed}}")
    if(file.exists(tf_bed))
        cli::cli_li("Transcription factors: {.file {tf_bed}}")
    if(file.exists(linker_bed))
        cli::cli_li("Linkers: {.file {linker_bed}}")
    cli::cli_end(ulid)
}
if(!opt$noenrich){
    cli::cli_li("{.field Average BG/TF scores and neg. log-p-values BigWig files in sliding windows of width {tile_width} with step {tile_step}:}")
    ulid <- cli::cli_ul()
    if(file.exists(n_frag_cover_bw))
        cli::cli_li("Data coverage: {.file {n_frag_cover_bw}}")

    cli::cli_li("Accessible positions (background):")
    bglid <- cli::cli_ul()
    if(file.exists(bg_score_bw))
        cli::cli_li("Average BG-scores: {.file {bg_score_bw}}")
    if(file.exists(bg_score_Zstat_bw))
        cli::cli_li("Z-statistics for BG-scores: {.file {bg_score_Zstat_bw}}")
    if(file.exists(bg_pval_bw))
        cli::cli_li("-log10(p-values): {.file {bg_pval_bw}}")
    if(file.exists(bg_fdr_bw))
        cli::cli_li("-log10(FDR): {.file {bg_fdr_bw}}")

    cli::cli_end(bglid)

    cli::cli_li("TF-covered positions:")
    tflid <- cli::cli_ul()
    if(file.exists(tf_score_bw))
        cli::cli_li("Average TF-scores: {.file {tf_score_bw}}")
    if(file.exists(tf_score_Zstat_bw))
        cli::cli_li("Z-statistics for TF-scores: {.file {tf_score_Zstat_bw}}")
    if(file.exists(tf_pval_bw))
        cli::cli_li("-log10(p-values): {.file {tf_pval_bw}}")
    if(file.exists(tf_fdr_bw))
        cli::cli_li("-log10(FDR): {.file {tf_fdr_bw}}")

    cli::cli_end(tflid)

    cli::cli_end(ulid)
}

cli::cli_end()
closeAllConnections()
