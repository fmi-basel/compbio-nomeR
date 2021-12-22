


#' Fetch statistics of duplicated fragments
#'
#' @param bamfiles \code{character} paths to bam files
#' @param samplenames \code{character} names of samples. If \code{collapseBySample} is \code{TRUE} counts are aggregated across bam files with same sample name.
#' @param regions \linkS4class{GRanges} object for regions of interest.
#' @param genome \code{character} path to fasta file which was used as reference for alignment.
#' @param collapseBySample \code{logical} indicating whether to collapse counts for bam files with same \code{samplenames}. If \code{FALSE} prefix \code{s} followed by index is added to \code{samplenames}.
#' @param mapqMin \code{integer} setting minimum MAPQ of reads to be included into analysis.
#' @param mapqMax \code{integer} setting maximum MAPQ of reads to be included into analysis.
#' @param max_read_size maximum read length.
#' @param ncores number of cores to use.
#'
#' @return \code{tibble} where each row corresponds to sample - region combination. 
#' Columns of the returned \code{tibble} represent:
#' \describe{
#'   \item{SampleName}{Name of a sample as defined by \code{samplenames}.}
#'   \item{seqnames, start, end}{seqnames (or chromosomes) and coordinates of corresponding regions.}
#'   \item{names}{the same as \code{names(regions)} if it is not \code{NULL} or index of a corresponding region in \code{regions}.}
#'   \item{nFragsFetched}{number of fetched fragments before filtering.}
#'   \item{nFragsNonUnique}{number of non-unique fragments as defined by coordinates and methylation profiles.}

#'   \item{duplStatsMatrix}{\code{matrix} with 2 columns. First column represents number of fragments and second column represents number of times it is duplicated. 
#'   For example, row with numbers 1 at the first column and 100 at the second column should be interpreted as following: 100 fragments in the data has exactly 1 unique match. In other words, 100 fragments are unique.
#'   Another example, row with numbers 5 at the first column and 200 at the second column means that 5 unique fragments are duplicated 200 times in the data.
#'   }
#'   \item{duplFragNamesList}{\code{list} containing names of duplicated fragments. Names of elements of this list are data encoded as ...
#'   }
#' }
#' @export
#'

get_dupl_stats_from_bams <- function(bamfiles,
                                     samplenames,
                                     regions,
                                     genome,
                                     #whichContext = c("GCH","WCG","bisC","otherC", "allC"),
                                     collapseBySample = TRUE,
                                     #remove_nonunique = TRUE,
                                     #clip_until_nbg  = 1L,
                                     #noclip_protect_frac_above = 0.90,
                                     #max_bisC_meth = 0.1,
                                     #min_bisC_size = 10,
                                     mapqMin = 0,
                                     mapqMax = 255,
                                     max_read_size = 1000L,
                                     ncores = 1L){
  if(!inherits(bamfiles,"character"))
    stop("bamfiles must be vector of character")
  if(!inherits(samplenames,"character"))
    stop("samplenames must be vector of character")
  if(any(!file.exists(bamfiles))){
    nonex <- which(!file.exists(bamfiles))
    nonex <- paste0(bamfiles[nonex],collapse = "\n")
    stop(paste0("Could not find following BAM files:\n",nonex))
  }
  
  if(!inherits(regions,"GRanges"))
    stop("regions must be Granges")
  
  
  if(!inherits(genome,"character") & length(genome) !=1)
    stop("genome must be character and contain only one element")
  
  if(!file.exists(genome))
    stop("Could not find genome:", genome)
  
  
  ## remove metadata from regions
  GenomicRanges::mcols(regions) <- NULL
  
  ## load refsequences for regions
  ## extend regions by max_read_size from left and right
  
  fa_index <- Rsamtools::scanFaIndex(genome)
  GenomeInfoDb::seqlevels(regions) <- GenomeInfoDb::seqlevels(fa_index)
  GenomeInfoDb::seqinfo(regions) <- suppressMessages(GenomeInfoDb::Seqinfo(seqnames = as.character(GenomeInfoDb::seqnames(fa_index)),
                                                                           seqlengths = GenomicRanges::width(fa_index)))
  refseq_gr <- suppressWarnings(GenomicRanges::resize(regions,width = GenomicRanges::width(regions) + 2*max_read_size,fix="center"))
  
  ## trim to remove parts out of chromosomes
  refseq_gr <- GenomicRanges::trim(refseq_gr)
  refseqs <- Rsamtools::scanFa(genome,refseq_gr)
  
  # 2. split vector by sample name. If collapseBySample is FALSE add to samplenames some suffixes
  if(!collapseBySample){
    samplenames <- paste0("s",seq_along(samplenames),"_",samplenames)
  }
  
  bamfiles <- split(bamfiles, f=samplenames)
  # 3. for each set of bam files call C++ function
  # 
  duplStats_list <- do.call(rbind,lapply(names(bamfiles),
                                         function(bamgroupname){
                                           ctbl <- do.call(rbind,parallel::mclapply(seq_along(regions),
                                                                                    function(regi){
                                                                                      
                                                                                      regdf <- GenomicRanges::as.data.frame(regions[regi])
                                                                                      regdf$seqnames <- as.character(regdf[,"seqnames"])
                                                                                      seqgr <- as.data.frame(refseq_gr[regi])
                                                                                      
                                                                                      outlist <- fetch_dupl_stats_from_bams_cpp(infiles = bamfiles[[bamgroupname]],
                                                                                                                                regionChr = regdf[1,"seqnames"],
                                                                                                                                regionStart = regdf[1,"start"] - 1, # convert to 0-based
                                                                                                                                regionEnd = regdf[1,"end"],
                                                                                                                                
                                                                                                                                seqstring = as.character(refseqs[regi]),
                                                                                                                                seqStart = seqgr[1,"start"] - 1,
                                                                                                                                seqEnd = seqgr[1,"end"],
                                                                                                                                
                                                                                                                                mapqMin = mapqMin,
                                                                                                                                mapqMax = mapqMax)
                                                                                      ## output from Rcpp is like
                                                                                      # Rcpp::Named("nFragsFetched") = Rcpp::wrap(n_fetched),
                                                                                      # Rcpp::Named("nFragsNonUnique") = Rcpp::wrap(n_nonUnique),
                                                                                      # Rcpp::Named("duplStats") = Rcpp:wrap(duplStats4R),
                                                                                      # Rcpp::Named("duplFragNames" = EncToQnames4R)
                                                                                      
                                                                                      outtbl <- tibble::tibble_row("SampleName" = bamgroupname,
                                                                                                                   "seqnames" = regdf[1,"seqnames"],
                                                                                                                   "start" = regdf[1,"start"],
                                                                                                                   "end" = regdf[1,"end"],
                                                                                                                   "strand" = regdf[1,"strand"],
                                                                                                                   "names" = ifelse(!is.null(names(regions)), names(regions)[regi],regi),
                                                                                                                   "nFragsFetched" = outlist[["nFragsFetched"]],
                                                                                                                   "nFragsNonUnique" = outlist[["nFragsNonUnique"]],
                                                                                                                   "duplStatsMatrix" = list(outlist[["duplStats"]]),
                                                                                                                   "duplFragNamesList" = list(outlist[["duplFragNames"]]))
                                                                                      
                                                                                      return(outtbl)
                                                                                    },mc.cores = ncores))
                                           
                                           return(ctbl)
                                         }))
  
}