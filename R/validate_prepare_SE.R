#' Checks Summarized experiment provided by footprointR
#' and prepares data structure for c++ run_cpp_nomeR function
#'
#' @param data \code{matrix} or \code{list} with NOMe-seq data
#'
#' @return \code{list} with slots - data_list and fragnames
#'
#' @keywords internal
#' @noRd
#' @importFrom SummarizedExperiment SummarizedExperiment rowRanges colData assay assayNames
#' @importFrom SparseArray NaArray nnawhich
#' @importFrom GenomicRanges GPos match seqnames start end
#' @importFrom IRanges subsetByOverlaps
#' @importFrom S4Vectors DataFrame SimpleList metadata
#' @import data.table
validate_prepare_SE <- function(se,
																assayName,
																threshUnmod,
																threshMod) {
	
	protect = mod_prob = fidx_sample = posidx_ref = refpos = fidx_glob = ftp_group = NULL # due to NSE notes in R CMD check
	
	
	### The code for checking the vailidity of se is copied
	### from the footprintR package developed by Charlotte Soneson and Michael Stadler
	
	stopifnot(is(se, "SummarizedExperiment"))
	
	stopifnot(!is.null(assayNames(se)) &&
							all(assayNames(se) != "") &&
							!any(duplicated(assayNames(se))))
	
	if (nrow(se) > 0) {
		stopifnot(!is.null(rownames(se)) &&
								!any(duplicated(rownames(se))))
	}
	stopifnot(!is.null(metadata(se)$readLevelData) &&
							is.list(metadata(se)$readLevelData) &&
							all(c("assayNames", "colDataColumns") %in%
										names(metadata(se)$readLevelData)))
	
	stopifnot("sample" %in% colnames(colData(se)))
	
	stopifnot(assayName %in% assayNames(se))
	
	## extract mod_prob and convert to binary
	
	mod_prob_assays <- assay(se,"mod_prob")
	bin_protect_data <- data.table::rbindlist(lapply(1:ncol(mod_prob_assays),
																									 function(sidx){
																									 	
																									 	read_naar <- mod_prob_assays[[sidx]]
																									 	## get M-indices of non-NAs
																									 	nonNA_data <- nnawhich(read_naar,arr.ind=TRUE)
																									 	colnames(nonNA_data) <- c("posidx_ref","fidx_sample")
																									 	## convert to data.table
																									 	nonNA_data <- as.data.table(nonNA_data)
																									 	## first column - positions (rows), second column - reads(columns)
																									 	
																									 	## add modprob
																									 	nonNA_data <- nonNA_data[,"mod_prob" := read_naar[as.matrix(nonNA_data)]]
																									 	
																									 	## convert to binary protection
																									 	nonNA_data <- nonNA_data[,protect := ifelse(mod_prob >= threshMod,0,
																									 																							ifelse(mod_prob < threshUnmod,1,
																									 																										 NA))]
																									 	## construct output
																									 	## sidx - index of sample in SE
																									 	## fidx_glob - unique index of fragment across all samples, as if they were cbinded
																									 	## fidx_sample - index of fragment for the current sample
																									 	## posidx_ref - index of rows in SE, corresponds to reference position stored in rowRanges(se)
																									 	## protect - binary protection data, 0 - accessible, 1 - protected
																									 	nonNA_data <- nonNA_data[,c("sidx","fidx_glob") := list(rep(sidx,nrow(nonNA_data)),
																									 																											 ncol(read_naar)*(sidx - 1) + fidx_sample)]
																									 	
																									 	## remove those positions which did not pass thresholding and return
																									 	nonNA_data <- nonNA_data[!is.na(protect)]
																									 	
																									 	return(nonNA_data)
																									 }))
	
	## add reference position
	bin_protect_data <- bin_protect_data[,"refpos" := start(rowRanges(se))[posidx_ref]]
	
	## add position within fragments
	## NOTE: the fragpos are 1 - based positions within fragments
	bin_protect_data <- bin_protect_data[,"fragpos" := refpos - min(refpos) + 1, list(fidx_glob)]
	
	## order by fidx_glob and fragpos by setting keyv
	setkeyv(bin_protect_data,cols = c("fidx_glob","fragpos"))
	
	## sidx - index of sample in SE
	## fidx_glob - unique index of fragment across all samples, as if they were cbinded
	## fidx_sample - index of fragment for the current sample
	## posidx_ref - index of rows in SE, corresponds to reference position stored in rowRanges(se)
	## protect - binary protection data, 0 - accessible, 1 - protected
	## refpos - genomic position within a reference
	## fragpos - position within a frament, 1 - based
	
	## reorder columns
	setcolorder(bin_protect_data, c("sidx", "fidx_glob", "fidx_sample","posidx_ref",
																	"refpos","fragpos",
																	"mod_prob","protect"))
	return(bin_protect_data)
}
