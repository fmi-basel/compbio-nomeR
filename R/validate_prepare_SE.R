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
validate_prepare_SE <- function(se,
																assayName,
																threshUnmod,
																threshMod) {
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
	bin_protect_data <- do.call(rbind,lapply(1:ncol(mod_prob_assays),
																					 function(sidx){
																					 	read_naar <- mod_prob_assays[[sidx]]
																					 	## get M-indices of non-NAs
																					 	nonNA_data <- nnawhich(read_naar,arr.ind=TRUE)
																					 	## first column - positions (rows), second column - reads(columns)

																					 	## convert to binary protection
																					 	protect_vec <- ifelse(read_naar[nonNA_data] >= threshMod,0,
																					 												ifelse(read_naar[nonNA_data] < threshUnmod,1,
																					 															 NA))

																					 	## construct output
																					 	## sidx - index of sample in SE
																					 	## fidx_glob - unique index of fragment across all samples, as if they were cbinded
																					 	## fidx_sample - index of fragment for the current sample
																					 	## posidx_ref - index of rows in SE, corresponds to reference position stored in rowRanges(se)
																					 	## protect - binary protection data, 0 - accessible, 1 - protected
																					 	nonNA_data <- cbind("sidx"=rep(sidx,nrow(nonNA_data)),
																					 											"fidx_glob"=ncol(read_naar)*(sidx - 1) + nonNA_data[,2],
																					 											"fidx_sample"=nonNA_data[,2],
																					 											"posidx_ref" = nonNA_data[,1],
																					 											"protect" = protect_vec)

																					 	## remove those positions which did not pass thresholding and return
																					 	return(nonNA_data[which(!is.na(nonNA_data[,"protect"])),,drop=F])
																					 }))
	## add reference position
	bin_protect_data <- cbind(bin_protect_data,
														"refpos" = start(rowRanges(se))[bin_protect_data[,"posidx_ref"]])

	## add position within fragments
	frag_rstart <-tapply(bin_protect_data[,"refpos"],bin_protect_data[,"fidx_glob"],min,na.rm=T)

	## NOTE: the fragpos are 0 - based positions within fragments
	bin_protect_data <- cbind(bin_protect_data,
														"fragpos" = bin_protect_data[,"refpos"] - frag_rstart[bin_protect_data[,"fidx_glob"]])

	## order 1) by global fragment index fidx_glob; 2) by positions within fragments fragpos
	bin_protect_data <- bin_protect_data[order(bin_protect_data[,"fidx_glob"],
																						 bin_protect_data[,"fragpos"]),,drop=F]
	## sidx - index of sample in SE
	## fidx_glob - unique index of fragment across all samples, as if they were cbinded
	## fidx_sample - index of fragment for the current sample
	## posidx_ref - index of rows in SE, corresponds to reference position stored in rowRanges(se)
	## protect - binary protection data, 0 - accessible, 1 - protected
	## refpos - genomic position within a reference
	## fragpos - position within a frament, 0 - based
	return(bin_protect_data)
}
