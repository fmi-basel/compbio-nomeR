#' Checks and prepares data structure for c++ run_cpp_nomeR function. version 2
#'
#' @param data \code{matrix} or \code{list} with NOMe-seq data
#'
#' @return \code{list} with slots - data_list and fragnames
#'
#' @keywords internal
#' @noRd
#' @importFrom checkmate test_matrix test_list
validate_data_for_predict2 <- function(data) {


	## the C++ function requires: vector with unique integer fragment IDs, 0-based positions within fragments and vector with binary protection.
	# out.list <- run_cpp_nomeR(protect_data[,"fidx_glob"], ## unique fragment ID (actually index in the original matrix)
	# 													protect_data[,"fragpos"],      ## position within fragment, 0 - based
	# 													protect_data[,"protect"],   ## binary protection data, 0 - accessible, 1 - protected
	# 													footprint_models,
	# 													bgprotectprob,
	# 													start_priors["BG"],
	# 													report_prediction_in_flanks,
	# 													ncpu,
	# 													verbose)


	if (test_matrix(data, mode = "integerish",
									any.missing = TRUE, all.missing = FALSE,
									min.rows = 1, min.cols = 1)) {

		## check if any values are not in c(0,1,NA)
		if (any(!(data %in% c(0, 1, NA)))) {
			stop("data must contain only 0, 1 or NA")
		}
		## check if row.names exist, if not set it to 1:nrow
		if (is.null(row.names(data))) {
			row.names(data) <- seq_len(nrow(data))
		}
		## find non-na elements
		nonNA_data <- which(!is.na(data),arr.ind=TRUE,useNames=F)
		## 1st column - rows (fragments); 2nd column - columns (positions)
		nonNA_data <- cbind(nonNA_data,
												data[nonNA_data])
		colnames(nonNA_data) <- c("fidx_glob","colidx","protect")
		## add position within fragments
		frag_rstart <-tapply(nonNA_data[,"colidx"],nonNA_data[,"fidx_glob"],min,na.rm=T)
		## NOTE: the fragpos are 0 - based positions within fragments
		nonNA_data <- cbind(nonNA_data,
												"fragpos" = nonNA_data[,"colidx"] - frag_rstart[nonNA_data[,"fidx_glob"]])

		## fidx_glob - unique index of fragment across all samples, as if they were cbinded
		## "colidx" - column index in the input matrix data
		## protect - binary protection data, 0 - accessible, 1 - protected
		## fragpos - position within a frament, 0 - based
		fragnames <- row.names(data)
	} else if (test_list(data,types = "integerish",
											 any.missing = TRUE, all.missing = TRUE,
											 min.len = 1)) {

		## check if any values are not in c(0,1,NA)
		if (any(!(unlist(data, recursive = TRUE, use.names = FALSE) %in%
							c(0, 1, NA)))) {
			stop("data must contain only 0, 1 or NA")
		}

		## check if names exist, if not set it to 1:length
		if (is.null(names(data))) {
			names(data) <- seq_len(length(data))
		}

		## create a matrix with indices where 1st column is and index in the input list; 2nd column is and index within the vector; 3rd column is protect
		nonNA_data <- do.call(rbind,lapply(seq_along(data),
												 function(i){
												 	dvec <- data[[i]]
												 	nna_idx <- which(!is.na(dvec))
												 	nonNA_data <- cbind("fidx_glob" = rep(i,length(nna_idx)),
												 											"colidx" = nna_idx,
												 											"protect" = dvec[nna_idx],
												 											"fragpos" = nna_idx - min(nna_idx,na.rm=T))
												 	return(nonNA_data)
												 }))
		fragnames <- names(data)
	} else {
		stop("'data' must be 'matrix' or 'list'")
	}

	## order 1) by global fragment index fidx_glob; 2) by positions within fragments fragpos
	nonNA_data <- nonNA_data[order(nonNA_data[,"fidx_glob"],
																 nonNA_data[,"fragpos"]),,drop=F]
	return(list("nonNA_data" = nonNA_data,
							"fragnames" = fragnames))
}
