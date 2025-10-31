#' Checks and prepares data structure for c++ run_cpp_nomeR function
#'
#' @param data \code{matrix} or \code{list} with NOMe-seq data
#'
#' @return \code{list} with slots - data_list and fragnames
#'
#' @keywords internal
#' @noRd
#' @importFrom checkmate test_matrix test_list
validate_prepare_listOrMat <- function(data) {
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
		## convert NAs to 2. For matrix we keep all positions, even with NAs as starts and ends must match dimensions of input matrix
		data[which(is.na(data), arr.ind = TRUE)] <- 2
		
		## find non-na elements
		nonNA_data <- which(!is.na(data),arr.ind=TRUE,useNames=F)
		## 1st column - rows (fragments); 2nd column - columns (positions)
		nonNA_data <- cbind(nonNA_data,
												data[nonNA_data])
		colnames(nonNA_data) <- c("fidx_glob","colidx","protect")
		
		## add position within fragments
		nonNA_data <- cbind(nonNA_data,
												"fragpos" = nonNA_data[,"colidx"])
		## fidx_glob - unique index of fragment across all samples, as if they were cbinded
		## "colidx" - column index in the input matrix data
		## protect - binary protection data, 0 - accessible, 1 - protected
		## fragpos - position within a frament, 1 - based
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
																			 	## convert NAs to 2
																			 	dvec[is.na(dvec)] <- 2
																			 	
																			 	nonNA_data <- cbind("fidx_glob" = rep(i,length(dvec)),
																			 											"colidx" = 1:length(dvec),
																			 											"protect" = dvec)
																			 	return(nonNA_data)
																			 }))
		## add position within fragments
		nonNA_data <- cbind(nonNA_data,
												"fragpos" = nonNA_data[,"colidx"])
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