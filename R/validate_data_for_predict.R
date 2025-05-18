#' Checks and prepares data structure for c++ run_cpp_nomeR function
#'
#' @param data \code{matrix} or \code{list} with NOMe-seq data
#'
#' @return \code{list} with slots - data_list and fragnames
#'
#' @keywords internal
#' @noRd
#' @importFrom checkmate test_matrix test_list
validate_data_for_predict <- function(data) {
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
        ## convert NAs to 2
        data[which(is.na(data), arr.ind = TRUE)] <- 2
        ## creating list of rows and fragnames
        data_list <- lapply(seq_len(nrow(data)), function(i) data[i, ])
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
        ## convert NAs to 2
        data <- lapply(data, function(x) {
            x[is.na(x)] <- 2
            return(x)
        })
        
        data_list <- data
        fragnames <- names(data)
    } else {
        stop("'data' must be 'matrix' or 'list'")
    }
    
    return(list("data_list" = data_list,
                "fragnames" = fragnames))
}