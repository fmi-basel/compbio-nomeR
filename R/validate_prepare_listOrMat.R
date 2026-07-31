#' Checks and prepares data structure for the nomeR C++ prediction functions
#'
#' @param data \code{matrix} or \code{list} with SMF protection data.
#'              1 - accessible position; 0 - inaccessible position; NA - missing data.
#'              If \code{matrix}, rows are fragments and columns are positions.
#'              If \code{list}, each element is a numeric vector of protection for a fragment,
#'              and the names of the list elements are used as fragment names in the output.
#'
#' @return \code{list} with slots - data_list and fragnames
#'
#' @keywords internal
#' @noRd
#' @importFrom checkmate test_matrix test_list
validate_prepare_listOrMat <- function(data) {
    if (test_matrix(data, mode = "numeric",
                    any.missing = TRUE, all.missing = FALSE,
                    min.rows = 1, min.cols = 1)) {

        ## check if any values are not in [0,1] or NA
        if (any(!is.na(data) & (data < 0 | data > 1))) {
            stop("data must contain values in [0, 1] or NA")
        }
        ## check if row.names exist, if not set it to 1:nrow
        if (is.null(row.names(data))) {
            row.names(data) <- seq_len(nrow(data))
        }
        ## encode NAs as -1.0 sentinel (C++ checks p < 0 to detect NA positions)
        data[which(is.na(data), arr.ind = TRUE)] <- -1.0

        ## collect all positions (including -1.0 NA sentinels); C++ skips p < 0
        nonNA_data <- which(!is.na(data), arr.ind = TRUE, useNames = FALSE)
        ## 1st column - rows (fragments); 2nd column - columns (positions)
        nonNA_data <- cbind(nonNA_data,
                            data[nonNA_data])
        colnames(nonNA_data) <- c("fidx_glob", "colidx", "mod_prob")

        ## ensure integer types required by C++ IntegerVector parameters
        storage.mode(nonNA_data[, "fidx_glob"]) <- "integer"
        storage.mode(nonNA_data[, "colidx"])    <- "integer"

        ## add position within fragments
        nonNA_data <- cbind(nonNA_data,
                            "fragpos" = nonNA_data[, "colidx"])
        ## fidx_glob - unique index of fragment across all samples, as if they were cbinded
        ## "colidx" - column index in the input matrix data
        ## mod_prob - modification probability in [0,1]; high values indicate accessible positions
        ## fragpos - position within a frament, 1 - based
        fragnames <- row.names(data)

    } else if (test_list(data, types = "numeric",
                         any.missing = TRUE, all.missing = TRUE,
                         min.len = 1)) {

        ## check if any values are not in [0,1] or NA
        if (any(!is.na(unlist(data, recursive = TRUE, use.names = FALSE)) &
                (unlist(data, recursive = TRUE, use.names = FALSE) < 0 |
                 unlist(data, recursive = TRUE, use.names = FALSE) > 1))) {
            stop("data must contain values in [0, 1] or NA")
        }

        ## check if names exist, if not set it to 1:length
        if (is.null(names(data))) {
            names(data) <- seq_len(length(data))
        }
        ## create a matrix with indices where 1st column is and index in the
        ## input list; 2nd column is an index within the vector; 3rd column is mod_prob
        nonNA_data <- do.call(rbind, lapply(
            seq_along(data),
            function(i){
                dvec <- data[[i]]
                ## encode NAs as -1.0 sentinel (C++ checks p < 0 to detect NA positions)
                dvec[is.na(dvec)] <- -1.0

                nonNA_data <- cbind("fidx_glob" = rep(i, length(dvec)),
                                    "colidx" = seq_len(length(dvec)),
                                    "mod_prob" = dvec)
                ## ensure integer types required by C++ IntegerVector parameters
                storage.mode(nonNA_data[, "fidx_glob"]) <- "integer"
                storage.mode(nonNA_data[, "colidx"])    <- "integer"
                return(nonNA_data)
            }))
        ## add position within fragments
        nonNA_data <- cbind(nonNA_data,
                            "fragpos" = nonNA_data[, "colidx"])
        fragnames <- names(data)

    } else {
        stop("'data' must be 'matrix' or 'list'")
    }

    ## order 1) by global fragment index fidx_glob; 2) by positions within fragments fragpos
    nonNA_data <- nonNA_data[order(nonNA_data[, "fidx_glob"],
                                   nonNA_data[, "fragpos"]), , drop = FALSE]
    return(list("nonNA_data" = nonNA_data,
                "fragnames" = fragnames))
}
