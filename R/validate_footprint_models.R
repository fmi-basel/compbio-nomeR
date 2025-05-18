#' @importFrom checkmate check_list assert_names assert_numeric assert_number
#'     assert_character
validate_footprint_models <- function(footprint_models,
                                      bgprotectprob,
                                      bgcoverprior,
                                      verbose,
                                      add = NULL) {
    ### validate binding_models
    if (check_list(footprint_models, types = "list",
                   any.missing = FALSE, min.len = 1)) {
        ftpnames <- vapply(footprint_models, function(x) {
            assert_names(
                names(x), 
                must.include = c("PROTECT_PROB", "COVER_PRIOR", "NAME"),
                add = add)
            assert_numeric(
                x = x[["PROTECT_PROB"]], lower = 0, upper = 1, 
                any.missing = FALSE, add = add)
            assert_number(
                x = x[["COVER_PRIOR"]], lower = 0, upper = 1, na.ok = FALSE, 
                add = add)
            return(x[["NAME"]])
        }, NA_character_)
        assert_character(x = ftpnames, any.missing = FALSE,
                         unique = TRUE, add = add)
    }
    
    ### validate bgprotectprob
    assert_number(x = bgprotectprob, na.ok = FALSE, lower = 0, 
                  upper = 1, add = add)
    ### validate bgcoverprior
    assert_number(x = bgcoverprior, na.ok = FALSE, lower = 0, 
                  upper = 1, add = add)
    
    ## convert cover priors to start priors and add them into binding_models
    if (verbose) {
        .message_timestamp(c("convert cover priors to start priors and add ", 
                             "them into binding_models..."))
    }
    cover_priors <- c("BG" = bgcoverprior,
                      vapply(footprint_models, function(x) { 
                          x[["COVER_PRIOR"]]
                      }, NA_real_))
    footpr_lens <- c("BG" = 1,
                     vapply(footprint_models, function(x) {
                         length(x[["PROTECT_PROB"]])
                     }, NA_integer_))
    start_priors <- .cover_prior2start_prior(cover_prior = cover_priors,
                                             footprint_len = footpr_lens)
    
    footprint_models <- sapply(seq_along(footprint_models),
                               function(i) {
                                   ft.model <- footprint_models[[i]]
                                   c(ft.model,
                                     list("PRIOR" = start_priors[i + 1]))
                               },
                               simplify = FALSE, USE.NAMES = TRUE)
    
    return(list("footprint_models" = footprint_models,
                "start_priors" = start_priors))
}
