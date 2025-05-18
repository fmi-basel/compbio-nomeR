#' Generate in-silico single-molecule footprinting dataset
#'
#' @param region_len length of the region.
#' @param n_reads number of molecule to generate.
#' @param footprint_models \code{list} describing footprints which should be 
#'     generated. Each footprint within footprint_models is also a \code{list} 
#'     that must contain:
#'     \describe{
#'     \item{\code{PROTECT_PROB}}{vector of emission probabilities 
#'     with same length and footprint length}
#'     \item{\code{COVER_PROB}}{\code{numeric} scalar
#'     that describes footprint coverage}
#'     \item{\code{NAME}}{\code{character} scalar with 
#'     the name of the footprint}
#'     }
#' @param ftp_emission_posprob \code{matrix} describing positional probabilities 
#'     for each footprint. Number of rows must be equal to number of footprints 
#'     and number of columns must be equal to \code{region_len}.
#' @param bgprotectprob \code{numeric} emission probability for background.
#' @param infposdens if scalar between 0 and 1 treated as percentage of 
#'     informative positions. If a vector of integers treated as predefined 
#'     informative positions.
#'
#' @return \code{list} that contains
#' \describe{
#' \item{"DATA"}{\code{matrix} with in-silico data}
#' \item{"TRUE_CONF"}{true positions for each footprint}
#' \item{"ftp_pos_prob"}{normalized probabilities of placing each footprint}
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' ## length of a chromosome
#' chr_len <- 1000
#'
#' ftp_cover_c2 <- c(0.05, 0.5)
#' ftp_len_c2 <- c(50, 150)
#' case2_pos_mat <- rbind(matrix(c(rep(0, 5150), 1, rep(0, 4849)),
#'                               nrow = 1,
#'                               ncol = chr_len),
#'                        matrix(1 / chr_len,
#'                               nrow = 1,
#'                               ncol = chr_len))
#'
#' ## sequencing depths
#' max_nreads <- 30
#' ## footprint model
#' case2_ftp_models <- list(
#' "ftp50" = list("PROTECT_PROB" = rep(1,ftp_len_c2[1]),
#'                "COVER_PRIOR" = ftp_cover_c2[1],
#'                "NAME" = "ftp50"),
#' "ftp150" = list("PROTECT_PROB" = rep(1,ftp_len_c2[2]),
#'                 "COVER_PRIOR" = ftp_cover_c2[2],
#'                 "NAME" = "ftp150"))
#'
#' smf_data <-
#'   generate_insilico_SMF_data(region_len = chr_len,
#'                              n_reads = max_nreads,
#'                              footprint_models = case2_ftp_models,
#'                              ftp_emission_posprob = case2_pos_mat,
#'                              bgprotectprob = 0,
#'                              infposdens = 1
#'                              )
#'}
#'
generate_insilico_SMF_data <- function(region_len, # length of the amplicon
                                       n_reads, # number of reads
                                       footprint_models,
                                       ftp_emission_posprob = NULL,
                                       bgprotectprob,
                                       infposdens = 1) {
    
    extend_ampl <- FALSE
    #### CHECK ARGUMENTS ####
    checkmate::assert_integerish(region_len, len = 1, any.missing = FALSE)
    checkmate::assert_integerish(n_reads, len = 1, any.missing = FALSE)
    if (checkmate::testNumber(infposdens, lower = 0.0, upper = 1.0, 
                              null.ok = FALSE)) {
        checkmate::assert_number(infposdens, lower = 0.0, upper = 1.0, 
                                 null.ok = FALSE)
    } else {
        checkmate::assert_integerish(infposdens, lower = 1, upper = region_len,
                                     any.missing = FALSE, min.len = 1,
                                     unique = TRUE, null.ok = FALSE)
    }
    checkmate::assert_list(footprint_models,types = "list",
                           any.missing = FALSE, min.len = 1)
    lapply(footprint_models,
           function(x) {
               checkmate::assert_numeric(x$PROTECT_PROB, lower = 0.0,
                                         upper = 1.0, finite = TRUE, 
                                         null.ok = FALSE, any.missing = FALSE,
                                         min.len = 2)
               checkmate::assert_number(x$COVER_PRIOR, lower = 0.0,
                                        upper = 1.0, finite = TRUE, 
                                        null.ok = FALSE)
               checkmate::assert_character(x$NAME, min.chars = 1,
                                           any.missing = FALSE, len = 1)
           })
    
    checkmate::assert_matrix(ftp_emission_posprob,
                             mode = "numeric",
                             any.missing = FALSE,
                             nrows = length(footprint_models),
                             ncol = region_len,
                             null.ok = TRUE)
    
    checkmate::assert_number(bgprotectprob, lower = 0.0, upper = 1.0,
                             null.ok = FALSE, finite = TRUE)
    
    #### GET INFORMATIVE POSITIONS ####
    
    if (length(infposdens) == 1 && all(infposdens > 0) && 
        all(infposdens <= 1)) {
        # replace = FALSE makes sure that positions will not repeat
        inf_pos <- sort(sample(seq_len(region_len), 
                               size = floor(infposdens * region_len),
                               replace = FALSE))
    } else if (is.vector(infposdens,"integer") || 
               is.vector(infposdens,"numeric")) {
        infposdens <- as.integer(infposdens)
        if (any(!infposdens %in% seq_len(region_len))) {
            warnings("infposdens contains positions outsite amplicon ", 
                     "length. Removing such positions")
        }
        inf_pos <- infposdens[infposdens %in% seq_len(region_len)]
    } else {
        stop("infposdens must be a scalar between 0 and 1 or a vector of ", 
             "integers or numerics.")
    }
    
    #### PROCESSING FOOTPRINT MODELS #####
    ## check whether any of the footprint has length 1 or name "BG"
    if (any(vapply(footprint_models, function(x) {
        x$NAME
    }, NA_character_) == "BG") || any(vapply(footprint_models, function(x) {
        length(x$PROTECT_PROB)
    }, NA_integer_) == 1)) {
        stop("Footprint with name BG or length 1 are reserved for background. ", 
             "Remove them from footprint models.")
    }
    
    
    ## adding BG
    bgcoverprior <- 1 - sum(vapply(footprint_models, function(x) {
        x$COVER_PRIOR
    }, FUN.VALUE = NA_real_, USE.NAMES = FALSE))
    
    footprint_models <- c(list(list("NAME" = "BG",
                                    "COVER_PRIOR" = bgcoverprior,
                                    "PROTECT_PROB" = bgprotectprob)),
                          footprint_models)
    
    model_names <- vapply(footprint_models,function(x) {
        x$NAME
    }, NA_character_)
    
    if (any(duplicated(model_names))) {
        stop("Footprint names must not be duplicated!")
    }
    
    names(footprint_models) <- model_names
    
    model_cover_priors <- vapply(footprint_models, function(x) {
        x$COVER_PRIOR
    }, FUN.VALUE = NA_real_, USE.NAMES = FALSE)
    
    model_lengths <- vapply(footprint_models, function(x) {
        length(x$PROTECT_PROB)
    }, FUN.VALUE = NA_integer_, USE.NAMES = FALSE)
    
    model_start_priors <- .cover_prior2start_prior(model_cover_priors,
                                                   model_lengths)
    
    ## decide whether start of footprints will be outside amplicon
    if (extend_ampl) {
        max.model.len <- 2 * max(model_lengths)
    } else {
        max.model.len <- 1
    }
    
    reglen_ext <- region_len + max.model.len - 1
    
    ### assign footprint start probabilities for each position of the amplicon
    if (is.null(ftp_emission_posprob)) {
        ftp_emission_posprob <- matrix(data = 1/reglen_ext,
                                       nrow = length(model_names),
                                       ncol = reglen_ext)
    } else {
        ## add first row for background
        ftp_emission_posprob <- rbind(rep(1 / region_len,region_len),
                                      ftp_emission_posprob)
        
        ## if start of ftps outside amplicon boundaries allowed, add 
        ## additional probs to the left and renormalize
        if (extend_ampl) {
            ext_mat <- matrix(data = ftp_emission_posprob[, 1],
                              nrow = nrow(ftp_emission_posprob),
                              ncol = max.model.len - 1, byrow = FALSE)
            ftp_emission_posprob <- cbind(ext_mat,
                                          ftp_emission_posprob)
        } else {
            ## if extention of amplicon is not allowed, set positional 
            ## probabilities of footprints which do not fit at some positions 
            ## to 0
            ftp_emission_posprob <- do.call(
                rbind, lapply(
                    seq_len(nrow(ftp_emission_posprob)),
                    function(ftp) {
                        mod_len <- model_lengths[ftp]
                        prvec <- ftp_emission_posprob[ftp, ]
                        zero_start_pos <- reglen_ext - mod_len + 2
                        if (zero_start_pos <= reglen_ext) {
                            prvec[seq(reglen_ext - mod_len + 1, reglen_ext)] <- 0
                        }
                        return(prvec)
                    }))
        }
    }
    
    ftp_emission_posprob <- ftp_emission_posprob / rowSums(ftp_emission_posprob)
    
    ftp_emission_posprob <- apply(
        ftp_emission_posprob, MARGIN = 2,
        function(x) {
            x * model_start_priors / sum(x * model_start_priors)
        })
    SMF_data <- lapply(seq_len(n_reads), function(i) {
        dat <- matrix(NA, nrow = 1, ncol = reglen_ext)
        true_conf_data <- data.frame(
            seq = i,
            pos = seq_len(reglen_ext) - max.model.len + 1,
            model = NA_character_,
            stringsAsFactors = FALSE)
        pos <- 1
        while (pos <= reglen_ext) {
            ## choose which model to emit
            browser()
            emit_model <- sample(x = model_names,
                                 size = 1,
                                 prob = ftp_emission_posprob[, pos])
            
            true_conf_data$model[pos] <- emit_model
            
            model_len <- model_lengths[emit_model]
            model_prot_prob <- footprint_models[[emit_model]][["PROTECT_PROB"]]
            
            ## draw random data using binomial distribution
            model_synth_data <- stats::rbinom(n = model_len, size = 1, 
                                              prob = model_prot_prob)
            fill.idx <- seq(pos, min(pos + model_len - 1, reglen_ext))
            dat[1, fill.idx] <- model_synth_data[seq_len(length(fill.idx))]
            pos <- pos + model_len
            names(pos) <- NULL
        }
        
        dat.out <- list(
            "DATA" = dat,
            "TRUE_CONF" = true_conf_data[!is.na(true_conf_data$model), ]
        )
        return(dat.out)
    })
    
    SMF_sim_data <- do.call(rbind, lapply(SMF_data,function(x) {
        x[["DATA"]]
    }))
    true_conf <- do.call(rbind, lapply(SMF_data,function(x) {
        x[["TRUE_CONF"]]
    }))
    ## remove additional flanks
    SMF_sim_data <- SMF_sim_data[, seq(max.model.len, reglen_ext), drop = FALSE]
    ## set non GpC positions to NA
    SMF_sim_data[, !(seq_len(ncol(SMF_sim_data)) %in% inf_pos)] <- NA
    
    ## return without additional flanking columns
    list("DATA" = SMF_sim_data,
         "TRUE_CONF" = true_conf,
         "ftp_pos_prob" = ftp_emission_posprob)
}
