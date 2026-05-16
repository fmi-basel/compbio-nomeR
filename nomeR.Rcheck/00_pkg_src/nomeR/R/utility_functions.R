

.is_profiling_enabled <- function(profile) {
    isTRUE(profile) || isTRUE(getOption("nomeR.profile", FALSE))
}


.time_block <- function(expr, label, timings_env, enabled, cli_report = TRUE) {
    if (!enabled) {
        return(force(expr))
    }

    t0 <- proc.time()[["elapsed"]]
    out <- force(expr)
    dt <- proc.time()[["elapsed"]] - t0

    timings_env[[label]] <- dt

    if (cli_report) {
        cli::cli_inform(
            "{.emph {label}} took {formatC(dt, digits = 3, format = 'f')} s"
        )
    }

    out
}

## function for converting vector of start priors to cover priors
.start_prior2cover_prior <- function(start_prior,
                                     footprint_len) {
    stopifnot(length(start_prior) == length(footprint_len))

    cover_prior <- start_prior * footprint_len
    cover_prior / sum(cover_prior)
}

## function for converting vector of covert priors to start priors
.cover_prior2start_prior <- function(cover_prior,
                                     footprint_len) {
    stopifnot(length(cover_prior) == length(footprint_len))
    stopifnot(all(footprint_len > 0))

    start_prior <- cover_prior / footprint_len
    start_prior / sum(start_prior)
}

.message_timestamp <- function(msg) {
    message("[", strftime(Sys.time()), "]: ", msg)
}

.warning_timestamp <- function(msg){
    warning("[", strftime(Sys.time()), "]: ", msg)
}

.onUnload <- function(libpath) {
    library.dynam.unload("nomeR", libpath)
}

.get_ftp_annotation <- function(ftp_len_mat,
                                      ftp_spec,
                                      bg_cover) {

    ftp_anno <- do.call(rbind,lapply(seq_len(nrow(ftp_len_mat)),
                       function(ridx) {

                           ## get middle points
                           ftplen <- seq(from = ftp_len_mat[ridx, 1],
                                         to = ftp_len_mat[ridx, 2],
                                         by = ftp_len_mat[ridx, 3])
                           if (ftp_len_mat[ridx, 3] == 1) {
                               ftpcov <- ftp_spec$mean[match(ftplen,
                                                             ftp_spec$ftp_length)]
                           } else {
                               intstep <- round(ftp_len_mat[ridx, 3] / 2)
                               curftp_spec <- ftp_spec[(ftp_spec$ftp_length >=
                                                            min(ftplen) - intstep + 1) &
                                                           (ftp_spec$ftp_length <=
                                                                max(ftplen) + intstep),
                                                       , drop = FALSE]
                               curftp_spec$group <- cut(curftp_spec$ftp_length,
                                                        breaks = length(ftplen),
                                                        labels = ftplen,
                                                        include.lowest = TRUE
                                                        )
                               ftpcov <- tapply(curftp_spec$mean, curftp_spec$group,
                                                sum)
                           }
                           ftp_anno <- data.frame(
                               ftp_name = paste0("ftp--", ftplen),
                               ftp_group = row.names(ftp_len_mat)[ridx],
                               ftp_length = ftplen,
                               ftp_cover = ftpcov)
                           return(ftp_anno)
                       }))
    missing_len <- ftp_anno$ftp_length[is.na(ftp_anno$ftp_cover)]
    if (length(missing_len) > 0) {
        stop("The following footprint lengths are missing from ftp_spectrum: ",
             paste(missing_len, collapse = ", "),
             ". Ensure ftp_len_mat only references lengths present in the spectrum.")
    }
    cover_sum <- sum(ftp_anno$ftp_cover)
    if (cover_sum == 0) {
        stop("Sum of footprint coverages in ftp_spectrum is 0. ",
             "Cannot normalise cover priors.")
    }
    ftp_anno$ftp_cover <- ftp_anno$ftp_cover / cover_sum * (1 - bg_cover)

    return(ftp_anno)
}
