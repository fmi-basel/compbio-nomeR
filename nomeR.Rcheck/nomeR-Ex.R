pkgname <- "nomeR"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('nomeR')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("generate_insilico_SMF_data")
### * generate_insilico_SMF_data

flush(stderr()); flush(stdout())

### Name: generate_insilico_SMF_data
### Title: Generate an in-silico single-molecule footprinting dataset
### Aliases: generate_insilico_SMF_data

### ** Examples

## length of a chromosome
chr_len <- 1000

ftp_cover_c2 <- c(0.05, 0.5)
ftp_len_c2 <- c(50, 150)
case2_pos_mat <- rbind(matrix(c(rep(0, 515), 1, rep(0, 484)),
                              nrow = 1,
                              ncol = chr_len),
                       matrix(1 / chr_len,
                              nrow = 1,
                              ncol = chr_len))

## sequencing depths
max_nreads <- 30
## footprint model
case2_ftp_models <- list(
"ftp50" = list("PROTECT_PROB" = rep(1,ftp_len_c2[1]),
               "COVER_PRIOR" = ftp_cover_c2[1],
               "NAME" = "ftp50"),
"ftp150" = list("PROTECT_PROB" = rep(1,ftp_len_c2[2]),
                "COVER_PRIOR" = ftp_cover_c2[2],
                "NAME" = "ftp150"))

smf_data <-
  generate_insilico_SMF_data(region_len = chr_len,
                             n_reads = max_nreads,
                             footprint_models = case2_ftp_models,
                             ftp_emission_posprob = case2_pos_mat,
                             bgprotectprob = 0,
                             infposdens = 1
                             )




cleanEx()
nameEx("get_ftp_inference_summary")
### * get_ftp_inference_summary

flush(stderr()); flush(stdout())

### Name: get_ftp_inference_summary
### Title: Extract and/or plot estimates from footprint inference and
###   optionally suggest footprints
### Aliases: get_ftp_inference_summary

### ** Examples


## Simple data with two footprints of lengths 5 and 10 bp.
## The table below is a count table of observed occurrences of
## 00, 01, 10, and 11 at spacings from 1 to 15.

ftp_5_10_data <- data.frame("S" = 1:15,
                            "N00" = c(1626964, 1508381,
                                      1420066, 1336897,
                                      1258679, 1185045,
                                      1136784, 1090029,
                                      1045066, 1001670,
                                      959927, 983020,
                                      1000410, 1012326,
                                      1019633),
                            "N01" = c(0, 113856,
                                      197657, 276539,
                                      350664, 420399,
                                      464873, 507988,
                                      549466, 589487,
                                      627997, 601629,
                                      580965, 565776,
                                      555217),
                            "N10" = c(0, 113921,
                                      197854, 276862,
                                      351147, 421042,
                                      465716, 509005,
                                      550645, 590837,
                                      629533, 603328,
                                      582814, 567737,
                                      557220),
                            "N11" = c(873036, 758842,
                                      674423, 594702,
                                      519510, 448514,
                                      402627, 357978,
                                      314823, 273006,
                                      232543, 257023,
                                      275811, 289161,
                                      297930))

## variational inference for footprints in the data
inf <- infer_footprints_vb(cooc_ctable = ftp_5_10_data,
                           ftp_lengths = 2:15)

## get estimates and plot footprint spectrum
inference_summary_list <- get_ftp_inference_summary(inf, plot = TRUE)




cleanEx()
nameEx("infer_footprints_optim")
### * infer_footprints_optim

flush(stderr()); flush(stdout())

### Name: infer_footprints_optim
### Title: Find point estimates for footprint abundance using Stan
###   optimization
### Aliases: infer_footprints_optim

### ** Examples


## Simple data with two footprints of lengths 5 and 10 bps.
## The table below is a count table of observed occurrences of
## 00, 01, 10, and 11 at spacings from 1 until 15.

ftp_5_10_data <- data.frame("S" = 1:15,
                            "N00" = c(1626964, 1508381,
                                      1420066, 1336897,
                                      1258679, 1185045,
                                      1136784, 1090029,
                                      1045066, 1001670,
                                      959927, 983020,
                                      1000410, 1012326,
                                      1019633),
                            "N01" = c(0, 113856,
                                      197657, 276539,
                                      350664, 420399,
                                      464873, 507988,
                                      549466, 589487,
                                      627997, 601629,
                                      580965, 565776,
                                      555217),
                            "N10" = c(0, 113921,
                                      197854, 276862,
                                      351147, 421042,
                                      465716, 509005,
                                      550645, 590837,
                                      629533, 603328,
                                      582814, 567737,
                                      557220),
                            "N11" = c(873036, 758842,
                                      674423, 594702,
                                      519510, 448514,
                                      402627, 357978,
                                      314823, 273006,
                                      232543, 257023,
                                      275811, 289161,
                                      297930))


## finding MAP estimate
inf_output <- infer_footprints_optim(cooc_ctable = ftp_5_10_data,
                                     ftp_lengths = 2:15)

## plot footprint spectrum
get_ftp_inference_summary(inf_output, plot = TRUE)




cleanEx()
nameEx("infer_footprints_sampling")
### * infer_footprints_sampling

flush(stderr()); flush(stdout())

### Name: infer_footprints_sampling
### Title: Footprint spectral analysis using the No-U-Turn Sampler (NUTS)
###   implemented in Stan
### Aliases: infer_footprints_sampling

### ** Examples


## Simple data with two footprints of lengths 5 and 10 bps.
## The table below is a count table of observed occurrences of
## 00, 01, 10, and 11 at spacings from 1 until 15.

ftp_5_10_data <- data.frame("S" = 1:15,
                            "N00" = c(1626964, 1508381,
                                      1420066, 1336897,
                                      1258679, 1185045,
                                      1136784, 1090029,
                                      1045066, 1001670,
                                      959927, 983020,
                                      1000410, 1012326,
                                      1019633),
                            "N01" = c(0, 113856,
                                      197657, 276539,
                                      350664, 420399,
                                      464873, 507988,
                                      549466, 589487,
                                      627997, 601629,
                                      580965, 565776,
                                      555217),
                            "N10" = c(0, 113921,
                                      197854, 276862,
                                      351147, 421042,
                                      465716, 509005,
                                      550645, 590837,
                                      629533, 603328,
                                      582814, 567737,
                                      557220),
                            "N11" = c(873036, 758842,
                                      674423, 594702,
                                      519510, 448514,
                                      402627, 357978,
                                      314823, 273006,
                                      232543, 257023,
                                      275811, 289161,
                                      297930))


## HMC inference for footprints in the data
inf_output <- infer_footprints_sampling(cooc_ctable = ftp_5_10_data,
                                        ftp_lengths = 2:15, ncpu = 2)

## plot footprint spectrum
get_ftp_inference_summary(inf_output, plot = TRUE)




cleanEx()
nameEx("infer_footprints_vb")
### * infer_footprints_vb

flush(stderr()); flush(stdout())

### Name: infer_footprints_vb
### Title: Footprint spectral analysis by Variational Bayes inference using
###   Stan Automatic Differentiation Variational Inference (ADVI)
### Aliases: infer_footprints_vb

### ** Examples


## Simple data with two footprints of lengths 5 and 10 bps.
## The table below is a count table of observed occurrences of
## 00, 01, 10, and 11 at spacings from 1 until 15.

ftp_5_10_data <- data.frame("S" = 1:15,
                            "N00" = c(1626964, 1508381,
                                      1420066, 1336897,
                                      1258679, 1185045,
                                      1136784, 1090029,
                                      1045066, 1001670,
                                      959927, 983020,
                                      1000410, 1012326,
                                      1019633),
                            "N01" = c(0, 113856,
                                      197657, 276539,
                                      350664, 420399,
                                      464873, 507988,
                                      549466, 589487,
                                      627997, 601629,
                                      580965, 565776,
                                      555217),
                            "N10" = c(0, 113921,
                                      197854, 276862,
                                      351147, 421042,
                                      465716, 509005,
                                      550645, 590837,
                                      629533, 603328,
                                      582814, 567737,
                                      557220),
                            "N11" = c(873036, 758842,
                                      674423, 594702,
                                      519510, 448514,
                                      402627, 357978,
                                      314823, 273006,
                                      232543, 257023,
                                      275811, 289161,
                                      297930))


## VB inference for footprints in the data
inf_output <- infer_footprints_vb(cooc_ctable = ftp_5_10_data,
                                  ftp_lengths = 2:15)

## plot footprint spectrum
get_ftp_inference_summary(inf_output, plot = TRUE)




cleanEx()
nameEx("predict_footprints")
### * predict_footprints

flush(stderr()); flush(stdout())

### Name: predict_footprints
### Title: Calculate posterior probabilities and predict footprints in
###   single-molecule footprinting (SMF) data
### Aliases: predict_footprints

### ** Examples

set.seed(3346)
nc <- 50
nr <- 50
rmatr <- matrix(data = as.integer(rnorm(nc * nr) >= 0.5),
                ncol = nc,nrow=nr)

## create dummy footprints
bg.pr <- 0.5
ft.pr <- 1-bg.pr
ft.len <- 15

## creating a list of binding models for nomeR
ftp.models <- list(list("PROTECT_PROB" = rep(0.99,ft.len),
                        "COVER_PRIOR" = ft.pr,
                        "NAME" = "FOOTPRINT"))

nomeR.out <- predict_footprints(data=rmatr,
                                footprint_models = ftp.models,
                                bgprotectprob = 0.05,
                                bgcoverprior = bg.pr)




### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
