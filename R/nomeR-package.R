#' The 'nomeR' package.
#'
#' @description
#'
#' Regulation of gene expression by DNA-binding proteins, such as
#' sequence-specific transcription factors (TFs), is a key process driving
#' many biological phenotypes. While bulk assays such as ChIP-seq and ATAC-seq
#' have advanced our understanding of TF binding, they only provide
#' population-averaged information and cannot resolve stochastic variation
#' at single-DNA-molecule resolution.
#'
#' Single-molecule footprinting (SMF) technologies — including NOMe-seq,
#' SAMOSA, Fiber-seq, DAF-seq, FOODIE, and others — overcome this limitation
#' by recording the methylation state of individual DNA molecules, enabling
#' inference of protein-DNA interactions at single-molecule resolution.
#'
#' nomeR provides a Bayesian statistical framework for the analysis of
#' SMF data from any such technology. The package implements:
#' \itemize{
#'   \item Footprint Spectral Analysis (FSA): Bayesian inference of the
#'     footprint length spectrum and emission probabilities from a sample
#'     of single molecules.
#'   \item Genome-wide footprint prediction: per-molecule footprint
#'     configuration decoding and per-tile enrichment scoring using
#'     Hidden Markov Models with posterior or Viterbi decoding.
#' }
#' Both components accept continuous modification probabilities directly,
#' without requiring binarization of the input data.
#'
#' @useDynLib nomeR, .registration = TRUE
#' @import methods
#' @import Rcpp
#' @importFrom rstantools rstan_config
#' @importFrom RcppParallel RcppParallelLibs
#'
#' @references
#' Stan Development Team (NA). RStan: the R interface to Stan. R package
#' version 2.32.6. https://mc-stan.org
#'
"_PACKAGE"
