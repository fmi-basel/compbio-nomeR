#' The 'nomeR' package.
#'
#' @description
#'
#' Regulation of gene expression by DNA binding proteins, such as sequence 
#' specific transcription factor (TF), is a key step which drives a plethora of 
#' biological processes.
#' During the past decade, remarkable achievements have been made in 
#' understanding how TFs regulate their targets  using ChIP-seq and ATAC-seq 
#' assays.
#' However, these assays provide only averaged information of larger cell 
#' populations. Hence, such assays are unable to detect states of binding at 
#' single-DNA-molecule resolution.
#' In contrast, emerging single-molecule footprinting (SMF) technologies, such 
#' as the Nucleosome Occupancy and Methylome Sequencing (NOMe-seq) assay, 
#' provide information about states of binding sites at single-DNA-molecule 
#' resolution and allow better molecular characterization of inherently 
#' stochastic processes of protein-DNA interactions.
#'
#' Despite the many advantages that single molecule footprinting technologies 
#' bring for studying mechanisms of gene regulation, rigorous statistical 
#' models for the analysis of datasets generated using these technologies are 
#' still missing.
#' Here, we present, to our knowledge, the first statistical model for 
#' inference and predictions of footprints for unbiased quantitative analysis 
#' and interpretation of SMF datasets.
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
