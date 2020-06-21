


#' Add to a \code{data.frame} with empirical joint frequencies columns with empirical joint probabilities as well as mutual information 
#'
#' @param emp_joint_frequencies count table which consists empirical joint frequencies. It must contain the columns with names \code{S, N00, N01, N10, N11}.
#'
#' @return \code{data.frame} with additional columns
#' @export
#'
#' @examples
#' 
#' \dontrun{
#' count_table_spacings <- add_joint_probs_and_MI(emp_joint_frequencies = count_table_spacings)
#' }
#' 
#' 
add_joint_probs_and_MI <- function(emp_joint_frequencies){
  stopifnot(all(c("S","N00","N01","N10","N11") %in% colnames(emp_joint_frequencies)) | 
              all(c("S","P00","P01","P10","P11") %in% colnames(emp_joint_frequencies)))
  
  # tot.counts <- rowSums(wind.count.tbl[,c("N_0_0","N_0_1","N_1_0","N_1_1",
  #                                         "N_0_NA","N_1_NA","N_NA_0","N_NA_1",
  #                                         "N_NA_NA")])
  if(!all(c("P00","P01","P10","P11") %in% colnames(emp_joint_frequencies))){
    tot.counts <- rowSums(emp_joint_frequencies[,c("N00","N01","N10","N11")])
    
    # wind.count.tbl$P_0_0 <- with(data = wind.count.tbl,
    #                              (N_0_0 + 0.5*N_0_NA + 0.5*N_NA_0 + 0.25 * N_NA_NA)/tot.counts)
    # wind.count.tbl$P_0_1 <- with(data = wind.count.tbl,
    #                              (N_0_1 + 0.5*N_0_NA + 0.5*N_NA_1 + 0.25 * N_NA_NA)/tot.counts)
    # wind.count.tbl$P_1_0 <- with(data = wind.count.tbl,
    #                              (N_1_0 + 0.5*N_1_NA + 0.5*N_NA_0 + 0.25 * N_NA_NA)/tot.counts)
    # wind.count.tbl$P_1_1 <- with(data = wind.count.tbl,
    #                              (N_1_1 + 0.5*N_1_NA + 0.5*N_NA_1 + 0.25 * N_NA_NA)/tot.counts)
    emp_joint_frequencies$P00 <- with(data = emp_joint_frequencies,
                                      (N00)/tot.counts)
    emp_joint_frequencies$P01 <- with(data = emp_joint_frequencies,
                                      (N01)/tot.counts)
    emp_joint_frequencies$P10 <- with(data = emp_joint_frequencies,
                                      (N10)/tot.counts)
    emp_joint_frequencies$P11 <- with(data = emp_joint_frequencies,
                                      (N11)/tot.counts)
  }
  indep.prob <- with(emp_joint_frequencies,
                     data.frame("P0any" = P00 + P01,
                                "P1any" = P10 + P11,
                                "Pany0" = P00 + P10,
                                "Pany1" = P01 + P11,
                                stringsAsFactors = F))
  
  emp_joint_frequencies <- cbind(emp_joint_frequencies,
                                 indep.prob)
  emp_joint_frequencies$MI <- with(emp_joint_frequencies,
                                   ifelse(P00 == 0,0,P00 * log(P00/(P0any * Pany0))) + 
                                     ifelse(P01 == 0,0,P01 * log(P01/(P0any * Pany1))) +
                                     ifelse(P10 == 0,0,P10 * log(P10/(P1any * Pany0))) +
                                     ifelse(P11 == 0,0,P11 * log(P11/(P1any * Pany1))))
  
  emp_joint_frequencies$H_ij <- -with(emp_joint_frequencies,
                                      ifelse(P00 == 0,0,P00 * log(P00)) + 
                                        ifelse(P01 == 0,0,P01 * log(P01)) +
                                        ifelse(P10 == 0,0,P10 * log(P10)) +
                                        ifelse(P11 == 0,0,P11 * log(P11)))
  
  emp_joint_frequencies$IQR <- with(emp_joint_frequencies,
                                    MI/H_ij)
  emp_joint_frequencies
}

