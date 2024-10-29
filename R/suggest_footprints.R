
#' Suggest footprints using inferred spectrum
#'
#' Function for getting potential footprints based on inferred abundances.
#' This function uses \code{stats::smooth.spline} to smooth spectrum and detect local extremums.
#' After local extremums has been detected it takes coordinates of local maxima and extends peaks
#' until signal drops by \code{max_abund_log2drop} or it reaches local minima or width of a peak becomes larger than \code{max_peak_width}.
#'
#'
#' @param S footprint lengths as returned by \code{\link{infer_footprints_vb}}.
#' Note that S=1 is reserved for background, therefore it is recommended to remove prior to running this function.
#' @param y mean or another measure of footprint abundance as returned \code{\link{infer_footprints_vb}}.
#' @param isLog \code{logical} indicating whether \code{y} is log-scaled.

#' @param spline_spar parameter controlling smoothness of spline ((0,1], the higher the smoother). Please see \code{\link[stats]{smooth.spline}}.
#' It is recommended to test different values for \code{spline_spar}, for example 0.1, 0.3, 0.5, 0.75 to check whether suggested footprints look as expected.
#'
#' @param max_abund_log2drop maximum decrease in abundance relative to value at local maxima until which peaks are extended.
#' @param max_peak_width maximum width of peaks in spectrum.
#' @param ... parameters for function \code{\link[stats]{smooth.spline}}.
#'
#' @return \code{list} containing following slots:
#'
#' \code{"ftp_ranges"} - \code{matrix} with suggested footprints.
#' \code{"smoothed_signal"} - smoothed spectrum as returned by  \code{\link[stats]{smooth.spline}} for each value in \code{S}.
#'
#' @export
#'
## @examples
suggest_footprints <- function(S,
                               y,
                               max_abund_log2drop = 1,
                               max_peak_width = Inf,
                               isLog = T,
                               spline_spar = 0.6,
                               ...){
  stopifnot(length(S) == length(y))

  ## find local extremums using spline.smooth ##
  ## smooth signal

  smspl <- stats::smooth.spline(x=S,y=y,
                                spar = spline_spar,
                                ...)
  ## smoothed signal
  smoothed_signal <- stats::predict(object = smspl,
                                    x = S)$y
  ## first derivative
  spline_deriv1 <- stats::predict(object = smspl,
                                  x = S,
                                  deriv=1)$y
  ## second derivative
  spline_deriv2 <- stats::predict(object = smspl,
                                  x = S,
                                  deriv=2)$y

  ## find x where spline intercepts 0, i.e. roots
  roots_x_idx <- which(spline_deriv1[1:(length(spline_deriv1)-1)] * spline_deriv1[2:length(spline_deriv1)] <=0)+1



  x_extremums <- S[roots_x_idx]
  y_extremums <- y[roots_x_idx]
  # type_extremums: 1 fox max, -1 for min
  type_extremums <- ifelse(spline_deriv2[roots_x_idx] <=0, 1,-1)

  # plot(x=S,y,type="p",pch=20)
  # lines(x=S,smoothed_signal,col="red")
  # points(x=x_extremums,y=y_extremums,col = ifelse(type_extremums==1,"green","blue"),pch=20)
  spectr_max_idx <- roots_x_idx[which(type_extremums == 1)]


  if(length(spectr_max_idx) > 0){
    ftp_ranges <- do.call(rbind,lapply(spectr_max_idx,
                                       function(idx){

                                         if(smoothed_signal[idx] == 0 & !isLog){
                                           stop(paste0("Incorrect value of local maximum at S=",S[idx],"; y=",smoothed_signal[idx]))
                                         }
                                         ## find left range
                                         left_range <- max(idx-1,1):1
                                         ## remove indices with negative derivative
                                         neg_der <- which(spline_deriv1[left_range] < 0)
                                         if(length(neg_der) > 0){
                                           left_range <- left_range[1:max(neg_der[1] - 1,1)]
                                         }

                                         ## remove indices where ratio to extremum is lower than max_abund_drop
                                         if(isLog){
                                           ratio_to_extrem <- smoothed_signal[left_range] - smoothed_signal[idx]
                                         } else {
                                           ratio_to_extrem <- log2(smoothed_signal[left_range]/smoothed_signal[idx])
                                         }

                                         below_cutoff <- which(ratio_to_extrem < -max_abund_log2drop)
                                         if(length(below_cutoff) > 0){
                                           left_range <- left_range[1:max(below_cutoff[1] - 1,1)]
                                         }

                                         ## keep indices where length is less than max_peak_width
                                         left_range <- left_range[which(S[left_range[1]] - S[left_range] + 1 <= max_peak_width/2)]
                                         # get most left index
                                         left_S <- S[min(left_range)]

                                         ## find right range
                                         right_range <- min(idx+1,length(S)):length(S)
                                         ## remove indices with positive derivative
                                         pos_der <- which(spline_deriv1[right_range] > 0)
                                         if(length(pos_der) > 0){
                                           right_range <- right_range[1:max(pos_der[1] - 1,1)]
                                         }

                                         ## remove indices where ratio to extremum is lower than max_abund_drop
                                         if(isLog){
                                           ratio_to_extrem <- smoothed_signal[right_range] - smoothed_signal[idx]
                                         } else {
                                           ratio_to_extrem <- log2(smoothed_signal[right_range]/smoothed_signal[idx])
                                         }

                                         below_cutoff <- which(ratio_to_extrem < -max_abund_log2drop)
                                         if(length(below_cutoff) > 0){
                                           right_range <- right_range[1:max(below_cutoff[1] - 1,1)]
                                         }
                                         ## keep indices where length is less than max_peak_width
                                         right_range <- right_range[which(S[right_range] - S[right_range[1]]  + 1 <= max_peak_width/2)]

                                         # get most right index
                                         right_S <- S[max(right_range)]

                                         return(c("min_ftp_length" = left_S,
                                                  "max_ftp_length" = right_S))
                                       }))
    row.names(ftp_ranges) <- paste0("ftp",1:nrow(ftp_ranges))
  } else {
    ftp_ranges <- matrix(nrow=0,ncol=2)
    colnames(ftp_ranges) <- c("min_ftp_length","max_ftp_length")
  }

  return(list("ftp_ranges" = ftp_ranges,
              "smoothed_signal" = smoothed_signal))
}
