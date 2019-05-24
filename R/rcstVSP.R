#'
#' @title Reconstruction points
#'
#' @description To reconstruct the V-Spline with \eqn{\theta}, which is gained from \code{getTheta()}.
#' @references Z. Cao, D. Bryant, C. Fox, T. Molten and M. Parry. "V-Spline: an Adaptive Smoothing Spline
#' for Trajectory Reconstruction" arXiv preprint arXiv:1803.07184 (2018).
#' @export
rcstVSP <- function(data_frame,theta,xout=NULL){
  sub_t <- data_frame$t
  rowlen<- length(sub_t)

  #tx   <- seq(sub_t[1],sub_t[rowlen],length = (rowlen*10))
  if(is.null(xout)) tx <- data_frame$t
  else tx <- xout
  
  ft   <- numeric(length(tx))
  dft  <- numeric(length(tx))
  d2ft <- numeric(length(tx))

  for(k in 1:length(tx)){
    x_est    <- tx[k]
    temp_x   <- 0
    temp_dx  <- 0
    temp_d2x <- 0

    temp_x   <- theta[1]*bN1(1,sub_t,x_est)+theta[2]*bN2(1,sub_t,x_est)
    temp_dx  <- theta[1]*dN1(1,sub_t,x_est)+theta[2]*dN2(1,sub_t,x_est)
    temp_d2x <- theta[1]*d2N1(1,sub_t,x_est)+theta[2]*d2N2(1,sub_t,x_est)

    for(l in 1:(length(sub_t)-2)){
      temp_x   <- temp_x+theta[2*l+1]*N2k1(l,sub_t,x_est)+theta[2*l+2]*N2k2(l,sub_t,x_est)
      temp_dx  <- temp_dx+theta[2*l+1]*dN2k1(l,sub_t,x_est)+theta[2*l+2]*dN2k2(l,sub_t,x_est)
      temp_d2x <- temp_d2x+theta[2*l+1]*d2N2k1(l,sub_t,x_est)+theta[2*l+2]*d2N2k2(l,sub_t,x_est)
     }
    temp_x   <- temp_x+theta[2*rowlen-1]*N2n1(rowlen-1,sub_t,x_est)+theta[2*rowlen]*N2n(rowlen,sub_t,x_est)
    temp_dx  <- temp_dx+theta[2*rowlen-1]*dN2n1(rowlen-1,sub_t,x_est)+theta[2*rowlen]*dN2n(rowlen,sub_t,x_est)
    temp_d2x <- temp_d2x+theta[2*rowlen-1]*d2N2n1(rowlen-1,sub_t,x_est)+theta[2*rowlen]*d2N2n(rowlen,sub_t,x_est)
    ft[k]    <- temp_x
    dft[k]   <- temp_dx
    d2ft[k]  <- temp_d2x
  }
  return(list(t = tx,y = ft,v = dft,av = d2ft))
}
