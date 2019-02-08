#'
#' @title Reconstruction points
#'
#' @description To reconstruct the V-Spline with \eqn{\theta}, which is gained from \code{getTheta()}.
#' @references Z. Cao, D. Bryant, C. Fox, T. Molten and M. Parry. "V-Spline: an Adaptive Smoothing Spline
#' for Trajectory Reconstruction" arXiv preprint arXiv:1803.07184 (2018).
#' @export
rcstVSP <- function(data_frame,theta){
  sub_t <- data_frame$t
  n     <- length(sub_t)

  tx   <- seq(sub_t[1],sub_t[n],length = (n*10))
  fx   <- numeric(length(tx))
  dfx  <- numeric(length(tx))
  d2fx <- numeric(length(tx))

  for(k in 1:length(tx)){
    x_est    <- tx[k]
    temp_x   <- 0
    temp_dx  <- 0
    temp_d2x <- 0

    temp_x   <- theta[1]*N1(1,sub_t,x_est)+theta[2]*N2(1,sub_t,x_est)
    temp_dx  <- theta[1]*dN1(1,sub_t,x_est)+theta[2]*dN2(1,sub_t,x_est)
    temp_d2x <- theta[1]*d2N1(1,sub_t,x_est)+theta[2]*d2N2(1,sub_t,x_est)

    for(l in 1:(length(sub_t)-2)){
      temp_x   <- temp_x+theta[2*l+1]*N2k1(l,sub_t,x_est)+theta[2*l+2]*N2k2(l,sub_t,x_est)
      temp_dx  <- temp_dx+theta[2*l+1]*dN2k1(l,sub_t,x_est)+theta[2*l+2]*dN2k2(l,sub_t,x_est)
      temp_d2x <- temp_d2x+theta[2*l+1]*d2N2k1(l,sub_t,x_est)+theta[2*l+2]*d2N2k2(l,sub_t,x_est)
     }
    temp_x   <- temp_x+theta[2*n-1]*N2n1(n-1,sub_t,x_est)+theta[2*n]*N2n(n,sub_t,x_est)
    temp_dx  <- temp_dx+theta[2*n-1]*dN2n1(n-1,sub_t,x_est)+theta[2*n]*dN2n(n,sub_t,x_est)
    temp_d2x <- temp_d2x+theta[2*n-1]*d2N2n1(n-1,sub_t,x_est)+theta[2*n]*d2N2n(n,sub_t,x_est)
    fx[k]    <- temp_x
    dfx[k]   <- temp_dx
    d2fx[k]  <- temp_d2x
  }
  return(list(t = tx,x = fx,vx = dfx,ax = d2fx))
}
