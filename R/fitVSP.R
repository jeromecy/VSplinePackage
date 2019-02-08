#'
#' @title Fit V-Spline model
#'
#' @description
#' @references Z. Cao, D. Bryant, C. Fox, T. Molten and M. Parry. "V-Spline: an Adaptive Smoothing Spline
#' for Trajectory Reconstruction" arXiv preprint arXiv:1803.07184 (2018).
#' @export
fitVSP <- function(dat,pa=c(1,1)){

  mats  <- iniVSPmat(dat,a=3,b=2)
  oppa  <- optim(fn = loocvVSP,dat=dat, pa = pa, matrix_list = mats)
  theta <- getTheta(dat,oppa$par,mats)
  rct   <- rcstVSP(dat,theta)

  return(rct)
}
