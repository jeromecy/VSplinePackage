#'
#' @title Fit V-Spline model
#'
#' @description
#' @references Z. Cao, D. Bryant, C. Fox, T. Molten and M. Parry. "V-Spline: an Adaptive Smoothing Spline
#' for Trajectory Reconstruction" arXiv preprint arXiv:1803.07184 (2018).
#' @seealso \code{\link{getTheta}}, \code{\link{loocvVSP}}
#' @examples
#' n <-100
#' s <- seq(1,2*pi,length=n)
#' y <- sin(s) + rnorm(n)
#' v <- cos(s) + rnorm(n)
#' simuData <- data.frame(t=s,x=y,mx=v)
#' fitted   <- fitVSP(simuData)
#' plot(s,y)
#' points(fitted$t,fitted$x,col="red",type="l")
#' @export
fitVSP <- function(dat,pa=c(1,1)){

  mats  <- iniVSPmat(dat,a=3,b=2)
  oppa  <- optim(fn = loocvVSP,dat=dat, pa = pa, matrix_list = mats)
  theta <- getTheta(dat,oppa$par,mats)
  rct   <- rcstVSP(dat,theta)

  return(rct)
}
