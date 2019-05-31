#' @title Kernel functions of Bayes estimate of V-Splines
#'
#' @description Kernel functions \eqn{R_1(s,t), \dot{R}_1(s,t), 
#' R'_1(s,t), \dot{R}'_1(s,t)}.
#'
#' @keywords Kernel functions
#' 
#' @seealso \code{\link{adapKernel}} for adaptive kernel functions.
#'
#' @references Z. Cao, D. Bryant, C. Fox, T. Molten and M. Parry. "V-Spline: an Adaptive Smoothing Spline
#' for Trajectory Reconstruction" arXiv preprint arXiv:1803.07184 (2018).
#' @references Z. Cao, D. Bryant, and M. Parry. "V-spline and its Bayes estimate" arXiv (2018).
#' @export
kernelR1<- function(x1,y1){
  return(abs(x1-y1)*min(x1,y1)^2/2+min(x1,y1)^3/3)  
}
#' @rdname kernelR1
dotR1<- function(x1,y1){  
  return(y1*min(x1,y1)-min(x1,y1)^2/2)  
}
#' @rdname kernelR1
dR1<- function(x1,y1){
  return(x1*min(x1,y1)-min(x1,y1)^2/2)
}
#' @rdname kernelR1
ddotR1<- function(x1,y1){
  return(min(x1,y1))
}
#' @rdname kernelR1
ddotdotR1<- function(x1,y1){
  if(y1>=x1) return(0) # was 1
  else if(y1<x1)  return(1) # was 0
}
#' @rdname kernelR1
dotdotR1<- function(x1,y1){
  if(y1>=x1) return(y1-x1)
  else if(y1<x1)  return(0)
}