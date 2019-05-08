#' @title Kernel functions of Bayes estimate of V-Splines
#'
#' @description Kernel functions \eqn{R_1(s,t), \dot{R}_1(s,t), 
#' R'_1(s,t), \dot{R}'_1(s,t)}.
#'
#' @keywords Kernel functions
#'
#' @references Z. Cao, D. Bryant, C. Fox, T. Molten and M. Parry. "V-Spline: an Adaptive Smoothing Spline
#' for Trajectory Reconstruction" arXiv preprint arXiv:1803.07184 (2018).
#' @references Z. Cao, D. Bryant, and M. Parry. "V-spline and its Bayes estimate" arXiv (2018).
#' @export
kernelR1<- function(x,y){
  return(abs(x-y)*min(x,y)^2/2+min(x,y)^3/3)  
}
#' @rdname kernelR1
dotR1<- function(x,y){  
  return(y*min(x,y)-min(x,y)^2/2)  
}
#' @rdname kernelR1
dR1<- function(x,y){
  return(x*min(x,y)-min(x,y)^2/2)
}
#' @rdname kernelR1
ddotR1<- function(x,y){
  return(min(x,y))
}
#' @rdname kernelR1
ddotdotR1<- function(x,y){
  if(y>=x) return(1)
  else if(y<x)  return(0)
}
#' @rdname kernelR1
dotdotR1<- function(x,y){
  if(y>=x) return(y-x)
  else if(y<x)  return(0)
}