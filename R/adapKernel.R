#' @title Kernel functions of adaptive V-splines
#'
#' @description adaptive Kernel functions \eqn{R_{\lambda 1}(s,t), \dot{R}_{\lambda 1}(s,t), 
#' R_{\lambda 1}'(s,t), \dot{R}_{\lambda 1}'(s,t)}.
#'
#' @keywords adaptive Kernel functions
#' 
#' @seealso \code{\link{kernels}} for generic kernels.
#'
#' @references Z. Cao, D. Bryant, C. Fox, T. Molten and M. Parry. "V-Spline: an Adaptive Smoothing Spline
#' for Trajectory Reconstruction" arXiv preprint arXiv:1803.07184 (2018).
#' @references Z. Cao, D. Bryant, and M. Parry. "V-spline and its Bayes estimate" arXiv (2018).
#' @export
kernelR1adap<- function(x1,y1,lambs,alt){
  miv   <- min(x1,y1)
  alt <- c(0,alt,1)
  lambs <- c(lambs,Inf)
  l   <- findInterval(miv,alt)
  ti1 <- alt[2:l]     ### t[i+1]
  ti  <- alt[1:(l-1)] ### t[i]
  
  p1<- sum((-(x1-ti1)*(y1-ti1)^2/2+(x1-ti)*(y1-ti)^2/2+
              (y1-ti1)^3/6-(y1-ti)^3/6)/lambs[1:(l-1)])
  p2<- (-(x1-miv)*(y1-miv)^2/2+(x1-alt[l])*(y1-alt[l])^2/2+
          (y1-miv)^3/6-(y1-alt[l])^3/6)/lambs[l]
  
  out <- p1+p2
  
  return(out)
}
#' @rdname kernelR1adap
dotR1adap<- function(x1,y1,lambs,alt){  
  miv   <- min(x1,y1)
  alt <- c(0,alt,1)
  lambs <- c(lambs,Inf)
  l   <- findInterval(miv,alt)
  
  tl  <- alt[l]
  ti1 <- alt[2:l]     ### t[i+1]
  ti  <- alt[1:(l-1)] ### t[i]
  
  p1<- sum((-(y1-ti1)^2/2+(y1-ti)^2/2)/lambs[1:(l-1)])
  p2<- (y1*miv-miv^2/2-y1*alt[l]+alt[l]^2/2)/lambs[l]
  
  out <- p1+p2
  return(out)
}
#' @rdname kernelR1adap
dR1adap<- function(x1,y1,lambs,alt){
  miv   <- min(x1,y1)
  alt <- c(0,alt,1)
  lambs <- c(lambs,Inf)
  l   <- findInterval(miv,alt)
  
  tl  <- alt[l]
  ti1 <- alt[2:l]     ### t[i+1]
  ti  <- alt[1:(l-1)] ### t[i]
  
  p1<- sum((-(x1-ti1)^2/2+(x1-ti)^2/2)/lambs[1:(l-1)])
  p2<- (x1*miv-miv^2/2-x1*alt[l]+alt[l]^2/2)/lambs[l]
  
  out <- p1+p2
  return(out)
}
#' @rdname kernelR1adap
ddotR1adap<- function(x1,y1,lambs,alt){
  miv   <- min(x1,y1)
  alt <- c(0,alt,1)
  lambs <- c(lambs,Inf)
  l   <- findInterval(miv,alt)
  
  tl  <- alt[l]
  ti1 <- alt[2:l]     ### t[i+1]
  ti  <- alt[1:(l-1)] ### t[i]
  
  out<- sum((ti1-ti)/lambs[1:(l-1)])+(miv-alt[l])/lambs[l]
  
  return(out)
}
#' @rdname kernelR1adap
ddotdotR1adap<- function(x1,y1,lambs,alt){
  miv   <- min(x1,y1)
  alt <- c(0,alt,1)
  lambs <- c(lambs,Inf)
  l   <- findInterval(miv,alt)
  if(y1>=x1) return(0)
  else if(y1<x1)  return(1/lambs[l])
}
