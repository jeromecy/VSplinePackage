#' @title Kernel functions of Bayes estimate of V-Splines
#'
#' @description adaptive Kernel functions \eqn{R_w^1(s,t), \dot{R}_w^1(s,t), 
#' R_w^{'1}(s,t), \dot{R}_w^{'1}(s,t)}.
#'
#' @keywords adaptive Kernel functions
#'
#' @references Z. Cao, D. Bryant, C. Fox, T. Molten and M. Parry. "V-Spline: an Adaptive Smoothing Spline
#' for Trajectory Reconstruction" arXiv preprint arXiv:1803.07184 (2018).
#' @references Z. Cao, D. Bryant, and M. Parry. "V-spline and its Bayes estimate" arXiv (2018).
#' @export
kernelR1adap<- function(x,y,lab,alt){
  v   <- min(x,y)
  alt <- c(0,alt,1)
  lab <- c(lab,Inf)
  l   <- findInterval(v,alt)
  ti1 <- alt[2:l]     ### t[i+1]
  ti  <- alt[1:(l-1)] ### t[i]
  
  p1<- sum((-(x-ti1)*(y-ti1)^2/2+(x-ti)*(y-ti)^2/2+
              (y-ti1)^3/6-(y-ti)^3/6)/lab[1:(l-1)])
  p2<- (-(x-v)*(y-v)^2/2+(x-alt[l])*(y-alt[l])^2/2+
          (y-v)^3/6-(y-alt[l])^3/6)/lab[l]
  
  out <- p1+p2
  
  return(out)
}
#' @rdname kernelR1adap
dotR1adap<- function(x,y,lab,alt){  
  v   <- min(x,y)
  alt <- c(0,alt,1)
  lab <- c(lab,Inf)
  l   <- findInterval(v,alt)
  
  tl  <- alt[l]
  ti1 <- alt[2:l]     ### t[i+1]
  ti  <- alt[1:(l-1)] ### t[i]
  
  p1<- sum((-(y-ti1)^2/2+(y-ti)^2/2)/lab[1:(l-1)])
  p2<- (y*v-v^2/2-y*alt[l]+alt[l]^2/2)/lab[l]
  
  out <- p1+p2
  return(out)
}
#' @rdname kernelR1adap
dR1adap<- function(x,y,lab,alt){
  v   <- min(x,y)
  alt <- c(0,alt,1)
  lab <- c(lab,Inf)
  l   <- findInterval(v,alt)
  
  tl  <- alt[l]
  ti1 <- alt[2:l]     ### t[i+1]
  ti  <- alt[1:(l-1)] ### t[i]
  
  p1<- sum((-(x-ti1)^2/2+(x-ti)^2/2)/lab[1:(l-1)])
  p2<- (x*v-v^2/2-x*alt[l]+alt[l]^2/2)/lab[l]
  
  out <- p1+p2
  return(out)
}
#' @rdname kernelR1adap
ddotR1adap<- function(x,y,lab,alt){
  v   <- min(x,y)
  alt <- c(0,alt,1)
  lab <- c(lab,Inf)
  l   <- findInterval(v,alt)
  
  tl  <- alt[l]
  ti1 <- alt[2:l]     ### t[i+1]
  ti  <- alt[1:(l-1)] ### t[i]
  
  out<- sum((ti1-ti)/lab[1:(l-1)])+(v-alt[l])/lab[l]
  
  return(out)
}
#' @rdname kernelR1adap
ddotdotR1adap<- function(x,y,lab,alt){
  v   <- min(x,y)
  alt <- c(0,alt,1)
  lab <- c(lab,Inf)
  l   <- findInterval(v,alt)
  if(y>=x) return(0)
  else if(y<x)  return(1/lab[l])
}
