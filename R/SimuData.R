#' @title Generate simulation data
#'
#' @description This function takes *noise-free* velocity data
#' computes the implied position data and then adds noise.
#' It requires function `forward1`.
#' @keywords V-Spline basis functions#'
#' @param v noise-free velocity data
#' @param rsnr signal to noise ratio
#' @example 
#' library(waveband)
#' library(VSPline)
#' set.seed(2016)
#' N = 100
#' velocity <- test.data(type = "blocks", n = N, signal = 1, rsnr = 7, plotfn = TRUE)
#' set.seed(2016)
#' position <- SimuData(velocity$y,7)
#' simulatedata<- data.frame(t=velocity$x,y=position$xnoise,v=velocity$ynoise,boom=0)
#' @export
SimuData = function(v,rsnr) {
  n  = length(v)
  dt = 1/n
  t  = (0:(n-1))*dt
  b  = list(x=numeric(n),xnoise=numeric(n))
  x  = numeric(n)
  x[1] = 0
  for (i in 1:(n-1)) {
    x[i+1] = forward1(t[i],t[i+1],v[i],v[i+1],x[i])
  }
  sig = sd(x)
  b$x = x
  b$xnoise = x + rnorm(n,0,sig/rsnr)
  return(b)
}
