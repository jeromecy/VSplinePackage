#'
#' @title Bayes estimate of V-Spline model
#' @description
#' @export
fitBayesVSpline <- function(dat,W=NULL,U=NULL,pa=c(1,1),xout=NULL){
  if(is.null(xout)) x_star <- dat$x
  else x_star <- xout
  
  X<- dat$x
  Y<- dat$y
  V<- dat$v
  GPmean_x  <- numeric(length(x_star))
  GPmean_v  <- numeric(length(x_star))
  GPsig_x   <- numeric(length(x_star))
  fitmu     <- numeric(length(x_star))
  fitvar    <- numeric(length(x_star))
  
  if(is.null(W)) W <- diag(length(X))
  if(is.null(U)) U <- diag(length(X))
  
  paraGPR <- optim(pa=pa,fn=remlScore_VSpline,X=X,Y=Y,V=V,W=W,U=U)
  coff    <- cofBayesVSpline(X,Y,V,W,U,paraGPR$par)
  for(l in 1:length(x_star)){
    gen        <- BayesVSpline(X,Y,V,coff,x_star[l])
    GPmean_x[l]<- gen$mu
    GPsig_x[l] <- gen$sig
    GPmean_v[l]<- gen$mv
  }
  return(list(x=x_star,y=GPmean_x,sig=GPsig_x,v=GPmean_v))
}



