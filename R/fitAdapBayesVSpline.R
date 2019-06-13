#'
#' @title Bayes estimate of adaptive V-Splines
#' @description Fitting adaptive V-splines with piecewise \eqn{\lambda}s.
#' @export
fitAdapBayesVSpline <- function(dat,W=NULL,U=NULL,pa=NULL,a=3,b=2,xout=NULL){
  if(is.null(xout)) x_star <- dat$t
  else x_star <- xout
  
  X<- dat$t
  Y<- dat$y
  V<- dat$v
  rowlen    <- nrow(dat)
  GPmean_x  <- numeric(length(x_star))
  GPmean_v  <- numeric(length(x_star))
  GPsig_x   <- numeric(length(x_star))
  fitmu     <- numeric(length(x_star))
  fitvar    <- numeric(length(x_star))
  
  if(is.null(W)) W <- diag(length(X))
  if(is.null(U)) U <- diag(length(X))
  
  d   <- X[2:rowlen]-X[1:(rowlen-1)]
  dist<- abs(diff(Y))
  dist[which(dist==0)]=1e-4
  
  lambda <- c(Inf,d^a/dist^b,Inf)
  
  if(is.null(pa)){
    pa      <- c(0,0)
    paraGPR <- optim(pa=pa,fn=remlAdaptVSpline,X=X,Y=Y,V=V,lambs=lambda)
    labs <- lambda*exp(paraGPR$par[1])
    gm   <- exp(paraGPR$par[2])
  }else{
    labs <- lambda*pa[1]
    gm   <- pa[2]
  }
  
  coff <- adapCofBayesVSpline(X=X,Y=Y,V=V,lambs=labs,gam=gm)
  
  for(l in 1:length(x_star)){
    gen        <- adapBayesVSpline(X=X,Y=Y,V=V,coff,lambs=labs,x_star[l])
    GPmean_x[l]<- gen$mu
    GPsig_x[l] <- gen$sig/rowlen
    GPmean_v[l]<- gen$mv
  }
  return(list(t=x_star,y=GPmean_x,adpSig=GPsig_x,v=GPmean_v,lambdas=labs,gam=gm))
}