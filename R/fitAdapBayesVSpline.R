#'
#' @title Bayes estimate of adaptive V-Splines
#' @description
#' @export
fitAdapBayesVSpline <- function(dat,W=NULL,U=NULL,
                                pa=c(1,1),xout=NULL){
  if(is.null(xout)) x_star <- dat$x
  else x_star <- xout
  
  X<- dat$x
  Y<- dat$y
  V<- dat$v
  n<- length(X)
  GPmean_x  <- numeric(length(x_star))
  GPmean_v  <- numeric(length(x_star))
  GPsig_x   <- numeric(length(x_star))
  fitmu     <- numeric(length(x_star))
  fitvar    <- numeric(length(x_star))
  
  if(is.null(W)) W <- diag(length(X))
  if(is.null(U)) U <- diag(length(X))
  
  d   <- X[2:n]-X[1:(n-1)]
  dist<- abs(diff(Y))
  dist[which(dist==0)]=1e-4
  
  lambda <- c(Inf,d^3/dist^2,Inf)
  
  paraGPR <- optim(pa=pa,fn=remlAdaptVSpline,X=X,Y=Y,V=V,w=lambda)
  
  w  <- lambda*exp(paraGPR$par[1])
  gm <- exp(paraGPR$par[2])
  coff <- adapCofBayesVSpline(X=s,Y=r,V=v1,w=w,gamma=gm)
  
  for(l in 1:length(x_star)){
    gen        <- adapBayesVSpline(X=s,Y=r,V=v1,coff,
                                   w=w,gamma=gm,x_star[l])
    GPmean_x[l]<- gen$mu
    GPsig_x[l] <- gen$sig
    GPmean_v[l]<- gen$mv
  }
  return(list(x=x_star,y=GPmean_x,sig=GPsig_x,v=GPmean_v))
}