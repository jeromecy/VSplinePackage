#'
#' @title Bayes estimate for adaptive V-splines
#'
#' @description Fitting the Bayes estimate of adaptive V-splines
#' @export
adapBayesVSpline<- function(X,Y,V,coff,w,gamma,est){
  alt <- X
  ob<- c(Y,V)
  
  bc<- coff$pbc%*%ob
  d <- coff$pd%*%ob
  
  n<- length(X)  
  
  phi<- matrix(c(1,est),nrow=2,ncol=1)
  
  xi<- matrix(0,nrow=n,ncol=1)
  for(i in 1:n)  xi[i,1]=kernelR1adap(X[i],est,w,alt)
  
  psi<- matrix(0,nrow=n,ncol=1)   
  for(i in 1:n)
    psi[i,1]=dotR1adap(X[i],est,w,alt)
  
  newPa <- c(xi,psi)
  
  mu <- t(phi)%*%d+t(newPa)%*%bc
  sig<- kernelR1adap(est,est,w,alt) + t(phi)%*%coff$inW%*%phi- t(newPa)%*%coff$pbc%*%newPa -
    t(phi)%*%coff$pd%*%newPa - t(newPa)%*%coff$pe%*%phi
  
  dphi=matrix(c(0,1),nrow=2,ncol=1)
  
  dxi<- matrix(0,nrow=n,ncol=1)
  for(i in 1:n) dxi[i,1]=dR1adap(X[i],est,w,alt)
  
  dpsi<- matrix(0,nrow=n,ncol=1)   
  for(i in 1:n) dpsi[i,1]=ddotR1adap(X[i],est,w,alt)
  
  dnewPa <- c(dxi,dpsi)
  mv <- t(dphi)%*%d+t(dnewPa)%*%bc
  
  return(list(mu=mu,sig=sig,mv = mv))
}

