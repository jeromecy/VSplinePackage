#'
#' @title Bayes estimate for adaptive V-splines
#'
#' @description Fitting the Bayes estimate of adaptive V-splines
#' @export
adapBayesVSpline<- function(X,Y,V,coff,lambs,est){
  alt <- X
  ob  <- c(Y,V)
  
  bc<- coff$pbc%*%ob
  d <- coff$pd%*%ob
  rowlen <- length(X)  
  
  phi<- matrix(c(1,est),nrow=2,ncol=1)
  
  xi <- matrix(0,nrow=rowlen,ncol=1)
  for(i in 1:rowlen)  xi[i,1]<- kernelR1adap(X[i],est,lambs,alt)
  
  psi<- matrix(0,nrow=rowlen,ncol=1)   
  for(i in 1:rowlen)  psi[i,1]<- dotR1adap(X[i],est,lambs,alt)
  
  newPa <- c(xi,psi)
  
  mu <- t(phi)%*%d+t(newPa)%*%bc
  sig<- kernelR1adap(est,est,lambs,alt) + t(phi)%*%coff$inW%*%phi- t(newPa)%*%coff$pbc%*%newPa -
    t(phi)%*%coff$pd%*%newPa - t(newPa)%*%coff$pe%*%phi
  
  dphi<-matrix(c(0,1),nrow=2,ncol=1)
  
  dxi <- matrix(0,nrow=rowlen,ncol=1)
  for(i in 1:rowlen) dxi[i,1] <- dR1adap(X[i],est,lambs,alt)
  
  dpsi<- matrix(0,nrow=rowlen,ncol=1)   
  for(i in 1:rowlen) dpsi[i,1]<- ddotR1adap(X[i],est,lambs,alt)
  
  dnewPa <- c(dxi,dpsi)
  mv     <- t(dphi)%*%d+t(dnewPa)%*%bc
  
  return(list(mu=mu,sig=sig,mv = mv))
}

