#'
#' @title REML score for adaptive V-splines
#'
#' @description Parameter estimation by REML
#' @export
remlAdaptVSpline<- function(X,Y,V,lambs,pa){
  # print(pa)
  lambs  <- lambs*exp(pa[1])
  gam    <- exp(pa[2])
  ob     <- c(Y,V)
  rowlen <- length(X)
  alt    <- X
  
  S     <- matrix(0,nrow=rowlen,ncol=2)  
  S[,1] <- 1
  S[,2] <- X
  
  Q<- matrix(0,nrow=rowlen,ncol=rowlen)   
  for(i in 1:rowlen)
    for(j in 1:rowlen)
      Q[i,j]<-kernelR1adap(X[j],X[i],lambs,alt)
  
  P<- matrix(0,nrow=rowlen,ncol=rowlen)   
  for(i in 1:rowlen)
    for(j in 1:rowlen)
      P[i,j]<-dotR1adap(X[j],X[i],lambs,alt)
  
  dS<- matrix(c(0,1),nrow=rowlen,ncol=2,byrow=TRUE)
  
  dQ<- matrix(0,nrow=rowlen,ncol=rowlen)   
  for(i in 1:rowlen)
    for(j in 1:rowlen)
      dQ[i,j]<-dR1adap(X[j],X[i],lambs,alt)
  
  dP<- matrix(0,nrow=rowlen,ncol=rowlen)   
  for(i in 1:rowlen)
    for(j in 1:rowlen)
      dP[i,j]<-ddotR1adap(X[j],X[i],lambs,alt)
  
  TT<- rbind(S,dS)
  
  var_y <- Q+rowlen*diag(rowlen)
  cov_yv<- P 
  cov_vy<- dQ   
  var_v <- dP+rowlen*diag(rowlen)/gam
  
  M  <- rbind(cbind(var_y,cov_yv),cbind(cov_vy,var_v))
  
  R  <- chol(M)
  inM<- solve(R)%*%solve(t(R))
  
  inW<- solve(t(TT)%*%inM%*%TT)
  pbc<- inM-inM%*%TT%*%solve(t(TT)%*%inM%*%TT)%*%t(TT)%*%inM
  
  matA <- diag(2*rowlen)-diag(rep(c(rowlen,rowlen/gam),each=rowlen))%*%pbc
  
  reml <- t(ob)%*%t(diag(2*rowlen)-t(matA))%*%(diag(2*rowlen)-t(matA))%*%ob/
    sum(diag(diag(2*rowlen)-matA))
  
  return(reml)
}