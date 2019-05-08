#'
#' @title Coefficients
#'
#' @description To reconstruct the coefficients.
#' @export
cofBayesVSpline <- function(X,Y,V,W,U,pa){
  lambda = exp(pa[1])
  gamma  = exp(pa[2])
  
  n=length(X)
  S=matrix(0,nrow=n,ncol=2)  
  S[,1]=1
  S[,2]=X
  
  # dis = diff(X)
  # dd  = abs(diff(Y))
  
  Q<- matrix(0,nrow=n,ncol=n)   
  for(i in 1:n)
    for(j in 1:n)
      Q[i,j]=kernelR1(X[j],X[i])
  
  P<- matrix(0,nrow=n,ncol=n)   
  for(i in 1:n)
    for(j in 1:n)
      P[i,j]=dotR1(X[j],X[i])
  
  dS<- matrix(c(0,1),nrow=n,ncol=2,byrow=TRUE)
  
  dQ<- matrix(0,nrow=n,ncol=n)   
  for(i in 1:n)
    for(j in 1:n)
      dQ[i,j]=dR1(X[j],X[i])
  
  dP<- matrix(0,nrow=n,ncol=n)   
  for(i in 1:n)
    for(j in 1:n)
      dP[i,j]=ddotR1(X[j],X[i])
  
  TT<- rbind(S,dS)
  
  var_y <- Q+n*lambda*W
  cov_yv<- P 
  cov_vy<- dQ   
  var_v <- dP+n*lambda/gamma*U
  
  M <- rbind(cbind(var_y,cov_yv),cbind(cov_vy,var_v))
  
  #inM=solve(M)
  R<- chol(M)
  inM<- solve(R)%*%solve(t(R))
  
  # ob<- c(Y,V)
  
  inW<- solve(t(TT)%*%inM%*%TT)
  pbc<- inM-inM%*%TT%*%solve(t(TT)%*%inM%*%TT)%*%t(TT)%*%inM
  pd <- inW %*%t(TT)%*%inM
  pe <- inM%*%TT%*%inW

  return(list(pbc=pbc,pd=pd,pe=pe,inW=inW))
}