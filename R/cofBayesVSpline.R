#'
#' @title Coefficients
#'
#' @description To reconstruct the coefficients.
#' @export
cofBayesVSpline <- function(X,Y,V,W=NULL,U=NULL,lamb,gam){
  # lamb = exp(pa[1])
  # gam  = exp(pa[2])
  
  rowlen=length(X)
  if(is.null(W)) W = diag(rowlen)
  if(is.null(U)) U = diag(rowlen)
  
  S=matrix(0,nrow=rowlen,ncol=2)  
  S[,1]=1
  S[,2]=X
  
  # dis = diff(X)
  # dd  = abs(diff(Y))
  
  Q<- matrix(0,nrow=rowlen,ncol=rowlen)   
  for(i in 1:rowlen)
    for(j in 1:rowlen)
      Q[i,j]=kernelR1(X[j],X[i])
  
  P<- matrix(0,nrow=rowlen,ncol=rowlen)   
  for(i in 1:rowlen)
    for(j in 1:rowlen)
      P[i,j]=dotR1(X[j],X[i])
  
  dS<- matrix(c(0,1),nrow=rowlen,ncol=2,byrow=TRUE)
  
  dQ<- matrix(0,nrow=rowlen,ncol=rowlen)   
  for(i in 1:rowlen)
    for(j in 1:rowlen)
      dQ[i,j]=dR1(X[j],X[i])
  
  dP<- matrix(0,nrow=rowlen,ncol=rowlen)   
  for(i in 1:rowlen)
    for(j in 1:rowlen)
      dP[i,j]=ddotR1(X[j],X[i])
  
  TT<- rbind(S,dS)
  
  var_y <- Q+rowlen*lamb*W
  cov_yv<- P 
  cov_vy<- dQ   
  var_v <- dP+rowlen*lamb/gam*U
  
  M <- rbind(cbind(var_y,cov_yv),cbind(cov_vy,var_v))
  
  #inM=solve(M)
  R  <- chol(M)
  inM<- solve(R)%*%solve(t(R))
  
  # ob<- c(Y,V)
  
  inW<- solve(t(TT)%*%inM%*%TT)
  pbc<- inM-inM%*%TT%*%solve(t(TT)%*%inM%*%TT)%*%t(TT)%*%inM
  pd <- inW %*%t(TT)%*%inM
  pe <- inM%*%TT%*%inW

  return(list(pbc=pbc,pd=pd,pe=pe,inW=inW))
}