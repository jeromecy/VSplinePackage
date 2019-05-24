#'
#' @title Coefficients for adaptive V-splines
#'
#' @description To reconstruct the coefficients for adaptive V-splines.
#' @export
adapCofBayesVSpline <- function(X,Y,V,lab,gamma){
  rowlen   = length(X)
  alt = X
  
  S     = matrix(0,nrow=rowlen,ncol=2)  
  S[,1] = 1
  S[,2] = X
  
  Q<- matrix(0,nrow=rowlen,ncol=rowlen)   
  for(i in 1:rowlen)
    for(j in 1:rowlen)
      Q[i,j]=kernelR1adap(X[j],X[i],lab,alt)
  
  P<- matrix(0,nrow=rowlen,ncol=rowlen)   
  for(i in 1:rowlen)
    for(j in 1:rowlen)
      P[i,j]=dotR1adap(X[j],X[i],lab,alt)
  
  dS<- matrix(c(0,1),nrow=rowlen,ncol=2,byrow=TRUE)
  
  dQ<- matrix(0,nrow=rowlen,ncol=rowlen)   
  for(i in 1:rowlen)
    for(j in 1:rowlen)
      dQ[i,j]=dR1adap(X[j],X[i],lab,alt)
  
  dP<- matrix(0,nrow=rowlen,ncol=rowlen)   
  for(i in 1:rowlen)
    for(j in 1:rowlen)
      dP[i,j]=ddotR1adap(X[j],X[i],lab,alt)
  
  TT<- rbind(S,dS)
  
  var_y <- Q+rowlen*diag(rowlen)
  cov_yv<- P 
  cov_vy<- dQ   
  var_v <- dP+rowlen*diag(rowlen)/gamma
  
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