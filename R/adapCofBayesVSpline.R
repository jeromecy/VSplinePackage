#'
#' @title Coefficients for adaptive V-splines
#'
#' @description To reconstruct the coefficients for adaptive V-splines.
#' @export
adapCofBayesVSpline <- function(X,Y,V,lambs,gam){
  rowlen = length(X)
  alt    = X
  
  adapS     = matrix(0,nrow=rowlen,ncol=2)  
  adapS[,1] = 1
  adapS[,2] = X
  
  adapQ<- matrix(0,nrow=rowlen,ncol=rowlen)   
  for(i in 1:rowlen)
    for(j in 1:rowlen)
      adapQ[i,j]=kernelR1adap(X[j],X[i],lambs,alt)
  
  adapP<- matrix(0,nrow=rowlen,ncol=rowlen)   
  for(i in 1:rowlen)
    for(j in 1:rowlen)
      adapP[i,j]=dotR1adap(X[j],X[i],lambs,alt)
  
  adapdS<- matrix(c(0,1),nrow=rowlen,ncol=2,byrow=TRUE)
  
  adapdQ<- matrix(0,nrow=rowlen,ncol=rowlen)   
  for(i in 1:rowlen)
    for(j in 1:rowlen)
      adapdQ[i,j]=dR1adap(X[j],X[i],lambs,alt)
  
  adapdP<- matrix(0,nrow=rowlen,ncol=rowlen)   
  for(i in 1:rowlen)
    for(j in 1:rowlen)
      adapdP[i,j]=ddotR1adap(X[j],X[i],lambs,alt)
  
  TT<- rbind(adapS,adapdS)
  
  var_y <- adapQ+rowlen*diag(rowlen)
  cov_yv<- adapP 
  cov_vy<- adapdQ   
  var_v <- adapdP+rowlen*diag(rowlen)/gam
  
  Mlab <- rbind(cbind(var_y,cov_yv),cbind(cov_vy,var_v))
  
  # inM=solve(Mlab)
  Rlab <- chol(Mlab)
  inM  <- solve(Rlab)%*%solve(t(Rlab))
  
  # ob<- c(Y,V)
  
  inW<- solve(t(TT)%*%inM%*%TT)
  pbc<- inM-inM%*%TT%*%solve(t(TT)%*%inM%*%TT)%*%t(TT)%*%inM
  pd <- inW %*%t(TT)%*%inM
  pe <- inM%*%TT%*%inW
  
  return(list(pbc=pbc,pd=pd,pe=pe,inW=inW))
}