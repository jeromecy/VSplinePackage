#'
#' @title REML score of GPR V-splines
#'
#' @description Parameter estimation by REML
#' @export
remlAdaptVSpline<- function(X,Y,V,w,pa){
  #print(pa)
  w     <- w*exp(pa[1])
  gamma <- exp(pa[2])
  
  ob<- c(Y,V)

  n   = length(X)
  alt = X
  
  S     = matrix(0,nrow=n,ncol=2)  
  S[,1] = 1
  S[,2] = X
  
  Q<- matrix(0,nrow=n,ncol=n)   
  for(i in 1:n)
    for(j in 1:n)
      Q[i,j]=kernelR1adap(X[j],X[i],w,alt)
  
  P<- matrix(0,nrow=n,ncol=n)   
  for(i in 1:n)
    for(j in 1:n)
      P[i,j]=dotR1adap(X[j],X[i],w,alt)
  
  dS<- matrix(c(0,1),nrow=n,ncol=2,byrow=TRUE)
  
  dQ<- matrix(0,nrow=n,ncol=n)   
  for(i in 1:n)
    for(j in 1:n)
      dQ[i,j]=dR1adap(X[j],X[i],w,alt)
  
  dP<- matrix(0,nrow=n,ncol=n)   
  for(i in 1:n)
    for(j in 1:n)
      dP[i,j]=ddotR1adap(X[j],X[i],w,alt)
  
  TT<- rbind(S,dS)
  
  var_y <- Q+n*diag(n)
  cov_yv<- P 
  cov_vy<- dQ   
  var_v <- dP+n*diag(n)/gamma
  
  M <- rbind(cbind(var_y,cov_yv),cbind(cov_vy,var_v))
  
  R  <- chol(M)
  inM<- solve(R)%*%solve(t(R))
  
  inW<- solve(t(TT)%*%inM%*%TT)
  pbc<- inM-inM%*%TT%*%solve(t(TT)%*%inM%*%TT)%*%t(TT)%*%inM
  
  matA <- diag(2*n)-diag(rep(c(n,n/gamma),each=n))%*%pbc
  
  reml <- t(ob)%*%t(diag(2*n)-t(matA))%*%(diag(2*n)-t(matA))%*%ob/
    sum(diag(diag(2*n)-matA))
  
  return(reml)
}