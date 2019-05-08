#'
#' @title REML score of GPR V-splines
#'
#' @description Parameter estimation by REML
#' @export
remlScoreVSpline<- function(X,Y,V,W,U,pa){
  #print(pa)
  lambda = exp(pa[1])
  gamma  = exp(pa[2])
  
  ob<- c(Y,V)
  
  n=length(X)
  S=matrix(0,nrow=n,ncol=2)  
  S[,1]=1
  S[,2]=X
  
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
  R <- chol(M)
  inM<- solve(R)%*%solve(t(R))
  
  inW<- solve(t(TT)%*%inM%*%TT)
  pbc<- inM-inM%*%TT%*%solve(t(TT)%*%inM%*%TT)%*%t(TT)%*%inM
  # pd <- inW %*%t(TT)%*%inM
  # pe <- inM%*%TT%*%inW
  # 
  # bc<- pbc%*%ob
  # d <- pd%*%ob
  
  # Q1 <- qr.Q(qr(TT),complete = TRUE)[,1:2]
  # Q2 <- qr.Q(qr(TT),complete = TRUE)[,3:94]
  
  # A <- diag(n) - n*lambda*Q2%*%solve(t(Q2)%*%(Q+diag(n*lambda,n))%*%Q2)%*%t(Q2)
  # B <- diag(n) - n*lambda/gamma*dQ2%*%solve(t(dQ2)%*%(dP+diag(n*lambda/gamma,n))%*%dQ2)%*%t(dQ2)
  # matA <- rbind(cbind(A,matrix(0,n,n)),cbind(matrix(0,n,n),B))
  
  matA <- diag(2*n)-diag(rep(c(lambda*n,lambda*n/gamma),each=n))%*%pbc
  
  reml <- t(ob)%*%t(diag(2*n)-t(matA))%*%(diag(2*n)-t(matA))%*%ob/
    sum(diag(diag(2*n)-matA))
  
  return(reml)
}