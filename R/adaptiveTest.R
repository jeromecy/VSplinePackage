



adapCofBayesVSpline <- function(X,Y,V,W,U,pa){
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
      Q[i,j]=kernelR1adap(X[j],X[i])
  
  P<- matrix(0,nrow=n,ncol=n)   
  for(i in 1:n)
    for(j in 1:n)
      P[i,j]=dotR1adap(X[j],X[i])
  
  dS<- matrix(c(0,1),nrow=n,ncol=2,byrow=TRUE)
  
  dQ<- matrix(0,nrow=n,ncol=n)   
  for(i in 1:n)
    for(j in 1:n)
      dQ[i,j]=dR1adap(X[j],X[i])
  
  dP<- matrix(0,nrow=n,ncol=n)   
  for(i in 1:n)
    for(j in 1:n)
      dP[i,j]=ddotR1adap(X[j],X[i])
  
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

adapBayesVSpline<- function(X,Y,V,coff,est){
  
  ob<- c(Y,V)
  
  bc<- coff$pbc%*%ob
  d <- coff$pd%*%ob
  
  n<- length(X)  
  
  phi=matrix(c(1,est),nrow=2,ncol=1)
  
  xi<- matrix(0,nrow=n,ncol=1)
  for(i in 1:n)  xi[i,1]=kernelR1adap(X[i],est)
  
  psi<- matrix(0,nrow=n,ncol=1)   
  for(i in 1:n)
    psi[i,1]=dotR1adap(X[i],est)
  
  newPa <- c(xi,psi)
  
  mu <- t(phi)%*%d+t(newPa)%*%bc
  sig<- R1adap(est,est) + t(phi)%*%coff$inW%*%phi- t(newPa)%*%coff$pbc%*%newPa -
    t(phi)%*%coff$pd%*%newPa - t(newPa)%*%coff$pe%*%phi
  
  dphi=matrix(c(0,1),nrow=2,ncol=1)
  
  dxi<- matrix(0,nrow=n,ncol=1)
  for(i in 1:n) dxi[i,1]=dR1adap(X[i],est)
  
  dpsi<- matrix(0,nrow=n,ncol=1)   
  for(i in 1:n) dpsi[i,1]=ddotR1adap(X[i],est)
  
  dnewPa <- c(dxi,dpsi)
  mv <- t(dphi)%*%d+t(dnewPa)%*%bc
  
  return(list(mu=mu,sig=sig,mv = mv))
}



s<- rep(1:6)/10
r<- c(0.6,0.7,1.3,0.8,0.9,0.8)
plot(s,r)

n <- length(s)
w<-rep(0.1,n)

adapKernel(0.3,0.3,w,s)

Q<- matrix(0,nrow=n,ncol=n)   
for(i in 1:n)
  for(j in 1:n)
    Q[i,j]=adapKernel(s[j],s[i],w,s)

S=matrix(0,nrow=n,ncol=2)  
S[,1]=1
S[,2]=s

M <- Q+n*1/w
inM <- solve(M)
inW<- solve(t(S)%*%inM%*%S)

pbc<- inM-inM%*%S%*%solve(t(S)%*%inM%*%S)%*%t(S)%*%inM
pd <- inW %*%t(S)%*%inM

bc<- pbc%*%r
d <- pd%*%r

s_star <- seq(0.1,0.6,length=10)
mu <- numeric(length(s_star))

for(i in 1:length(mu)){
  phi=matrix(c(1,s_star[i]),nrow=2,ncol=1)
  xi<- matrix(0,nrow=n,ncol=1)
  for(k in 1:n)  xi[k,1]=adapKernel(s[k],s_star[i],w,s)
  mu[i] <- t(phi)%*%d+t(xi)%*%bc
}

points(s_star,mu,pch=20,type="b")








adapKernel <- function(x,y,w,alt){
  
  if(y>=x){
    l   <- findInterval(x,alt)
    ti1 <- alt[2:l]     ### t[i+1]
    ti  <- alt[1:(l-1)] ### t[i]
    
    p1<- sum((-(x-ti1)^2/2*(y-ti1)+(x-ti)^2/2*(y-ti)+
                (x-ti1)^3/6-(x-ti)^3/6)/w[1:(l-1)])
    p2<- ((x-alt[l])^2/2*(y-alt[l])-(x-alt[l])^3/6)/w[l]
    
    adpR1 <- p1+p2
  }else if(y<x){
    
    l   <- findInterval(y,alt)
    ti1 <- alt[2:l]     ### t[i+1]
    ti  <- alt[1:(l-1)] ### t[i]
    
    p1<- sum((-(x-ti1)^2/2*(y-ti1)+(x-ti)^2/2*(y-ti)+
                (x-ti1)^3/6-(x-ti)^3/6)/w[1:(l-1)])
    p2<- ((x-alt[l])^2/2*(y-alt[l])+(x-y)^3/6-(x-alt[l])^3/6)/w[l]
    
    adpR1 <- p1+p2
  }
  
  return(adpR1)
  
}

kernelR1adap<- function(x,y,w,alt){
  if(y>=x){
    l   <- findInterval(x,alt)
    ti1 <- alt[2:l]     ### t[i+1]
    ti  <- alt[1:(l-1)] ### t[i]
    
    p1<- sum((-(x-ti1)^2/2*(y-ti1)+(x-ti)^2/2*(y-ti)+
                (x-ti1)^3/6-(x-ti)^3/6)/w[1:(l-1)])
    p2<- ((x-alt[l])^2/2*(y-alt[l])-(x-alt[l])^3/6)/w[l]
    
    out <- p1+p2
  }else if(y<x){
    l   <- findInterval(y,alt)
    ti1 <- alt[2:l]     ### t[i+1]
    ti  <- alt[1:(l-1)] ### t[i]
    
    p1<- sum((-(x-ti1)^2/2*(y-ti1)+(x-ti)^2/2*(y-ti)+
                (x-ti1)^3/6-(x-ti)^3/6)/w[1:(l-1)])
    p2<- ((x-alt[l])^2/2*(y-alt[l])+(x-y)^3/6-(x-alt[l])^3/6)/w[l]
    
    out <- p1+p2
  }
  return(out)
}
dotR1adap<- function(x,y,w,alt){  
  l <- findInterval(min(x,y),alt)
  v <- min(x,alt[l],y)
  
  ti1 <- alt[2:l]     ### t[i+1]
  ti  <- alt[1:(l-1)] ### t[i]
  
  p1<- sum((-(y-ti1)^2/2+(y-ti)^2/2)/w[1:(l-1)])
  p2<- (y*v-v^2/2+y*alt[l]-alt[l]^2/2)/w[l]
  
  out <- p1+p2
  return(out)
}
dR1adap<- function(x,y,w,alt){
  l <- findInterval(min(x,y),alt)
  v <- min(x,alt[l],y)

  ti1 <- alt[2:l]     ### t[i+1]
  ti  <- alt[1:(l-1)] ### t[i]
  
  p1<- sum((-(x-ti1)^2/2+(x-ti)^2/2)/w[1:(l-1)])
  p2<- (x*v-v^2/2+x*alt[l]-alt[l]^2/2)/w[l]
  
  out <- p1+p2
  return(out)
}
ddotR1adap<- function(x,y,w,alt){
  l <- findInterval(min(x,y),alt)
  v <- min(x,alt[l],y)
  
  ti1 <- alt[2:l]     ### t[i+1]
  ti  <- alt[1:(l-1)] ### t[i]
  
  out<- sum((ti1-ti)/w[1:(l-1)])+(v-alt[l])/w[l]

  return(out)
}
ddotdotR1adap<- function(x,y,w,alt){
  l <- findInterval(min(x,y),alt)

  if(x>=y) return(1/w[l])
  else return(0)
}
# dotdotR1adap<- function(x,y,w,alt){
#   if(y>=x) return(y-x)
#   else if(y<x)  return(0)
# }



