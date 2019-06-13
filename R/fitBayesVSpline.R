#'
#' @title Bayes estimate of V-Spline model
#' @description
#' @export
fitBayesVSpline <- function(dat,W=NULL,U=NULL,pa=NULL,xout=NULL){
  if(is.null(xout)) x_star <- dat$t
  else x_star <- xout
 
  X<- dat$t
  Y<- dat$y
  V<- dat$v
  GPmean_x  <- numeric(length(x_star))
  GPmean_v  <- numeric(length(x_star))
  GPsig_x   <- numeric(length(x_star))
  fitmu     <- numeric(length(x_star))
  fitvar    <- numeric(length(x_star))
  
  if(is.null(W)) W <- diag(length(X))
  if(is.null(U)) U <- diag(length(X))
  
  if(is.null(pa)){
    pa      <- c(0,0)
    paraGPR <- optim(pa=pa,fn=remlScoreVSpline,X=X,Y=Y,V=V,W=W,U=U)
    pp      <- exp(paraGPR$par)
    coff    <- cofBayesVSpline(X,Y,V,W,U,lamb=pp[1],gam=pp[2])
    for(l in 1:length(x_star)){
      gen        <- BayesVSpline(X,Y,V,coff,x_star[l])
      GPmean_x[l]<- gen$mu
      GPsig_x[l] <- gen$sig/pp[1]/length(X)
      GPmean_v[l]<- gen$mv
    }
    return(list(t=x_star,y=GPmean_x,sig=GPsig_x,v=GPmean_v,par=pp,cv=paraGPR$value))
  }else{
    pp      <- pa
    coff    <- cofBayesVSpline(X,Y,V,W,U,lamb=pp[1],gam=pp[2])
    for(l in 1:length(x_star)){
      gen        <- BayesVSpline(X,Y,V,coff,x_star[l])
      GPmean_x[l]<- gen$mu
      GPsig_x[l] <- gen$sig/pp[1]/length(X)
      GPmean_v[l]<- gen$mv
    }
    return(list(t=x_star,y=GPmean_x,sig=GPsig_x,v=GPmean_v))
  }
 
}



