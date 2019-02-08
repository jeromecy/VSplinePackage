#' @title Initialise V-Spline matrices \eqn{B}, \eqn{C} and \eqn{\Omega}
#'
#' @description This function is utilised for initializing the V-Spline matreces \eqn{B},
#' \eqn{C} and \eqn{\Omega}, with whitch it potentially saves computating time.
#'
#' The default values of \eqn{a} and \eqn{b} are 3 and 2 respectively.
#'
#' @keywords V-Spline matrices
#'
#' @export
iniVSPmat <- function(dat,a=3,b=2){
  rowlen   <- nrow(dat)
  loo_t    <- dat$t # /dat$t[rowlen]
  loo_x    <- dat$x

  B      <- matrix(0,nrow=rowlen,ncol=2*rowlen)
  for(k in 1:(rowlen-2)) B[k+1,2*k+1] <- 1
  B[1,1] <- 1
  B[rowlen,2*rowlen-1] <- 1

  C      <- matrix(0,nrow=rowlen,ncol=2*rowlen)
  for(k in 1:(rowlen-2)) C[k+1,(2*k+2)] <- 1
  C[1,2] <- 1
  C[rowlen,2*rowlen] <- 1

  d    <- loo_t[2:rowlen]-loo_t[1:(rowlen-1)]
  dist <- abs(diff(loo_x))
  dist[which(dist==0)] <- 1e-4

  Omega <- matrix(0,nrow=2*rowlen,ncol=2*rowlen)

  for(i in 1:(rowlen-1)){
    OmegaTemp <- matrix(0,nrow=2*rowlen,ncol=2*rowlen)

    OmegaTemp[2*i-1,2*i-1] <- 12/d[i]^3
    OmegaTemp[2*i-1,2*i]   <- OmegaTemp[2*i,2*i-1]   <- 6/d[i]^2
    OmegaTemp[2*i-1,2*i+1] <- OmegaTemp[2*i+1,2*i-1] <- -12/d[i]^3
    OmegaTemp[2*i-1,2*i+2] <- OmegaTemp[2*i+2,2*i-1] <- 6/d[i]^2

    OmegaTemp[2*i,2*i]     <- 4/d[i]
    OmegaTemp[2*i,2*i+1]   <- OmegaTemp[2*i+1,2*i]   <- -6/d[i]^2
    OmegaTemp[2*i,2*i+2]   <- OmegaTemp[2*i+2,2*i]   <- 2/d[i]

    OmegaTemp[2*i+1,2*i+1] <- 12/d[i]^3
    OmegaTemp[2*i+1,2*i+2] <- OmegaTemp[2*i+2,2*i+1] <- -6/d[i]^2

    OmegaTemp[2*i+2,2*i+2] <- 4/d[i]

    Omega <- Omega+OmegaTemp*(d[i]^a/dist[i]^b)
  }
  return(list(B=B,C=C,Omega=Omega))
}
