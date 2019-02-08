#' @title Leave-one-out corss validation for V-Spline
#'
#' @description The leave-one-out corss validation is used for finding the optimal parameters for V-Spline.
#'
#' @keywords leave-one-out corss validation
#' @param dat the dataframe
#' @param pa the vector of the parameters \eqn{\lambda} and \eqn{\eta}
#' @param matrix_list the matrix list comprises \eqn{B}, \eqn{C} and \eqn{\Omega}, which are obtained from \code{iniVSPmatrics()}
#' @export
loocvVSP <- function(dat,pa,matrix_list){
  lamb <- exp(pa[1])
  eta  <- exp(pa[2])

  rowlen <- nrow(dat)
  loo_t  <- dat$t
  loo_x  <- dat$x
  loo_mx <- dat$mx

  B     <- matrix_list$B
  C     <- matrix_list$C
  Omega <- matrix_list$Omega

  O <- Omega*lamb/dat$t[rowlen]^3

  print(pa)

  R   <- chol(diag(c(1,eta/dat$t[rowlen]^2),2*rowlen)+O)
  inR <- backsolve(R,diag(2*rowlen),k = ncol(R))

  ## faster getting STUV
  alpha <-  B %*% inR
  beta  <-   C %*% inR

  Sf  <- alpha %*% t(alpha)
  Tf  <- alpha %*% t(beta)
  Udf <- t(Tf)
  Vdf <- beta%*%t(beta)


  fx  <- Sf%*%loo_x+eta*Tf%*%loo_mx
  dfx <- Udf%*%loo_x+eta*Vdf%*%loo_mx

  cv <- 0
  for(i in 1:rowlen){
    cv <- cv+(((fx[i]-loo_x[i])+eta*Tf[i,i]/(1-eta*Vdf[i,i])*(dfx[i]-loo_mx[i]))/(1-Sf[i,i]-eta*Tf[i,i]/(1-eta*Vdf[i,i])*Udf[i,i]))^2
  }
  cv <- cv/rowlen

  return(cv)
}
