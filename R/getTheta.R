#' @title The estimated parameter \eqn{\theta}
#'
#' @description to get the estimated parameter \eqn{\hat{\theta}}.
#' \deqn{f(x; \Phi) = \sum_{j = 1}^g \pi_j f_j(x; \theta_j)}
#'
#' @keywords theta
#' @param datall the data frame
#' @param pa the parameters \eqn{\lambda} and \eqn{\gamma}, which were found by \code{loocvVSP()} function.
#' @param matrix_list the initialized matrices found by \code{iniVSPmatrices()}.
#' @references Z. Cao, D. Bryant, C. Fox, T. Molten and M. Parry. "V-Spline: an Adaptive Smoothing Spline
#' for Trajectory Reconstruction" arXiv preprint arXiv:1803.07184 (2018).
#' @export
getTheta <- function(datall,pa,matrix_list){

  lamb  <- exp(pa[1])
  gam   <- exp(pa[2])

  rowlen<- nrow(datall)
  sub_t <- datall$t # /datall$t[rowlen]

  sub_y <- datall$y
  sub_v <- datall$v

  B     <- matrix_list$B
  C     <- matrix_list$C
  Omega <- matrix_list$Omega

  R     <- solve(t(B)%*%B+gam*t(C)%*%C+Omega*lamb*rowlen)
  theta <- R%*%(t(B)%*%sub_y+gam*t(C)%*%sub_v)

  return(theta)
}
