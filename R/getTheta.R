#' @title The estimated parameter \eqn{\theta}
#'
#' @description to get the estimated parameter \eqn{\hat{\theta}}.
#' \deqn{f(x; \Phi) = \sum_{j = 1}^g \pi_j f_j(x; \theta_j)}
#'
#' @keywords theta
#' @param data_frame the data frame
#' @param pa the parameters \eqn{\lambda} and \eqn{\gamma}, which were found by \code{loocvVSP()} function.
#' @param matrix_list the initialized matrices found by \code{iniVSPmatrices()}.
#' @references Z. Cao, D. Bryant, C. Fox, T. Molten and M. Parry. "V-Spline: an Adaptive Smoothing Spline
#' for Trajectory Reconstruction" arXiv preprint arXiv:1803.07184 (2018).
#' @export
getTheta <- function(data_frame,pa,matrix_list){

  lamb   <- exp(pa[1])
  eta    <- exp(pa[2])

  n      <- nrow(data_frame)
  sub_t  <- data_frame$t # /data_frame$t[n]

  sub_x  <- data_frame$x
  sub_mx <- data_frame$mx

  B     <- matrix_list$B
  C     <- matrix_list$C
  Omega <- matrix_list$Omega

  # O     <- Omega*lamb #/data_frame$t[n]^3
  # R     <- chol(diag(c(1,eta/data_frame$t[n]^2),2*n)+O)
  # theta <- backsolve(R,forwardsolve(t(R),t(B)%*%sub_x/n+eta/data_frame$t[n]^2*t(C)%*%sub_mx/n))
  R     <- solve(t(B)%*%B+eta*t(C)%*%C+Omega*n)
  theta <- R%*%(t(B)%*%sub_x+eta*t(C)%*%sub_mx)
  return(theta)
}
