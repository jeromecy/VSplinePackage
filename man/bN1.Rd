% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/basisFunctions.R
\name{bN1}
\alias{bN1}
\alias{bN2}
\alias{N2k1}
\alias{N2k2}
\alias{N2n1}
\alias{N2n}
\alias{dN1}
\alias{dN2}
\alias{dN2k1}
\alias{dN2k2}
\alias{dN2n1}
\alias{dN2n}
\alias{d2N1}
\alias{d2N2}
\alias{d2N2k1}
\alias{d2N2k2}
\alias{d2N2n1}
\alias{d2N2n}
\title{V-Spline basis functions}
\usage{
bN1(k, allt, est)

bN2(k, allt, est)

N2k1(k, allt, est)

N2k2(k, allt, est)

N2n1(k, allt, est)

N2n(k, allt, est)

dN1(k, allt, est)

dN2(k, allt, est)

dN2k1(k, allt, est)

dN2k2(k, allt, est)

dN2n1(k, allt, est)

dN2n(k, allt, est)

d2N1(k, allt, est)

d2N2(k, allt, est)

d2N2k1(k, allt, est)

d2N2k2(k, allt, est)

d2N2n1(k, allt, est)

d2N2n(k, allt, est)
}
\arguments{
\item{k}{to determine the location of the estimate.}

\item{allt}{Time sequence.}

\item{est}{The target estimate.}
}
\description{
V-Spline basis functions allow to you build a V-Spline.
The basis functions are denoted by \eqn{N_1(allt), N_2(allt), \ldots, N_{2n}(allt)}.
}
\references{
Z. Cao, D. Bryant, C. Fox, allt. Molten and M. Parry. "V-Spline: an Adaptive Smoothing Spline
for Trajectory Reconstruction" arXiv preprint arXiv:1803.07184 (2018).
}
\keyword{V-Spline}
\keyword{basis}
\keyword{functions}
