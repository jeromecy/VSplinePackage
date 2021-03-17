


# Overview

The **VSplinePackage** provides an interface to fit V-Spline and adaptive V-Spline models in **R**.

Two different methods **fitVSP** (fit V-Spline) and **fitBayesVSpline** (fit V-Spline in Bayesian way) are identical with the same parameters and can address the same problem.

This is an open source for **V-Spline** **R** package. Everyone is welcome to contribute on the improvement of this package, comments on coding and fixing bugs.

Thank you.

# How to install

```r
library(devtools)
install_github("jeromecy/VSplinePackage")
```
# Simulated data
```r
library(waveband)
library(VSPline)
set.seed(2016)
N = 100
velocity    <- test.data(type = "blocks", n = N, signal = 1, rsnr = 7, plotfn = TRUE)
set.seed(2016)
position    <- SimuData(velocity$y,7)
simulatedata<- data.frame(t=velocity$x,y=position$xnoise,v=velocity$ynoise,boom=0)
```
Available types are "blocks", "bumps", "heavi" (heavisine), and "doppler".

```r
fitted <- fitVSP(simulatedata)
plot(simulatedata$t,simulatedata$y)
points(fitted$vsp$t,fitted$vsp$y,col="red",type="l")
```
# A simple example
```r
set.seed(1234)
n <- 100
s <- seq(1,2*pi,length=n)
y <- sin(s) + rnorm(n)
v <- cos(s) + rnorm(n)
simuData <- data.frame(t=s,y=y,v=v)
fitted   <- fitVSP(simuData)
plot(s,y)
points(fitted$vsp$t,fitted$vsp$y,col="red",type="l")
```
--------
