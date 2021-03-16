


# Overview

The **VSplinePackage** provides an interface to fit V-Spline and adaptive V-Spline models in **R**.

Two different methods **fitVSP** (fit V-Spline) and **fitBayesVSpline** (fit V-Spline in a Bayesian way) are identical in some way and could address the same problem. 

This is an open source for V-Spline R package. Everyone is welcome to contribute on the improvement of this package, comments on coding and fixing bugs.

Thank you.

# How to install

```r
library(devtools)
install_github("jeromecy/VSplinePackage")
```
# Simulated data
```r
library(waveband)
set.seed(2016)
velocity17<- test.data(type = "blocks", n = 1024, signal = 1, rsnr = 7, plotfn = TRUE)
set.seed(2016)
position17=SimuData(velocity17$y,7)
simulate17data<- data.frame(t=velocity17$x,y=position17$xnoise,v=velocity17$ynoise,boom=0)
```
Available types are "blocks", "bumps", "heavi" (heavisine), and "doppler".

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
