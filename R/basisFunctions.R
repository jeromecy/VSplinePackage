#' @title V-Spline basis functions
#'
#' @description V-Spline basis functions allow to you build a V-Spline.
#' The basis functions are denoted by \eqn{N_1(allt), N_2(allt), \ldots, N_{2n}(allt)}.
#'
#' @keywords V-Spline basis functions
#'
#' @param k to determine the location of the estimate.
#' @param allt Time sequence.
#' @param est The target estimate.
#'
#' @references Z. Cao, D. Bryant, C. Fox, allt. Molten and M. Parry. "V-Spline: an Adaptive Smoothing Spline
#' for Trajectory Reconstruction" arXiv preprint arXiv:1803.07184 (2018).
#' @export
bN1= function(k,allt,est){

  a=allt[1]
  b=allt[2]

  if(est>=allt[1] & est<allt[2])
    fx=2*(est-a)^3/(b-a)^3-3*(est-a)^2/(b-a)^2+1
  else
    fx=0

  return(f=fx)
}
#' @rdname bN1
bN2= function(k,allt,est){

  a=allt[1]
  b=allt[2]
  if(est>=allt[1] & est<allt[2])
    fx=(est-a)^3/(b-a)^2-2*(est-a)^2/(b-a)+(est-a)
  else
    fx=0

  return(f=fx)
}
#' @rdname bN1
N2k1= function(k,allt,est){

  a=allt[k]
  b=allt[k+1]
  c=allt[k+2]

  if(est==allt[k+1]) fx=1
  else{

    if(est>=allt[k] & est<allt[k+1])
      fx=-2*(est-a)^3/(b-a)^3+3*(est-a)^2/(b-a)^2
    else if(est>=allt[k+1] & est<allt[k+2])
      fx=2*(est-b)^3/(c-b)^3-3*(est-b)^2/(c-b)^2+1
    else
      fx=0
  }

  return(f=fx)
}
#' @rdname bN1
N2k2= function(k,allt,est){

  a=allt[k]
  b=allt[k+1]
  c=allt[k+2]

  if(est>=allt[k] & est<allt[k+1])
    fx=(est-a)^3/(b-a)^2-(est-a)^2/(b-a)
  else if(est>=allt[k+1] & est<allt[k+2])
    fx=(est-b)^3/(c-b)^2-2*(est-b)^2/(c-b)+(est-b)
  else
    fx=0

  return(f=fx)
}
#' @rdname bN1
N2n1= function(k,allt,est){

  n=length(allt)

  a=allt[n-1]
  b=allt[n]

  if(est>=a & est<=b)
    fx=-2*(est-a)^3/(b-a)^3+3*(est-a)^2/(b-a)^2
  else
    fx=0

  return(f=fx)
}
#' @rdname bN1
N2n= function(k,allt,est){

  n=length(allt)

  a=allt[n-1]
  b=allt[n]

  if(est>=a & est<=b)
    fx=(est-a)^3/(b-a)^2-(est-a)^2/(b-a)
  else
    fx=0

  return(f=fx)
}
#' @rdname bN1
dN1= function(k,allt,est){

  a=allt[1]
  b=allt[2]

  if(est>=allt[1] & est<allt[2])
    fx=6*(est-a)^2/(b-a)^3-6*(est-a)/(b-a)^2
  else
    fx=0

  return(f=fx)
}
#' @rdname bN1
dN2= function(k,allt,est){

  a=allt[1]
  b=allt[2]
  if(est>=allt[1] & est<allt[2])
    fx=3*(est-a)^2/(b-a)^2-4*(est-a)/(b-a)+1
  else
    fx=0

  return(f=fx)
}
#' @rdname bN1
dN2k1= function(k,allt,est){

  a=allt[k]
  b=allt[k+1]
  c=allt[k+2]

  if(est==allt[k+1]) fx=0
  else{

    if(est>=allt[k] & est<allt[k+1])
      fx=-6*(est-a)^2/(b-a)^3+6*(est-a)/(b-a)^2
    else if(est>=allt[k+1] & est<allt[k+2])
      fx=6*(est-b)^2/(c-b)^3-6*(est-b)/(c-b)^2
    else
      fx=0
  }

  return(f=fx)
}
#' @rdname bN1
dN2k2= function(k,allt,est){

  a=allt[k]
  b=allt[k+1]
  c=allt[k+2]

  if(est>=allt[k] & est<allt[k+1])
    fx=3*(est-a)^2/(b-a)^2-2*(est-a)/(b-a)
  else if(est>=allt[k+1] & est<allt[k+2])
    fx=3*(est-b)^2/(c-b)^2-4*(est-b)/(c-b)+1
  else
    fx=0

  return(f=fx)
}
#' @rdname bN1
dN2n1= function(k,allt,est){

  n=length(allt)

  a=allt[n-1]
  b=allt[n]

  if(est>=a & est<=b)
    fx=-6*(est-a)^2/(b-a)^3+6*(est-a)/(b-a)^2
  else
    fx=0

  return(f=fx)
}
#' @rdname bN1
dN2n= function(k,allt,est){

  n=length(allt)

  a=allt[n-1]
  b=allt[n]

  if(est>=a & est<=b)
    fx=3*(est-a)^2/(b-a)^2-2*(est-a)/(b-a)
  else
    fx=0

  return(f=fx)
}
#' @rdname bN1
d2N1= function(k,allt,est){

  a=allt[1]
  b=allt[2]

  if(est>=allt[1] & est<allt[2])
    fx=12*(est-a)/(b-a)^3-12/(b-a)^2
  else
    fx=0

  return(f=fx)
}
#' @rdname bN1
d2N2= function(k,allt,est){

  a=allt[1]
  b=allt[2]
  if(est>=allt[1] & est<allt[2])
    fx=6*(est-a)/(b-a)^2-4/(b-a)
  else
    fx=0

  return(f=fx)
}
#' @rdname bN1
d2N2k1= function(k,allt,est){

  a=allt[k]
  b=allt[k+1]
  c=allt[k+2]

  if(est==allt[k+1]) fx=0
  else{

    if(est>=allt[k] & est<allt[k+1])
      fx=-12*(est-a)/(b-a)^3+6/(b-a)^2
    else if(est>=allt[k+1] & est<allt[k+2])
      fx=12*(est-b)/(c-b)^3-6/(c-b)^2
    else
      fx=0
  }

  return(f=fx)
}
#' @rdname bN1
d2N2k2= function(k,allt,est){

  a=allt[k]
  b=allt[k+1]
  c=allt[k+2]

  if(est>=allt[k] & est<allt[k+1])
    fx=6*(est-a)/(b-a)^2-2/(b-a)
  else if(est>=allt[k+1] & est<allt[k+2])
    fx=6*(est-b)/(c-b)^2-4/(c-b)
  else
    fx=0

  return(f=fx)
}
#' @rdname bN1
d2N2n1= function(k,allt,est){

  n=length(allt)

  a=allt[n-1]
  b=allt[n]

  if(est>=a & est<=b)
    fx=-12*(est-a)/(b-a)^3+6/(b-a)^2
  else
    fx=0

  return(f=fx)
}
#' @rdname bN1
d2N2n= function(k,allt,est){

  n=length(allt)

  a=allt[n-1]
  b=allt[n]

  if(est>=a & est<=b)
    fx=6*(est-a)/(b-a)^2-2/(b-a)
  else
    fx=0

  return(f=fx)
}
