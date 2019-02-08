#' @title V-Spline basis functions
#'
#' @description V-Spline basis functions allow to you build a V-Spline.
#' The basis functions are denoted by \eqn{N_1(t), N_2(t), \ldots, N_{2n}(t)}.
#'
#' @keywords V-Spline basis functions
#'
#' @param k to determine the location of the estimate.
#' @param t Time sequence.
#' @param est The target estimate.
#'
#' @references Z. Cao, D. Bryant, C. Fox, T. Molten and M. Parry. "V-Spline: an Adaptive Smoothing Spline
#' for Trajectory Reconstruction" arXiv preprint arXiv:1803.07184 (2018).
#' @export
N1= function(k,t,est){

  a=t[1]
  b=t[2]

  if(est>=t[1] & est<t[2])
    fx=2*(est-a)^3/(b-a)^3-3*(est-a)^2/(b-a)^2+1
  else
    fx=0

  return(f=fx)
}
#' @rdname N1
N2= function(k,t,est){

  a=t[1]
  b=t[2]
  if(est>=t[1] & est<t[2])
    fx=(est-a)^3/(b-a)^2-2*(est-a)^2/(b-a)+(est-a)
  else
    fx=0

  return(f=fx)
}
#' @rdname N1
N2k1= function(k,t,est){

  a=t[k]
  b=t[k+1]
  c=t[k+2]

  if(est==t[k+1]) fx=1
  else{

    if(est>=t[k] & est<t[k+1])
      fx=-2*(est-a)^3/(b-a)^3+3*(est-a)^2/(b-a)^2
    else if(est>=t[k+1] & est<t[k+2])
      fx=2*(est-b)^3/(c-b)^3-3*(est-b)^2/(c-b)^2+1
    else
      fx=0
  }

  return(f=fx)
}
#' @rdname N1
N2k2= function(k,t,est){

  a=t[k]
  b=t[k+1]
  c=t[k+2]

  if(est>=t[k] & est<t[k+1])
    fx=(est-a)^3/(b-a)^2-(est-a)^2/(b-a)
  else if(est>=t[k+1] & est<t[k+2])
    fx=(est-b)^3/(c-b)^2-2*(est-b)^2/(c-b)+(est-b)
  else
    fx=0

  return(f=fx)
}
#' @rdname N1
N2n1= function(k,t,est){

  n=length(t)

  a=t[n-1]
  b=t[n]

  if(est>=a & est<=b)
    fx=-2*(est-a)^3/(b-a)^3+3*(est-a)^2/(b-a)^2
  else
    fx=0

  return(f=fx)
}
#' @rdname N1
N2n= function(k,t,est){

  n=length(t)

  a=t[n-1]
  b=t[n]

  if(est>=a & est<=b)
    fx=(est-a)^3/(b-a)^2-(est-a)^2/(b-a)
  else
    fx=0

  return(f=fx)
}
#' @rdname N1
dN1= function(k,t,est){

  a=t[1]
  b=t[2]

  if(est>=t[1] & est<t[2])
    fx=6*(est-a)^2/(b-a)^3-6*(est-a)/(b-a)^2
  else
    fx=0

  return(f=fx)
}
#' @rdname N1
dN2= function(k,t,est){

  a=t[1]
  b=t[2]
  if(est>=t[1] & est<t[2])
    fx=3*(est-a)^2/(b-a)^2-4*(est-a)/(b-a)+1
  else
    fx=0

  return(f=fx)
}
#' @rdname N1
dN2k1= function(k,t,est){

  a=t[k]
  b=t[k+1]
  c=t[k+2]

  if(est==t[k+1]) fx=0
  else{

    if(est>=t[k] & est<t[k+1])
      fx=-6*(est-a)^2/(b-a)^3+6*(est-a)/(b-a)^2
    else if(est>=t[k+1] & est<t[k+2])
      fx=6*(est-b)^2/(c-b)^3-6*(est-b)/(c-b)^2
    else
      fx=0
  }

  return(f=fx)
}
#' @rdname N1
dN2k2= function(k,t,est){

  a=t[k]
  b=t[k+1]
  c=t[k+2]

  if(est>=t[k] & est<t[k+1])
    fx=3*(est-a)^2/(b-a)^2-2*(est-a)/(b-a)
  else if(est>=t[k+1] & est<t[k+2])
    fx=3*(est-b)^2/(c-b)^2-4*(est-b)/(c-b)+1
  else
    fx=0

  return(f=fx)
}
#' @rdname N1
dN2n1= function(k,t,est){

  n=length(t)

  a=t[n-1]
  b=t[n]

  if(est>=a & est<=b)
    fx=-6*(est-a)^2/(b-a)^3+6*(est-a)/(b-a)^2
  else
    fx=0

  return(f=fx)
}
#' @rdname N1
dN2n= function(k,t,est){

  n=length(t)

  a=t[n-1]
  b=t[n]

  if(est>=a & est<=b)
    fx=3*(est-a)^2/(b-a)^2-2*(est-a)/(b-a)
  else
    fx=0

  return(f=fx)
}
#' @rdname N1
d2N1= function(k,t,est){

  a=t[1]
  b=t[2]

  if(est>=t[1] & est<t[2])
    fx=12*(est-a)/(b-a)^3-12/(b-a)^2
  else
    fx=0

  return(f=fx)
}
#' @rdname N1
d2N2= function(k,t,est){

  a=t[1]
  b=t[2]
  if(est>=t[1] & est<t[2])
    fx=6*(est-a)/(b-a)^2-4/(b-a)
  else
    fx=0

  return(f=fx)
}
#' @rdname N1
d2N2k1= function(k,t,est){

  a=t[k]
  b=t[k+1]
  c=t[k+2]

  if(est==t[k+1]) fx=0
  else{

    if(est>=t[k] & est<t[k+1])
      fx=-12*(est-a)/(b-a)^3+6/(b-a)^2
    else if(est>=t[k+1] & est<t[k+2])
      fx=12*(est-b)/(c-b)^3-6/(c-b)^2
    else
      fx=0
  }

  return(f=fx)
}
#' @rdname N1
d2N2k2= function(k,t,est){

  a=t[k]
  b=t[k+1]
  c=t[k+2]

  if(est>=t[k] & est<t[k+1])
    fx=6*(est-a)/(b-a)^2-2/(b-a)
  else if(est>=t[k+1] & est<t[k+2])
    fx=6*(est-b)/(c-b)^2-4/(c-b)
  else
    fx=0

  return(f=fx)
}
#' @rdname N1
d2N2n1= function(k,t,est){

  n=length(t)

  a=t[n-1]
  b=t[n]

  if(est>=a & est<=b)
    fx=-12*(est-a)/(b-a)^3+6/(b-a)^2
  else
    fx=0

  return(f=fx)
}
#' @rdname N1
d2N2n= function(k,t,est){

  n=length(t)

  a=t[n-1]
  b=t[n]

  if(est>=a & est<=b)
    fx=6*(est-a)/(b-a)^2-2/(b-a)
  else
    fx=0

  return(f=fx)
}
