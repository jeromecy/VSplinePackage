#' @title Forward step for generating simulation data
#' @description Forward step.
#' forward1()
forward1 = function(t1,t2,v1,v2,x1){
  return(x1+(v1+v2)*(t2-t1)/2)
}
