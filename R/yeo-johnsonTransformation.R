
#' @title Logarithmic transformation function.
#'
#' @description Computes the logarithm of a number.
#'
#' @param y Numerical value of which the logarithm is computed.
#'
#' @return This function returns the logarithm of the provided argument y
#' if it is greater than zero. If y is smaller than zero, it will return 0.
#'
#' @noRd
#'

log_transform = function(y) {
  transform_y = (y>0)*y+(y<=0)*1
  return(log(transform_y))
}

#' @title Power transformation function.
#'
#' @description Computes a given power of a number.
#'
#' @param y The number which one wants to raise to a certain power \code{pw}.
#' @param pw The power to which to raise \code{y}.
#'
#' @return This function returns the result of raising \code{y} to the power
#' \code{pw} when \code{y > 0}. Otherwise, it will return 1.
#'
#' @noRd
#'

power_transform = function(y,pw) {
  transform_y = (y>0)*y+(y<=0)*1
  return(transform_y^pw)
}

#' @title Yeo-Johnson transformation function
#'
#' @description Computes the Yeo-Johnson transformation of the provided argument.
#'
#' @param y The argument to be supplied to the Yeo-Johnson transformation.
#' @param theta The parameter of the Yeo-Johnson transformation. This should be
#' a number in the range [0,2].
#'
#' @return The transformed value of y.
#'
#' @noRd
#'

YJtrans = function(y,theta) {
  sg = y>=0
  if (theta==0) {temp = log_transform(y+1)*sg+(1-sg)*(0.5-0.5*(y-1)^2)}
  if (theta==2) {temp = sg*(-0.5+0.5*(y+1)^2)-log_transform(-y+1)*(1-sg)}
  if ((theta!=0) & (theta!=2)) {temp = sg*(power_transform(y+1,theta)-1)/theta+(1-sg)*(1-power_transform(-y+1,2-theta))/(2-theta)}
  return(temp)
}

#' @title Inverse Yeo-Johnson transformation function
#'
#' @description Computes the inverse Yeo-Johnson transformation of the provided
#' argument.
#'
#' @param y The argument to be supplied to the inverse Yeo-Johnson transformation.
#' @param theta The parameter of the inverted Yeo-Johnson transformation. This
#' should be a number in the range [0,2].
#'
#' @return The transformed value of y.
#'
#' @noRd
#'

IYJtrans = function(y,theta) {
  sg = y>=0
  if (theta==0) {temp =(exp(y)-1)*sg+(1-sg)*(1-power_transform(-2*y+1,0.5))}
  if (theta==2) {temp = sg*(-1+power_transform(2*y+1,0.5))+(1-exp(-y))*(1-sg)}
  if ((theta!=0) & (theta!=2)) {temp = sg*(power_transform(abs(theta)*y+1,1/theta)-1)+(1-sg)*(1-power_transform(1-(2-theta)*y,1/(2-theta)))}
  return(temp)
}

#' @title Derivative of the Yeo-Johnson transformation function
#'
#' @description Evaluates the derivative of the Yeo-Johnson transformation at
#' the provided argument.
#'
#' @param y The argument to be supplied to the derivative of the Yeo-Johnson
#' transformation.
#' @param theta The parameter of the Yeo-Johnson transformation. This should be
#' a number in the range [0,2].
#'
#' @return The transformed value of y.
#'
#' @noRd
#'

DYJtrans = function(y,theta) {
  sg = y>=0
  temp = power_transform(y+1,theta-1)*sg+power_transform(-y+1,1-theta)*(1-sg)
  return(temp)
}

