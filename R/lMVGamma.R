#' Logarithm of the Multivariate Gamma Function
#'
#' This function computes the Multivariate Gamma Function.
#'
#' @param p Dimension parameter.
#' @param a Real number to evaluate the function.
#'
#' @return The Logarithm of the Multivariate Gamma Function.
#'
#' @author Rene Gutierrez Marquez

###############################################################################
###
### Log Multivariate Gamma Function
###
### Inputs:
### p: Dimension
### a: Value to be evaluated
###
### Outputs:
### l: log multivariate gamma function
###
###############################################################################

lMVGamma <- function(p = 1, a = 1){
  ### Computes the Value according to Wikipedia
  l <- (p * (p -1) / 4) * log(pi)
  for(i in 1:p){
    l <- l + lgamma(a + (1 - i) / 2)
  }
  ### Returns the log multivariate gamma function
  return(l)
}
