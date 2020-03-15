#' Logarithm of the Normalizing Constant of the Wishart Distribution
#'
#' This function computes the Logarithm of the Normalizing Constant of the
#' Wishart Distribution
#'
#' @param S Scale matrix parameter.
#' @param v Degrees of freedom.
#'
#' @return The Logarithm of the Normalizing Constant of the Wishart
#'   Distribution
#'
#' @author Rene Gutierrez Marquez

###############################################################################
###
### Log Normalizing Constant of the Wishart Distribution
###
### Inputs:
###   S: Scale Matrix
###   v: degrees of Freedom (note v is always)
###
### Outputs:
###   c: log-normalizing Constant
###
###############################################################################

logConsW <- function(v, S){
  ### Checks if S is a Matrix
  if(!is.matrix(S)){
    S <- matrix(S)
  }
  ### Size of the Matrix S
  p <- ncol(S)
  n <- v + p - 1 # Degrees of Freedom in the Wikipedia Notation (Cleaner Computations)
  ### Computes the Value of the normalizing constant
  c <- - n * p * log(2) / 2
  c <- c + n * log(det(S)) / 2
  c <- c - lMVGamma(p, n / 2)
  ### Returns the Normalizing Constant
  return(c)
}
