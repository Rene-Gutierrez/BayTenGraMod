#' Logarithm of J
#'
#' This function computes the function J that appears in the Metropolis Odds as
#' propose in Wang and Li, in Efficient Gaussian graphical Model determination
#' under G-Wishart prior distributions.
#'
#' @param h   Degrees of Freedom for the Wishart Distribution
#' @param B   2 by 2 sub-matrix of the precision matrix.
#' @param a11 Entry 1,1 from the A matrix.
#'
#' @return The Logarithm of J.
#'
#' @author Rene Gutierrez Marquez

###############################################################################
###
### log-constant J
###
### Check documentation for Input and Output Values
### Inputs:
###   h:
###   B:
###   a11:
### Outputs:
###   c: log-constant J
###
###############################################################################

logJ <- function(h, B, a11){
  c <- 1 / 2 * log(2 * pi)
  c <- c - 1 / 2 * log(B[2,2])
  c <- c + (h - 1) / 2 * log(a11)
  c <- c - 1 / 2 * (B[1,1] - B[1,2] / B[2,2] * B[1,2]) * a11
  c <- c - logConsW(h, B[2,2])
  ### Returns the constant
  return(c)
}
