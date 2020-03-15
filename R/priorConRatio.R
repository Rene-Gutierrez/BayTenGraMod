#' Logarithm of the Prior Constant Ratio Approximation
#'
#' This function computes the Prior Constant Ratio Approximation that appears in
#' the Metropolis Odds as propose in Wang and Li, in Efficient Gaussian
#' graphical Model determination under G-Wishart prior distributions.
#'
#' @param C      Initial Precision Matrix for the MCMC.
#' @param bPrior Prior dedrees of freedom. Defaults to 3.
#' @param DPrior Prior location SPD Matrix. Defaults to the identity matrix.
#' @param i      Row index.
#' @param j      Column index.
#'
#' @return Logarithm of the Prior Constant Ratio Approximation.
#'
#' @author Rene Gutierrez Marquez

#' @export

###############################################################################
###
### Log Prior Constant Ratio Approximation
###
### Inputs:
###   C:      Covariance Matrix
###   bPrior: Prior dedrees of freedom
###   DPrior: Prior location SPD Matrix
###   i:      Row Index
###   j:      Column Index
###
### Outputs:
###   h:      Log Prior Constant Ratio Approximation
###
###############################################################################

priorConRatio <- function(C,
                       bPrior = 100,
                       DPrior = diag(ncol(C)),
                       i      = 1,
                       j      = 1){
  ### Case Entry (i,j) = 0
  e        <- j # Only one entry needs to be evaluated
  ### We have to compute a similar precision matrix to C but:
  ### entries zero at (i,j) and (j,i) and c at (j,j)
  C0       <- C
  C0[i, e] <- 0
  C0[e, i] <- 0
  ### Auxilary Variables to compute c
  C12 <- C0[-e, e] # Column of C0 without e
  C22 <- C0[-e,-e] # Matrix with row e and column e
  ### Computes c
  c        <- t(C12) %*% solve(C22, C12)
  ### Creates C0
  C0[e, e] <- c

  ### Entry Case (i,j) = 1
  e <- c(i,j)
  ### We have to compute a similar precision matrix to C but:
  ### entries in (i, i), (i, j), (j, i) and (i, i)
  ### Auxilary variables to compute all the new entries toghether
  C1  <- C
  C12 <- C[-e, e]
  C22 <- C[-e,-e]
  ### Computes an auxilary Variable
  if(length(C22) == 1){
    C12 <- matrix(C12, c(1,2))
    k <- t(C12) %*% C12 / C22
  } else {
    k <- t(C12) %*% solve(C22, C12)
  }
  ### Creates C1
  C1[e,e] <- k
  ### Creates A
  A <- C[e,e] - k

  ### Computes log H
  h <- -logConsW(v = bPrior,
                 S = DPrior[j, j])
  h <- h - logJ(h   = bPrior,
                B   = DPrior[e,e],
                a11 = A[1,1])
  h <- h + (bPrior -2) / 2 * log(A[1,1])
  h <- h - sum(diag(DPrior[e, e] %*% (C0[e, e] - C1[e, e]))) / 2
  ### Returns the value of the constant
  return(-h)
}
