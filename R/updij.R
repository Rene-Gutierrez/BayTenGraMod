#' Updates the Entries of the Precision Matrix
#'
#' This function updates the Entries (j,j), (i,j) and (j,i) of the Precision Matrix according to the proposed Edge
#'
#' @param C      Initial Precision Matrix for the MCMC.
#' @param bPrior Prior dedrees of freedom. Defaults to 3.
#' @param DPrior Prior location SPD Matrix. Defaults to the identity matrix.
#' @param i      Row index.
#' @param j      Column index.
#' @param proij  Boolean indicatting an edgeif TRUE and no edge if FALSE.
#'
#' @return An updated Precision matrix.
#'
#' @author Rene Gutierrez Marquez

#' @export

###############################################################################
###
### Updates the Entries of the Covariance Matrix according to the proposed edge
###
### Inputs:
###   C:     Covariance Matrix
###   bPost: Posterior Degrees of Freedom
###   DPost: Posterior Scale Matrix
###   i:     Row
###   j:     Column
###   proij: Indicator for the edge between i and j
###
### Outputs:
###   C: Updated Covariance Matrix
###
###############################################################################

updij <- function(C, bPost, DPost, i, j, proij){
  ### Checks if the proposed edge is there or not
  if(proij){ # There is an edge
    if(length(C[-c(i,j),-c(i,j)]) == 1){
      Cij <- matrix(data = C[c(i,j),-c(i,j)], c(1,2))
      FF <- t(Cij) %*% Cij / C[-c(i,j),-c(i,j)]
    } else {
      FF <- C[c(i,j),-c(i,j)] %*% solve(C[-c(i,j),-c(i,j)], C[-c(i,j),c(i,j)])
    }
    A     <- C[c(i,j), c(i,j)] - FF
    u     <- rnorm(n    = 1,
                   mean = -DPost[i,j] * A[1,1] / DPost[j,j],
                   sd   = sqrt(A[1,1] / DPost[j,j]))
    v     <- BDgraph::rwish(n = 1,
                      p = 1,
                      b = bPost,
                      D = as.matrix(DPost[j,j]))
    C[i,j] <- u + FF[1,2]
    C[j,i] <- C[i,j]
    C[j,j] <- v + u^2 / A[1,1] + FF[2,2]
  } else { # No edge
    C[i,j] <- 0
    C[j,i] <- 0
    c <- C[j,-j] %*% solve(C[-j,-j], C[-j,j])
    u <- BDgraph::rwish(n = 1,
                  p = 1,
                  b = bPost,
                  D = as.matrix(DPost[j,j]))
    C[j,j] <- u + c
  }
  ### Return the whole Covariance Matrix
  return(C)
}
