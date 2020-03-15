#' Generates Sparse Symmetric Positive Definite Random Matrices.
#'
#' This function generates Sparse Symmetric Positive Definitive Matrices
#' according to the adjacency matrix provided and certain degrees of freedom.
#' The whole procedure follows, ...
#'
#' @param b Degrees of Freedom for Sampling from the Gwishart distribution.
#' @param E The Edge or Adjacency matrix for sampling from the GWishart
#'          distribution.
#'
#' @return A random Sparse Symmetric Positive Definite matrix that has zero
#'         entries whenever the adjacency matrix has a zero entry.
#'
#' @author Rene Gutierrez Marquez
#'
#' @export

###############################################################################
###
### Sparse Positve Definite Matrix Creator according to case 2 of
### Efficient Gaussian graphical model determination under G-Wishart prior
### distributions
###
### Requieres DBGraph library
###
### Inputs:
### b: GWishart Degrees of Freedom
### E: Edge Matrix
###
###############################################################################

spaPDM <- function(b, E){
  ### Matrix Size
  p <- ncol(E)
  ### Checks if there is a off diagonal entry
  if(sum(E) == p){
    D <- diag(p)
  } else {
    ### Auxilary Matrices
    B       <- (E + t(E)) * 0.5
    diag(B) <- rep(0, p)
    J1 <- B + 10 * diag(p)
    ### Biggest Eigenvalue of J1
    eigVal <- eigen(J1)$values

    lmax   <- eigVal[1]
    lmin   <- eigVal[p]
    del    <- (lmax - p * lmin) / (p - 1)
    ### New Auxilary Matrix
    J <- J1 + del * diag(p)
    D <- diag(p) + 100*solve(J)
  }
  ### Generates the Random Matrix
  M <- BDgraph::rgwish(n   = 1,
                       adj = E,
                       b   = b,
                       D   = D)
  M <- round(M, 7)
  ### Returns the Matrix
  return(M)
}
