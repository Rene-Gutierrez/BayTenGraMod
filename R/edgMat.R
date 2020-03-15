#' Generates a Random Adjacency or Edge Matrix
#'
#' This function generates Random Adjacency or Edge Matrix of a specific size
#' and degree of sparsity.
#'
#' @param matSiz A scalar specifying the size of the matrix. As the result is
#'               a square matrix only a scalar is requiered.
#' @param spa    The degree of sparsity specified with a number from 0 to 1.
#' @param upp    An indicator indicating if the result should be an Upper
#'               Triangular Matrix. The default is TRUE.
#'
#' @return A Random Square Matrix encoding the Adjacency of vertices on a
#'         graph. If upp is TRUE (defualt) the result would by an upper
#'         triangular matrix.
#'
#' @author Rene Gutierrez Marquez
#'
#' @export

###############################################################################
###
### Random Edge Matrix Creator
### Inputs:
### matSiz: Matrix Size
### spa:    Sparsity Level
### upp:    Indicator if only the upper triangular part is desired
###
###############################################################################

edgMat <- function(matSiz, spa, upp = TRUE){
  ### Generates a random number of edges based on the sparsity level
  edg <- rbinom(n = matSiz * (matSiz - 1) / 2,
                size = 1,
                prob = spa)
  ### Obtains the Indexes of Upper part of the Marix
  indSeq <- seq(from = 1, to = matSiz, by = 1)
  indSeq <- indSeq + (indSeq - 1) * matSiz
  uppMat <- c()
  for(j in 1:(matSiz - 1)){
    uppMat <- c(uppMat, seq(from = indSeq[j] + 1, to = j * matSiz, by = 1))
  }
  ### Creates a vector to fill the matrix with zeros and ones
  auxVec <- rep(0, matSiz * matSiz)
  auxVec[uppMat[edg == 1]] = 1
  ### Fills the Upper Part of a matrix
  spaMat <- matrix(data  = auxVec,
                   nrow  = matSiz,
                   ncol  = matSiz,
                   byrow = TRUE)
  ### Makes it symmetric and Fills the Diagonal (If necessary)
  if(!upp){
    spaMat <- spaMat + t(spaMat)
  }
  diag(spaMat) <- rep(1, matSiz)
  ### Returns the Matrix
  return(spaMat)
}
