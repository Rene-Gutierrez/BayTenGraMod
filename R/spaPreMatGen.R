#' Creates Sparse Precision Matrices
#'
#' This Funtion Creates Sparse Precision Matrices with a shape specified by the
#' user.
#'
#' @param p        A vector specifying the dimension size of each precision
#'   graph, it defaults to a 3 dimesional Tensor with each dimension os size 10.
#' @param type     A charachter strin specifying the matrix structure:
#' \describe{
#'   \item{"R"}{For a Random Matrix, as described in Wang & Li.}
#'   \item{"C"}{For a "Cyclic Graph" Matrix.}
#' }
#' It defaults to a random matrix "R".
#' @param sparsity A vector with the requiered sparsity for the Random Matrix.
#'   If NULL it defaults to a 10\% sparsity level for each precision matrix.
#'
#' @return A list of lists containing the Presion Matrices, the Covariance
#'   Matrices and the Adjacency Matrices.Each list is as follows:
#' \describe{
#' \item{P}{A list of Precision Matrices.}
#' \item{S}{A list of Covariance Matrices.}
#' \item{E}{A list of Adjacency Matrices.}
#' }
#' @author Rene Gutierrez Marquez
#'
#' @export

spaPreMatGen <- function(p        = c(10, 10, 10),
                         type     = "R",
                         sparsity = NULL){
  ### Computes some Global Parameters
  d <- length(p)
  ### Checks defaults
  if(type == "R"){
    if(is.null(sparsity)){
      sparsity <- rep(0.1, d)
    }
    ### Generates the Edge Matrices
    E <- list()
    for(i in 1:d){ # For Every Dimension
      E[[i]] <- edgMat(matSiz = p[i],
                       spa    = sparsity[i])
    }
    ### Sets some parameters to Compute the Random Matrices
    b <- 103  # Degrees of Freedom for Precision Matrices
    ### Generates the Precision Matrices and Covariance Matrices
    P <- list()
    S <- list()
    for(i in 1:d){
      P[[i]] <- Matrix::Matrix(data   = cov2cor(spaPDM(b = b,
                                                       E = E[[i]])),
                               sparse = TRUE)
      S[[i]] <- Matrix::Matrix(data   = solve(P[[i]]),
                               sparse = TRUE)
    }
  }
  ### Returns the Matrices
  return(list(P = P, S = S, E = E))
}

