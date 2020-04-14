#' Performs Tensor Lasso
#'
#' This function performs Tensor Lasso under different specifications.
#'
#' @param t      List of Sample Tensors
#' @param iter   Number of iterations
#' @param lambda A vector to use as the penalties for each precision matrix.
#' @param identi A method to obtain identifiability.
#'   'F' for Frobenius norm normalization.
#'   'C' for correlation normalization.
#'   '1' for dividing by the first entry of the precision matrix.
#' @param norm   A boolena indicating if the norm of the precision matrices
#'   should remain constant through the iterations.
#' @param search A boolean indicating if a penalty search based on likelihood
#'   should be performed.
#'
#' @return List containg the precision matrices.
#'
#' @author Rene Gutierrez Marquez

#' @export

TLasso <- function(t,
                   iter   = 1,
                   lambda = NULL,
                   identi = 'F',
                   norm   = FALSE,
                   search = FALSE){
  ### Computes some Global Parameters
  ### Number of Samples
  numSam <- length(t)
  ### Dimensions of the Tensors
  p      <- dim(t[[1]])
  ### Number of dimensions of the Tensors
  numDim <- length(p)
  ### Tlasso Replication
  Ohat <- list()
  for (i in 1:d) {
    Ohat[[i]] <- diag(p[i])
  }
  dNor <- norm(Ohat[[1]], 'F')
  for (j in 1:iter) {
    for (i in 1:d) {
      ### Computes the Conditional Sufficient Statistic of the current
      ### presicion matrix
      sK <- Sk(tensors    = t,
               precisions = Ohat,
               k          = i)
      sK <- sK / n / prod(p[-i])
      ### Computes the presicion matrix using Huge
      out <- huge::huge(x = sK,
                        lambda = lambda[i],
                        method = "glasso",
                        verbose = FALSE)
      ### Adjusts the presicion matrix for identifiability purposes
      if(identi == "F"){
        Ohat[[i]] <- as.matrix(out$icov[[1]]) / norm(as.matrix(out$icov[[1]]),
                                                     type = "F")
      } else if(identi == "C") {
        Ohat[[i]] <- cov2cor(out$icov[[1]])
      } else if(identi == "1") {
        Ohat[[i]] <- as.matrix(out$icov[[1]])
        Ohat[[i]] <- Ohat[[i]] / Ohat[[i]][1,1]
      }
      if(norm){
        Ohat[[i]] <- Ohat[[i]] / norm(Ohat[[i]], 'F') * dNor
      }
    }
  }
  return(Ohat)
}
