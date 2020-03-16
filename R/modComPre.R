#' Computes Several performance statistics for the Tensor Gaussian Graphical
#' Model
#'
#' This function computes the following statistics for each model provided.
#' \itemize{
#'   \item Frobenious Norm of the Kronecker product.
#'   \item Infinity Norm of the Kronecker product.
#'   \item Average Frobenious Norm of each matrix.
#'   \item Average Infinity Norm of each matrix.
#' }
#'
#' @param modelList A list of lists. The outer list is indexed by model. The
#'   inner list is indexed by Matrix.
#' @param trueModel A list with the true Matrices
#'
#' @return A list with different performance statistics in comparison with the
#'   true model.
#' \describe{
#'   \item{kroFro}{Frobenious Norm of the Kronecker product for each model.}
#'   \item{kroInf}{Infinity Norm of the Kronecker product for each model.}
#'   \item{sinFro}{Average Frobenious Norm of each matrix for each model.}
#'   \item{sinInf}{Average Infinity Norm of each matrix for each model.}
#' }
#'
#' @author Rene Gutierrez Marquez
#'
#' @export

###############################################################################
###
### Performance Statistics Precison Matrix
###
###
###############################################################################

modComPre <- function(modelList, trueModel){
  ### Number of Models
  numMod <- length(modelList)
  ### Number of Matrices
  numMat <- length(modelList[[1]])
  ### Frobenius and Infinity Norm of the Kronecker Product
  kroFro    <- numeric(length = numMod)
  kroInf    <- numeric(length = numMod)
  kroProTru <- mkronecker(trueModel)
  for(i in 1:numMod){
    kroProMod <- mkronecker(modelList[[i]])
    froNor    <- Matrix::norm(kroProMod - kroProTru, "F")
    infNor    <- Matrix::norm(kroProMod - kroProTru, "I")
    kroFro[i] <- froNor
    kroInf[i] <- infNor
  }
  ### Average Frobenius and Infinity Norm of each Matrix
  sinFro    <- numeric(length = numMod)
  sinInf    <- numeric(length = numMod)
  for(i in 1:numMod){
    sinFro[i] <- 0
    sinInf[i] <- 0
    for(j in 1:numMat){
      froNor    <- Matrix::norm(modelList[[i]][[j]] - trueModel[[j]], "F")
      infNor    <- Matrix::norm(modelList[[i]][[j]] - trueModel[[j]], "I")
      sinFro[i] <- sinFro[i] + froNor
      sinInf[i] <- sinInf[i] + infNor
    }
  }
  sinFro <- sinFro / numMat
  sinInf <- sinInf / numMat

  ### Return
  return(list(kroFro = kroFro, kroInf = kroInf,
              sinFro = sinFro, sinInf = sinInf))
}
