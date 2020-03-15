#' Computes Several performance statistics for the Tensor Gaussian Graphical
#' Model
#'
#' This function computes the following statistics for each model provided.
#' \itemize{
#'   \item Accuracy.
#'   \item True Positive Rate.
#'   \item True Negative Rate.
#' }
#'
#' @param modelList A list of lists. The outer list is indexed by model. The
#'   inner list is indexed by Matrix.
#' @param trueModel A list with the true Matrices
#'
#' @return A list with different performance statistics in comparison with the
#'   true model.
#' \describe{
#'   \item{acc}{Accuracy for each model.}
#'   \item{tpr}{True positive Rate for each model.}
#'   \item{tnr}{True Negative Rate for each model.}
#' }
#'
#' @author Rene Gutierrez Marquez

###############################################################################
###
### Performance Statistics Graph
###
###
###############################################################################

### Number

modComAdj <- function(modelList, trueModel){
  ### Number of Models
  numMod <- length(modelList)
  ### Number of Matrices
  numMat <- length(modelList[[1]])
  ### Dimensions of Matrices
  p <- numeric(numMat)
  for(i in 1:numMat){
    p[i] <- ncol(modelList[[1]][[i]])
  }
  ### Average Accuracy for each Matrix
  acc <- numeric(numMat)
  tpr <- numeric(numMat)
  tnr <- numeric(numMat)
  for(i in 1:numMod){
    accMod <- 0
    truPos <- 0
    conPos <- 0
    truNeg <- 0
    conNeg <- 0
    for(j in 1:numMat){
      accMod <- accMod + sum(modelList[[i]][[j]] == trueModel[[j]])
      truPos <- truPos +
        sum((modelList[[i]][[j]] == trueModel[[j]])[trueModel[[j]] == TRUE])
      conPos <- conPos + sum(trueModel[[j]] == TRUE)
      truNeg <- truNeg +
        sum((modelList[[i]][[j]] == trueModel[[j]])[trueModel[[j]] == FALSE])
      conNeg <- conNeg + sum(trueModel[[j]] == FALSE)
    }
    acc[i] <- (accMod - sum(p)) / 2 /(sum(p * (p - 1) / 2))
    tpr[i] <- (truPos - p[j]) / (conPos - p[j])
    tnr[i] <- (truNeg) / (conNeg)
  }
  ### Return
  return(list(acc = acc, tpr = tpr, tnr = tnr))
}
