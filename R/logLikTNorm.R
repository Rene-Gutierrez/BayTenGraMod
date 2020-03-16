#' Computes the log-likelihood of the Tesor Normal Distribution
#'
#' Computes the log-likelihood of the Tesor Normal Distribution given a sample
#' of tensors.
#'
#' @param tensors    A list with the sample of tensors
#' @param precisions A list of precision matrices.
#'
#' @return The log-likelihood for the Tensor Normal distribution precisions and
#'   sample specified.
#'
#' @author Rene Gutierrez Marquez
#'
#' @export

###############################################################################
###
### Log-Likelihood for the Tensor Normal Distribution of 1 Tensor
###
### Inputs:
### tensor:     Tensor
### precisions: List of Precision Matrices
###
### Outputs:
### logL: Log Likelihood of the
###
###############################################################################

logLikTNorm <- function(tensors, precisions){
  ### Gets the Number of Dimensions of the Tensor
  d <- length(precisions)
  ### Gets the Size of each Matrix
  p <- numeric()
  for(i in 1:d){
    p <- c(p, ncol(precisions[[i]]))
  }
  ### Computes the conditional Suffiecent Statistic Sk, for Precision Matrix k,
  ### where k is the Precision matrix of highest dimension
  ### ### Finds the precision matrix of highest dimension
  k <- which.max(p)
  ### ### Obtains Sk
  sS <- Sk(tensors    = tensors,
           precisions = precisions,
           k          = k)
  ### Computes the Likelihood
  logL <- - 1 / 2 * sum((sS %*% precisions[[k]])[diag(p[k]) == 1])
  logL <- logL +
    sum(log(unlist(lapply(precisions, Matrix::det))) * (n * prod(p) / p / 2))
  ### Returns the Value
  return(logL)
}
