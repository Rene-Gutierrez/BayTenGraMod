#' Computes the Conditional Statistic S_k
#'
#' This function computes the Conditional Statistic S_k for the profile
#' likelihood of the precision matrix C_k of a Tensor Normal distribution. That
#' is the quadratic form of the k-mode matricization of the sample tensor with
#' the Kronecker product of the remaining Presision matrices.
#'
#' @param tensors    List of Sample Tensors
#' @param precisions A list of the precision matrices
#' @param k          An integer specifying the precision matrix for wich the S_k
#'   statistic is desired.
#'
#' @return A sqaure matrix of the same dimension as the precision matrix k in
#'   the list with the Conditional Statistic S_k for the profile likelihood of
#'   the precision matrix C_k
#'
#' @author Rene Gutierrez Marquez
#'
#' @export

###############################################################################
###
### Computes the Conditional Statistic S_k for the profile likelihood of the
### precision matrix C_k
###
### Inputs:
### tensors: Array containing the sample tensors, last dimension is the sample
###          size
### precisions:
### k:
### Outputs:
### S_k:
###
###############################################################################

Sk <- function(tensors, precisions, k){
  numSam <- length(tensors)
  ### Computes the Kronecker Product of the precision matrices
  K <- mkronecker(precisions[-k])
  ### Initializes the COnditional Statistic
  Sk <-  0
  ### Loops through the sample
  for(i in 1:numSam){
    ### Computes the kmod Matricization of sample i
    matTensor <- kModMat(tensors[[i]], k)
    Sk        <- Sk + t(matTensor) %*% K %*% matTensor
  }
  ### Return
  return(Sk)
}
