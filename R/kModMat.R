#' k-mode Tesnor Matricization
#'
#' This function computes the k-mode Tesnor Matricization of a tensor.
#'
#' @param tensor An array for the Tensor to matricize.
#' @param k      An integer specifying the mode to matricize.
#'
#' @return A matrix that is the k-mode matricization of the tensor provided.
#'
#' @author Rene Gutierrez Marquez

###############################################################################
###
### k_mode Tesnor Matricization
###
### Inputs:
### tensor: Tensor to be Matricized
### mode:   Mode to perform the Matricization
###
### Output:
### matTen: The matricized tensor
###
###############################################################################

kModMat <- function(tensor, mode = 1){
  ### Length of Mode k
  p <- dim(tensor)[mode]
  ### Number of Rows of the Matricized Tensor
  n <- prod(dim(tensor)) / p
  ### Createst the Matrizization
  matTen = matrix(NA, n, p)
  for(i in 1:p){ # Loops through every "row" of the mode
    matTen[,i] <- tensor[slice.index(tensor, mode) == i]
  }
  ### Output
  return(matTen)
}
