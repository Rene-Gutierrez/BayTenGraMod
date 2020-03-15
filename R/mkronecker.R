#' Perfroms the Kronecker Product from a list of Matrices
#'
#' This function performs the kronecker product of a list of matrices in the
#' order in which they appear on the list with the first element of the list
#' beign the rightmost matrix in the product and the last element of the list
#' beign the leftmost matrix in the product.
#'
#' @param C The list of matrices on which the Kronecker product is to be
#'           performed
#'
#' @return A matrix which is the kronecker product of the list of matrices
#'         with number of rows eaul to the product of the number of
#'         rows of every matrix, and number of columns equal to the number of
#'         columns of every matrix in the list.
#'
#' @author Rene Gutierrez Marquez

#' @export

###############################################################################
###
### Multiple Kronecker Products C[m] o C[m-1] o ... C[1]
###
### Input:
### C: List of matrices
###
### Output:
### K: Kronecker Product of the list of Matrices
###
###############################################################################

mkronecker <- function(C){
  ### Initializes the Kronecker Product
  K <- 1
  ### Loops through every matrix
  for(M in C){
    K <- kronecker(M, K)
  }
  ### Returns the Value
  return(K)
}
