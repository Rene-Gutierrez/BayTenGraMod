#' Generates random samples from the Tensor Normal Distribution
#'
#' This function generates n samples from the Tensor Normal distribution as
#' specified by its covariances matrices.
#'
#' @param n         Number of samples generated.
#' @param Tmean     Tensor mean. Must be of the same dimensions as the
#'                  specified by the list of Covariance matrices provided.
#'                  If NULL (default) a zero tensor mean is generated.
#' @param SigmaList A list of the covariances matrices.
#' @param listInd   A boolean specifying if the samples should be provided in
#'                  a list or an array. If TRUE (default) a list is provided.
#'
#' @return Either a list with the samples or an array whose last dimension is
#'         the index of the sample.
#'
#' @author Rene Gutierrez Marquez

#' @export

###############################################################################
###
### Tensor Normal Sampler
###
###############################################################################

TNormalSampler <- function(n         = 1,
                           Tmean     = NULL,
                           SigmaList,
                           listInd   = TRUE){
  ### Obtains the tensor dimension
  d <- length(SigmaList)
  ### Obtains the Dimension of Each Matrix
  p <- numeric(d)
  for(i in 1:d){
    p[i] <- ncol(SigmaList[[i]])
  }
  ### Creates the Tensor Mean if Mean is Null
  if(is.null(Tmean)){
    Tmean <- array(0, p)
  }
  ### Obtains the Cholesky Decomposition of each Covariance Matrix
  L <- list()
  for(i in 1:d){
    L[[i]] <- t(chol(SigmaList[[i]]))
  }
  ### Obtains the Cholesky Factorization of the Kronecker Product
  kroL <- Matrix::Matrix(data   = mkronecker(L),
                 sparse = TRUE)
  ### Generates the samples
  sam <- rnorm(n    = n * prod(p),
               mean = 0,
               sd   = 1)
  sam <- matrix(data = sam,
                nrow = prod(p),
                ncol = n)
  sam <- kroL %*% sam
  ### Packe the Samples as a list
  samList <- list()
  for(i in 1:n){
    samList[[i]] <- array(data = sam[,i],
                          dim  = p)
  }
  ### Checks if it returns a list or an array
  if(listInd){
    out <- samList
  } else {
    out <- sam
  }
  ### Returns the Sample
  return(out)
}
