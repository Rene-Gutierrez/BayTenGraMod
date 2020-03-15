#' Computes the log-likelihood of the Tesor Normal for each model estimate.
#'
#' Computes the log-likelihood of the Tesor Normal for each model estimate
#' given a sample of tensors.
#'
#' @param modelList A list of lists. The outer list is indexed by model. The
#'   inner list is indexed by Matrix.
#' @param tensors   A list with the sample of tensors
#'
#' @return A vector with the log-likelihood for each model.
#'
#' @author Rene Gutierrez Marquez

###############################################################################
###
### Performance Statistics Graph
###
###
###############################################################################

### Number

modComLik <- function(modelList, tensors){
  ### Number of Models
  numMod <- length(modelList)
  logLik <- numeric()
  for(model in modelList){
    logLik <- c(logLik, logLikTNorm(tensors    = tensors,
                                    precisions = model))
  }
  ### Return
  return(logLik)
}
