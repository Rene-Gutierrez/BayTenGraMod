conInt <- function(precisions, data){
  ### Number of Precision Matrices
  d <- length(precisions)
  ### Dimension of the precision Matrix
  p <- numeric(d)
  for (i in 1:d) {
    p[i] <- nrow(precisions[[i]])
  }
  ### Number of Observations
  n <- length(data)
  ### Transforms the data into an Array
  testTT <- array(data = unlist(data),
                  dim  =c(p, n))
  ### For every Presicion Matrix
  tesSta <- list()
  for(k in 1:d){
    ### rho matrix (Sample Covariance of residuals, biased estimate)
    rho <- Tlasso::covres(data       = data,
                          Omega.list = precisions,
                          k          = k)
    ### Variance Correction Term
    varpi2 <- varcor(data       = data,
                     Omega.list = precisions,
                     k          = k)
    ### Bias Corrected
    bias.rho <- Tlasso::biascor(rho        = rho,
                                Omega.list = precisions,
                                k          = k)
    ### Bias Corrected test statistic
    tau.Test <- matrix(0, p[k], p[k])
    for(i in 1:(p[k] - 1)){
      for(j in (i + 1):p[k]){
        tau.Test[j,i] <- sqrt((n - 1) * prod(p[-k])) * bias_rho[i,j] / sqrt(varpi2 * rho[i,i] * rho[j,j])
      }
    }
    ### Saves the Test Statistic for every Matrix
    tesSta[[k]] <- tau.Test
  }
  ### Returns the Test Statistic
  return(tesSta)
}
