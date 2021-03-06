#' Performs Bayesian Tensor Graphical Model Estimation
#'
#' This function generates nmcmc MCMC samples for a Bayesian Tensor Graphical
#' Model. Each precision matrix is modeled according to a GWishart distribution.
#' The output, returns the samples for each presion matrix as well as the
#' adjacency matrix samples.
#'
#' @param t      List of Sample Tensors
#' @param b      A vector with the Prior degrees of Freedom for each precision
#'   matrix. If NULL it sets the prior for each presicion at 3.
#' @param D      List of Prior Scale matrices for each precision matrix. If NULL
#'   it defaults to the corresponding Identity matrices for each presicion
#'   matrix.
#' @param C      List of Initial Values of the Precicion Matrices to start the
#'   MC markov chain. If NULL each presicion matrix is initialized as the
#'   corresponding identity matrix.
#' @param beta   A vector containing the prior probability of having and edge
#'   between to vertices for each precision. If NULL it sets the value at 0.5
#'   for every edge matrix.
#' @param burnin Number of sample to burn in in the MCMC.
#' @param nmcmc  Number of MCMC samples desired as output. That is without
#'   considering the burn-in period.
#' @param method Either 'E' for Exchange or 'DMH' for Double Metropolis
#'   Hastings. By default is set to 'DMH'.
#' @param norInd. 2-dimensional index to normalize for Identifiability. The
#'   first entry indicates the Precision Matrix, while the second entry
#'   indicates the diagonal element. By default it normilizes the first
#'   diagonal element of the first precision matrix.
#' @param progress Boolean indicating the use of the progress bar. By default
#'   it is set true.
#'
#' @return List containg two other lists. One for the Precison matrices and
#'   another one for the adjacency matrices.
#' \describe{
#'   \item{samC}{A list of arrays for Precicion matrices, the list goes through
#'   every precision matrix.}
#'   \item{samE}{A list of arrays for Adjacency matrices, the list goes through
#'   every adjacency matrix.}
#' }
#'
#' @author Rene Gutierrez Marquez

#' @export

###############################################################################
###
### Bayesian Tensor Lasso
###
### Inputs:
### t:      List of Sample Tensors
### b:      Vector of prior degrees of freedom parameters
### D:      List of Prior Scale Matrices
### C:      List of Initial Presicion Matrices
### beta:   Vector of edge inclusion probability
### burnin: Burn-in for MCMC
### nmcmc:  Number of MCMC samples in the output
###
### Outputs:
### C: List of sample Precisions
### E: LIst of sample Adjacency Matrices
###
###############################################################################

BTGM <- function(t,
                 b        = NULL,
                 D        = NULL,
                 C        = NULL,
                 beta     = NULL,
                 burnin   = 0,
                 nmcmc    = 1,
                 method   = 'DMH',
                 norInd   = c(1, 1),
                 progress = TRUE){
  ### Computes some Global Parameters
  ### Number of Samples
  numSam <- length(t)
  ### Dimensions of the Tensors
  p      <- dim(t[[1]])
  ### Number of dimensions of the Tensors
  numDim <- length(p)
  ### Sets Values for the Null Entries
  ### For b
  if(is.null(b)){
    b <- rep(3, numDim)
  }
  ### For D
  if(is.null(D)){
    D <- list()
    for(i in 1:numDim){
      D[[i]] <- diag(p[i])
    }
  }
  ### For C
  if(is.null(C)){
    C <- list()
    for(i in 1:numDim){
      C[[i]] <- diag(p[i])
    }
  }
  ### For beta
  if(is.null(beta)){
    beta <- rep(0.5, numDim)
  }
  ### Begins the Gibbs Sampler
  samC <- list()
  for(i in 1:numDim){
    samC[[i]] <- array(data = NA,
                       dim  = c(p[i], p[i], nmcmc))
  }
  samE <- list()
  for(i in 1:numDim){
    samE[[i]] <- array(data = NA,
                       dim  = c(p[i], p[i], nmcmc))
  }
  ### Progress Bar
  if(progress){
    pb <- txtProgressBar(min     = 0,
                         max     = 1,
                         initial = 0,
                         style   = 3,
                         width   = 72)
  }
  ### For Every Gibbs Iteration
  for(s in 1:(burnin + nmcmc)){
    ### For every tensor dimension
    for(j in 1:numDim){
      ### Computes some Auxilary Variables
      bPrior <- b[j]
      DPrior <- D[[j]]
      Cj     <- C[[j]]
      n      <- numSam * prod(p) / p[j]
      sS     <- Sk(tensors    = t,
                   precisions = C,
                   k          = j)
      ### Computes the New Sample for Precision j
      outCj <- BGGMSampler(n      = n,
                           S      = sS,
                           C      = Cj,
                           beta   = 0.5,
                           bPrior = bPrior,
                           DPrior = DPrior,
                           burnin = 0,
                           nmcmc  = 1)
      ### Updates the Initial Point for the Sampler
      C[[j]] <- outCj$C[,,1]
      ### Saves the Samples
      if(s > burnin){
        samC[[j]][,,s - burnin] <- C[[j]]
        samE[[j]][,,s - burnin] <- outCj$E[,,1]
      }
    }
    ### Progress Bar Display
    if(progress){
      setTxtProgressBar(pb    = pb,
                        value = s / (burnin + nmcmc))
    }
  }
  ### Closes the Progress Bar
  if(progress){
    close(pb)
  }
  return(list(samC = samC, samE = samE))
}
