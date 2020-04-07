#' Bayesian Gaussian Graphical Model Sampler
#'
#' This function generates nmcmc MCMC samples for a Bayesian Gaussian Graphical
#' Model following the exchange method proposed by Wang and Li, in Efficient
#' Gaussian graphical Model determination under G-Wishart prior distributions.
#'
#' @param n      Number of observations.
#' @param S      Sample Covariance  multiplied by n.
#' @param C      Initial Precision Matrix for the MCMC.
#' @param beta   Prior probability for each edge.
#' @param bPrior Prior dedrees of freedom. Defaults to 3.
#' @param DPrior Prior location SPD Matrix. Defaults to the identity matrix.
#' @param burnin Number of sample to burn in in the MCMC.
#' @param nmcmc  Number of MCMC samples desired as output. That is without
#'   considering the burn-in period.
#' @param method Either 'E' for Exchange or 'DMH' for Double Metropolis
#'   Hastings. By default is set to 'DMH'.
#'
#' @return List containg two arrays. One for the Precison matrices and another
#'   one for the adjacency matrices.
#'   \describe{
#'     \item{samC}{An Array of Precicion matrices, in which the last dimesion
#'       goes through every Precision matrix.}
#'     \item{samE}{An Array of Adjacency matrices, in which the last dimesion
#'       goes through every adjacency matrix.}
#'  }
#'
#' @author Rene Gutierrez Marquez

#' @export

###############################################################################
###
### Bayesian Gaussian Graphical Model Sampler
###
### Inputs:
###   bPrior: Prior dedrees of freedom
###   DPrior: Prior location SPD Matrix
###   n:      Number of observations
###   S:      Sample Cross-Product
###   C:      Initial Partial Covariance Matrix
###   beta:   Prior probability for each edge
###   burnin: Number of Initial Samples to Burn-In
###   nmcmc:  Number of MCMC samples saved
###   method: Either 'E' for Exchange or 'DMH' for Double Metropolis Hastings
### Outputs:
###   Csample: Samples Partial Covariance Matrices
###   Esample: Samples Adjacency Matrices
###
###############################################################################

BGGMSampler <- function(n,
                        S,
                        C,
                        beta   = 0.5,
                        bPrior = 3,
                        DPrior = diag(ncol(C)),
                        burnin = 0,
                        nmcmc  = 1,
                        method = 'DMH'){
  ### Matrix dimension
  p <- ncol(DPrior)
  ### Posterior Parameters of the G-Wishart
  bPost <- bPrior + n
  DPost <- DPrior + S
  ### Samples Arrays
  CSample <- array(data = NA,
                   dim  = c(p, p, nmcmc))
  ESample <- CSample

  ### Initialization
  ### Initial Adjacency Matrix Guess
  E      <- abs(C) > 1e-05
  ### Number of Edges
  numEdg <- (sum(E) - p) / 2
  ### Update C
  C <- C * E

  ### Samples Arrays
  for(s in 1:(nmcmc + burnin)){ # For every burn-in and numberof samples saved
    ### Samples Off Diagonal Elements
    for(i in 1:(p - 1)){ # For every row
      for(j in (i + 1):p){ # For every column (Only the upper Triangular part)
        ### Number of Edges
        numEdg <- (sum(E) - p) / 2
        ### Proposes a New Graph E' that can differ from E only on edge (i,j)
        ### ### Computes the log-odds in favor of no edge as proposal
        w <- logH(n = n,
                  S = S,
                  C = C,
                  bPrior = bPrior,
                  DPrior = DPrior,
                  i      = i,
                  j      = j)
        #print(paste0("Node ", i, ",",j))
        w <- w + log(1 - beta) - log(beta)
        ### ### Computes the probability of proposing no edge
        ### Samples A New Auxilary C
        if(method == 'DMH'){ # By Double Metropolis Hastings Algorithm
          Caux <- updij(C     = C,
                        bPost = bPrior,
                        DPost = DPrior,
                        i     = i,
                        j     = j,
                        proij = TRUE)
        } else { # By Exchange Algorithm
          Caux <- BDgraph::rgwish(n   = 1,
                                  adj = E,
                                  b   = bPrior,
                                  D   = DPrior)
        }
        ### Estimates the log Prior Constant Ratio Using 1 Observation
        ### it always assumes the top part has no edge
        norConRat <- priorConRatio(Caux,
                                   bPrior = bPrior,
                                   DPrior = DPrior,
                                   i      = i,
                                   j      = j)
        ### Samples a new Graph
        w     <- w + norConRat
        w     <- 1 / (exp(w) + 1)
        #print(paste0("The probability of an edge between ",i," and ",j," is ",round(w * 100)))
        proij <- runif(1) < w
        ### Updates the Edge Matrix
        E[i,j] = proij
        E[j,i] = proij
        ### Updates C(i,j) and C(j,j) according to algorithm 1 part 1.c
        C <- updij(C = C,
                   bPost = bPost,
                   DPost = DPost,
                   i     = i,
                   j     = j,
                   proij = proij)
      }
    }
    ### Samples C based on the Updated E
    C <- BDgraph::rgwish(n   = 1,
                         adj = E,
                         b   = bPost,
                         D   = DPost)
    C <- round(C, 7)
    ### Saves the Samples
    if(s > burnin){
      CSample[,,s - burnin] <- C
      ESample[,,s - burnin] <- E
    }
  }

  ### Returns the Samples
  return(list(C = CSample, E = ESample))
}
