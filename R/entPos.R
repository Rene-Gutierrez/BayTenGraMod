#' Computes Graphic outputs of the Marginal Posterior of the Precision Matrices
#'
#' Computes the trace plot and density plot of one entry of a Presicion Matrix.
#'
#' @param posPre A list of arrays for Precicion matrices, the list goes through
#'   every precision matrix.
#' @param i      Row number of the entry.
#' @param j      Column number of the entry
#'
#' @return A data frame used to plot the trace plot and density of an entry. It
#'   specifies the matrix is comming from and the entry. Test
#'
#' @author Rene Gutierrez Marquez
#'
#' @export

entPos <- function(posPre,
                   truth  = NULL,
                   matNum = 1,
                   i      = 1,
                   j      = 2){
  ### Gets the number of iterations of the MCMC
  n   <- dim(posPre)[3]
  dat <- data.frame(sample = numeric(length = n),
                    truth  = numeric(length = n),
                    iter   = integer(length = n),
                    entry  = character(length = n),
                    stringsAsFactors = FALSE)
  ### Obtains the Sample (Correlation)
  dat$sample <- posPre[i, j, ] / sqrt(posPre[i, i, ] * posPre[j, j, ])
  ### Creates a Data Frame
  dat$truth  <- truth
  dat$iter   <- 1:n
  dat$entry  <- paste0(matNum,",",i,",",j)
  ### Returns the data frame
  return(dat)
}
