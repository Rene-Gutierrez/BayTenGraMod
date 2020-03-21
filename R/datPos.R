#' Computes Graphic outputs of the Marginal Posterior of the Precision Matrices
#'
#' Computes the trace plot and density plot of one entry of a Presicion Matrix.
#'
#' @param posPre A list of arrays for Precicion matrices, the list goes through
#'   every precision matrix.
#' @param reaPre A list of the real precision matrices. By default it is set to
#'   NULL, in which case no real values are ploted for reference.
#' @param sel    A list of selected entries for each Matrix
#'
#' @return A data frame used to plot the trace plot and density of each entry.
#'  It specifies the matrix is comming from and the entry.
#'
#' @author Rene Gutierrez Marquez
#'
#' @export

datPos <- function(posPre,
                   reaPre = NULL,
                   sel){
  ### Number of Matrices
  d <- length(reaPre)
  ### Initialices the Data Frame
  dat <- data.frame(sample = numeric(length = 0),
                    truth  = numeric(length = 0),
                    iter   = integer(length = 0),
                    entry  = character(length = 0),
                    stringsAsFactors = FALSE)
  ### For Each Matrix
  for (i in 1:d) {
    ### Obtains the matrix size
    p <- dim(posPre[[i]])[1]
    ### For Each Entry
    for (entry in sel[[i]]) {
      ### Obtains the row and column number
      rowInd <- entry %/% p + 1
      colInd <- entry %%  p
      if(colInd == 0){
        colInd = p
      }
      ### Computes a temporary data frame
      tem <- entPos(posPre = posPre[[i]],
                    truth  = reaPre[[i]][rowInd, colInd],
                    matNum = i,
                    i      = rowInd,
                    j      = colInd)
      ### Updates the data Frame
      dat <- rbind(dat, tem)
    }
  }
  ### Returns the data frame
  return(dat)
}
