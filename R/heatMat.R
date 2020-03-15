#' Creates a Heatmap Database given a Matrix
#'
#' This function creates a data frame to create Heatmaps for each model and for
#' each Matrix
#'
#' @param ll         A list of lists. The outer list is indexed by model. The
#'   inner list is indexed by matrix.
#' @param modelNames The names of the models to use in the data Frame.
#'
#' @return A data Frame to create heatmaps in ggplot. With columns indicating
#'   the Matrix and Model
#'
#' @author Rene Gutierrez Marquez

###############################################################################
###
### Plots the Compares Precision Matrices with a Heat Map
###
### Inputs:
### ll:         A list of list of Precion Matrices
###             (possibly including the Truth) or a list of the
###             Adjacency Matrix
### modelNames: A list of names for the models (or True Parameters)
###
###############################################################################

heatMat <- function(ll,
                    modelNames = NULL){
  ### Number of Candidates
  n <- length(ll)
  ### If model names is Null creates a list of Names for Each Model
  if(is.null(modelNames)){
    modelNames <- character()
    for(i in 1:n){
      modelNames <- c(modelNames, paste0("Model ",i))
    }
  }
  ### Gets the Number of Precision Matrices to be evaluated
  numMat <- length(ll[[1]])
  ### Creates the Matrix Names
  matNam <- character()
  for(i in 1:numMat){
    matNam <- c(matNam, paste0("Matrix ", i))
  }
  ### Creates a list of Data Frames with the respecting
  dataList <- list()
  ### Initializes the Data Frame
  tab <- data.frame(x      = numeric(),
                    y      = numeric(),
                    value  = numeric(),
                    model  = character(),
                    matrix = character())
  for(i in 1:numMat){
    ### Obtains the Matrix Size
    p <- ncol(ll[[1]][[i]])
    ### Creates the Grid
    x <- seq(1,p)
    y <- seq(p,1)
    ### For every Model
    for(j in 1:n){
      ### Creates a Table for Model i
      tem           <- cbind(expand.grid(x = x,
                                         y = y),
                             as.vector(ll[[j]][[i]]),
                             rep(x     = modelNames[j],
                                 times = p),
                             rep(x     = matNam[i],
                                 times = p))
      colnames(tem) <- c("x", "y", "value", "model", "matrix")
      ### Joins the tables
      tab <- rbind(tab, tem)
    }
  }
  ### Returns Something?
  return(tab)
}
