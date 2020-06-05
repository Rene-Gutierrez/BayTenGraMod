### Creates Precision Matrices

preMatCre <- function(p        = 10,
                      type     = 'random',
                      ar       = 1,
                      sparsity = 0.1){
  ### Creates the Matrix
  O = matrix(data = 0,
             nrow = p,
             ncol = p)
  ### Sets the Diagonal
  diag(O) <- rep(1, p)
  if(type == "circular"){
    ### Sets the secon diagonal
    for(i in 2:p){
      O[i - 1, i] <- 0.5
      O[i, i - 1] <- 0.5
    }
    ### Sets the Corner
    O[1, p] <- 0.45
    O[p, 1] <- 0.45
  } else if(type == "ar"){
    ### Sets the diagonals
    for(j in 2:(ar + 1)){
      for(i in j:p){
        O[i - j + 1, i] <- (1 / 2)^(j - 1)
        O[i, i - j + 1] <- (1 / 2)^(j - 1)
      }
    }
  } else if(type == "star"){
    ### Sets the edges
    for(i in 2:p){
      O[i, 1] <- 0.1
      O[1, i] <- 0.1
    }
  }
  ### Creates the Edge Matrix
  E = O
  E[O != 0] = 1
  ### Returns the Graph
  return(list(P = O, S = solve(O), E = E))
}
