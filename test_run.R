###############################################################################
###
### Simulated Data Generation for Tensor Graphical Models
###
###############################################################################

### Set-up
set.seed(240220)
devtools::load_all(".")
library(Matrix)
library(Tlasso)

### Parameter Settings
p <- c(10, 10, 10)    # Matrix size for each dimension
d <- length(p)           # Number of tensor dimensions
r <- c(0.10, 0.20, 0.30) # Sparity level
b <- 103                 # Degrees of Freedom for Precision Matrices
n <- 10                  # Number of Observations

### Generates the Edge Matrices
E <- list()
for(i in 1:d){ # For Every Dimension
  E[[i]] <- edgMat(matSiz = p[i],
                   spa    = r[i])
}

### Generates the Precision Matrices and Covariance Matrices
O <- list()
S <- list()
for(i in 1:d){
  O[[i]] <- Matrix(data   = cov2cor(spaPDM(b = b,
                                           E = E[[i]])),
                   sparse = TRUE)
  S[[i]] <- Matrix(data   = solve(O[[i]]),
                   sparse = TRUE)
}

### Generates the Samples
print("Creating Samples")
### Estimation Sample
### As A List
lT <- TNormalSampler(n         = n,
                     SigmaList = S)
### Cleans
gc()
### As an Array
TT <- array(data = unlist(lT),
            dim  =c(p, n))
### Cleans
gc()

### Test Sample
### As A List
testlT <- TNormalSampler(n         = n,
                         SigmaList = S)
### Cleans
gc()

### As an Array
testTT <- array(data = unlist(testlT),
                dim  =c(p, n))
### Cleans
gc()

### Runs the Algorithm
print("Performing Bayesian Tensor Gaussian Graphical Model")
timBay <- Sys.time()
BFit <- BTGM(lT,
                burnin = 0,
                nmcmc  = 10)
timBay <- difftime(time1 = Sys.time(),
                   time2 = timBay,
                   units = "secs")

### Runs TLasso
### Computes the Penalizing Parameter
CC     <- 20
lambda <- CC * sqrt(log(p) / (n * prod(p) / p))
print("Performing Tensor Lasso Model")
timFre <- Sys.time()
CFit <- Tlasso.fit(data       = TT,
                   T          = 1,
                   lambda.vec = lambda)
timFre <- difftime(time1 = Sys.time(),
                   time2 = timFre,
                   units = "secs")

BMFit <- list()
for(i in 1:d){
  BMFit[[i]] <- apply(BFit$samC[[i]], c(1,2), median)
  BMFit[[i]] <- cov2cor(BMFit[[i]])
  BMFit[[i]] <- Matrix(data   = BMFit[[i]],
                       sparse = TRUE)
}

for(i in 1:d){
  CFit[[i]] <- cov2cor(CFit[[i]])
  CFit[[i]] <- Matrix(data   = CFit[[i]],
                      sparse = TRUE)
}

### Model Lists
models      <- list()
models[[1]] <- CFit
models[[2]] <- BMFit
models[[3]] <- O

print("Computing Precision Statistics")
preSta <- modComPre(modelList = models,
                    trueModel = models[[3]])

gc()

print("Computing Likelihood")
likSta <- modComLik(modelList = models,
                    tensors   = testlT)
gc()

BMFit <- list()
for(i in 1:d){
  BMFit[[i]] <- apply(BFit$samE[[i]], c(1,2), median)
  BMFit[[i]] <- Matrix(data   = BMFit[[i]],
                       sparse = TRUE)
}

gc()

for(i in 1:3){
  if(i != 2){
    for(j in 1:d){
      models[[i]][[j]] <- models[[i]][[j]] != 0
    }
  }
}

models[[2]] <- BMFit

gc()

print("Computing Adjacency Statistics")
adjSta <- modComAdj(modelList = models,
                    trueModel = models[[3]])

gc()

timSta <- c(timFre, timBay, 0)

sta <- cbind(preSta$kroFro, preSta$kroInf, preSta$sinFro, preSta$sinInf,
             adjSta$acc, adjSta$tpr, adjSta$tnr, likSta, timSta)

sta <- data.frame(sta)
colnames(sta) <- c("Norm K-F", "Norm K-I", "Avg. Norm F", "Avg. Norm I",
                   "Accuracy", "TPR", "TNR", "Log_Lik", "Time")
row.names(sta) <- c("TLasso", "TBGGM", "Truth")

filNam <- paste0("./out/res-",p[1],"-",p[2],"-",p[3],"-",r[1],"-",r[2],"-",r[3],".txt")

print(sta)

write.table(x    = sta,
            file = filNam)

### Removes Variables
remove(b)
remove(d)
remove(i)
remove(n)
remove(p)
remove(r)
remove(CC)
remove(lambda)
remove(j)
