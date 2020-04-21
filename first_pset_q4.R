library(MASS)
library(parallel)
coresAllowed <- ceiling(detectCores() / 2)
n = 500
t = 5
explanSize = 1
set.seed(2020)

genRawErrors <- function(rawDataRow){
  #vcov <- matrix(abs(rawDataRow), nrow = 1) %*% abs(rawDataRow)
  rnorm(1, sd = sum(rawDataRow^2))
}

calcPersonMeans <- function(rawData){
  personMeans <- apply(rawData, c(1,3), mean)
}

calcPersonDeviations <- function(data, personMeans){
  data - personMeans
}

calcSXX <- function(rawData2){
  personExplMeans <- calcPersonMeans(rawData2)
  personXDeviations <- rawData2 - aperm(array(personExplMeans, dim = c(dim(rawData2)[1], dim(rawData2)[3],dim(rawData2)[2])), c(1, 3,2))
  Reduce(`+`, lapply(1:dim(rawData2)[1], function(i){
    Reduce(`+`,lapply(1:dim(rawData2)[2], function(j){
      matrix(personXDeviations[i,j,], ncol = 1) %*% personXDeviations[i,j,]
    }))
  }))
}

calcTimeCov <- function(singleXDev, singleErrDev){
  if(is.vector(singleXDev)){
    return(sum(singleXDev * singleErrDev))
  }
  return(colSums(singleXDev * singleErrDev))
}
  
calcPersonXU <- function(rawData2, Errors2){
  personExplMeans <- calcPersonMeans(rawData2)
  personErrMeans <- calcPersonMeans(Errors2)
  personErrDeviations <- Errors2 - array(personErrMeans, dim = c(dim(rawData2)[1], dim(rawData2)[2], 1))
  personXDeviations <- rawData2 - aperm(array(personExplMeans, dim = c(dim(rawData2)[1], dim(rawData2)[3],dim(rawData2)[2])), c(1, 3,2))
  personXU <- t(sapply(1:dim(rawData2)[1], function(i){
    return(calcTimeCov(personXDeviations[i,,], personErrDeviations[i,,]))
  }))
  if(dim(personXU)[1] == 1 & dim(rawData2)[1] != 1){
    personXU <- aperm(personXU, c(2, 1))
  }
  return(personXU)
}

calcSigHatSq <- function(SXXInv, XU){
  perPerson <- lapply(1:NROW(XU), function(i){
    matrix(XU[i,], ncol = 1) %*% XU[i,]
  })
  SXXInv %*% SXXInv %*% Reduce(`+`, perPerson)
}

calcSigBarSq <- function(SXXInv, rawData2, errorArray){
  personExplMeans <- calcPersonMeans(rawData2)
  personXDeviations <- rawData2 - aperm(array(personExplMeans, dim = c(dim(rawData2)[1], dim(rawData2)[3],dim(rawData2)[2])), c(1, 3,2))
  errorMeans <- calcPersonMeans(errorArray)
  summandTerm <- Reduce(`+`, lapply(1:dim(rawData2)[1], function(i){
    Reduce(`+`,lapply(1:dim(rawData2)[2], function(j){
      matrix(personXDeviations[i,j,], ncol = 1) %*% personXDeviations[i,j,] * (errorArray[i, j, 1] - errorMeans[i,1])^2
    }))
  }))
  
  SXXInv %*% SXXInv %*% summandTerm
}

calcSyntheticErrors <- function(personRawData, betaGuess, personRawY){
  array(personRawY - personRawData %*% betaGuess, dim = c(1, t, 1))
}

runSimulation <- function(i){
  rawData <- mvrnorm(n * t, mu = rep(0, explanSize), Sigma = diag(explanSize))
  rawErrors <- t(apply(rawData, 1, genRawErrors))
  #rawY <- t(rowSums(rawData) + rawErrors)
  
  xData <- array(rawData, dim = c(n, t, explanSize))
  trueErrors <- array(rawErrors, dim = c(n, t, 1))
  y <- array(apply(xData, c(1,2), sum) + trueErrors[,,1], dim = c(n, t, 1))
  trueXU <- calcPersonXU(xData, trueErrors)
  SXX <- calcSXX(xData)
  SXXInv <- solve(SXX)
  
  betaGuess <- SXXInv %*% colSums(trueXU) + rep(1, explanSize)
  betaGuess
  
  syntheticErrors <- array(do.call(rbind, lapply(1:n, function(i){
    calcSyntheticErrors(xData[i,,], betaGuess, y[i,,])
  })), dim = c(n, t, 1))
  
  syntheticXU <- calcPersonXU(xData, syntheticErrors)
  
  sigHatSq <- calcSigHatSq(SXXInv, syntheticXU)
  sigBarSq <- calcSigBarSq(SXXInv, xData, syntheticErrors)
  return(list(beta = betaGuess, sigHat = sqrt(diag(sigHatSq)), sigBar = sqrt(diag(sigBarSq))))
}

a <- Sys.time()
simRuns <- mclapply(1:1000, runSimulation, mc.cores = coresAllowed)
print(Sys.time() - a)

betaMat <- do.call(rbind, lapply(1:1000, function(i){
  return(t(simRuns[[i]]$beta))
}))

sigHatMat <- do.call(rbind, lapply(1:1000, function(i){
  return(simRuns[[i]]$sigHat)
}))

sigBarMat <- do.call(rbind, lapply(1:1000, function(i){
  return(simRuns[[i]]$sigBar)
}))

betaMatAvg <- colMeans(betaMat)
sdBeta <- sqrt(diag(var(betaMat)))

sqrt(mean((sigHatMat - sdBeta)^2))
mean(sigHatMat) - sdBeta

sqrt(mean((sigBarMat - sdBeta)^2))
mean(sigBarMat) - sdBeta

sd(sigHatMat)
sd(sigBarMat)
