library(data.table)
library(AER)
library(parallel)
library(ggplot2)
coresAllowed <- ceiling(detectCores() / 2)

runRhoSim <- function(rho, N, periods, sims){
  #initialize y0
  print(rho)
  alpha <- rnorm(sims * N)
  e <- rnorm(sims * N)
  y0 <- .5 * alpha + e
  yMat <- matrix(y0, ncol = 1)
  
  for(i in 1:periods){
    yNext <- rho * yMat[,i] + alpha + rnorm(sims * N)
    yMat <- cbind(yMat, matrix(yNext, ncol = 1))
    
  }
  
  lst.sims <- Map(function(u,v) yMat[u:v,], seq(1,nrow(yMat),N), seq(N,nrow(yMat),N))
  
  lst.results <- mclapply(lst.sims, function(mat){
    dt.pooledOLS <- data.table(y = as.vector(mat[, 2:(periods+1)]),
                               yLag = as.vector(mat[, 1:(periods)]))
    
    fdMeans <- rowMeans(mat)
    
    rhoOLS <- coef(lm(y ~ yLag, dt.pooledOLS))[2]
    
    dt.fixedEffects <- copy(dt.pooledOLS)[, `:=`(y = y - fdMeans, yLag = yLag - fdMeans)]
    rhoFE <- coef(lm(y ~ 0 + yLag, dt.fixedEffects))[1]
    
    dt.firstDifference <- data.table(y = as.vector(mat[, 3:(periods+1)]) - as.vector(mat[, 2:(periods)]),
                                     yLag = as.vector(mat[, 2:(periods)]) - as.vector(mat[, 1:(periods-1)]),
                                     instrument = as.vector(mat[, 1:(periods-1)]))
    rhoFD <- coef(lm(y ~ 0 + yLag, dt.firstDifference))[1]
    
    rhoAH <- coef(ivreg(y ~ 0 + yLag | instrument, data = dt.firstDifference))[1]
    
    return(matrix(c(rhoOLS, rhoFE, rhoFD, rhoAH), nrow = 1))
  }, mc.cores = coresAllowed)
  
  mat.results <- Reduce(rbind, lst.results)
  
  bias <- colMeans(mat.results) - rho
  bias
  
  se1 <- apply(mat.results, 2, sd)
  se1
  
  rmse <- (bias^2 + se1^2)^(1/2)
  rmse
  
  lst.output <- list(estimators = mat.results, bias = bias, se = se1, rmse = rmse)
}

#rhoHalf <- runRhoSim(.5, 100, 6, 1000)

rhoSims <- lapply(seq(0, 1, by= .1), runRhoSim, N = 100, periods = 6, sims = 1000)

dt.bias <- data.table(trueRho = numeric(), rhoOLS = numeric(), rhoFE = numeric(), rhoFD = numeric(), rhoAH = numeric())
dt.se <- data.table(trueRho = numeric(), rhoOLS = numeric(), rhoFE = numeric(), rhoFD = numeric(), rhoAH = numeric())
dt.rmse <- data.table(trueRho = numeric(), rhoOLS = numeric(), rhoFE = numeric(), rhoFD = numeric(), rhoAH = numeric())
rhoSeq <- seq(0, 1, by = .1)

for(i in 1:length(seq(0, 1, by = .1))){
  dt.bias <- rbind(dt.bias, matrix(c(rhoSeq[i], rhoSims[[i]]$bias), nrow = 1), use.names = F)
  dt.se <- rbind(dt.se, matrix(c(rhoSeq[i], rhoSims[[i]]$se), nrow = 1), use.names = F)
  dt.rmse <- rbind(dt.rmse, matrix(c(rhoSeq[i], rhoSims[[i]]$rmse), nrow = 1), use.names = F)
  
}

ggplot(melt(dt.bias, id.vars = "trueRho"), aes(x = trueRho, y = value, color = variable)) + geom_line() + 
  ylab("Value of Bias as Function of Rho") + xlab("True Rho")

ggplot(melt(dt.se, id.vars = "trueRho"), aes(x = trueRho, y = value, color = variable)) + geom_line() + 
  ylab("Value of SE as Function of Rho") + xlab("True Rho")

ggplot(melt(dt.rmse, id.vars = "trueRho"), aes(x = trueRho, y = value, color = variable)) + geom_line() + 
  ylab("Value of RMSE as Function of Rho") + xlab("True Rho")

rho7 <- rhoSims[[8]]$estimators
ggplot(data = data.table(rho7), aes(x = V4)) + geom_histogram(aes(y = ..density..)) + geom_line(aes(x = seq(from = -5, to = 5, length.out = 1000), y = dnorm(x = seq(from = -5, to = 5, length.out = 1000), mean = .7, sd = sd(rho7[,4]))), color = "red") +  
  ylab("Simulated Density") + xlab("AH Estimator") + coord_cartesian(xlim = c(-5, 5))

ggplot(data = data.table(rho7), aes(x = V3)) + geom_histogram(aes(y = ..density..)) + geom_line(aes(x = seq(from = -1, to = 1, length.out = 1000), y = dnorm(x = seq(from = -1, to = 1, length.out = 1000), mean = .7, sd = sd(rho7[,3]))), color = "red") +  
  ylab("Simulated Density") + xlab("FD Estimator") + coord_cartesian(xlim = c(-1, 1))

ggplot(data = data.table(rho7), aes(x = V2)) + geom_histogram(aes(y = ..density..)) + geom_line(aes(x = seq(from = -0, to = 1, length.out = 1000), y = dnorm(x = seq(from = -0, to = 1, length.out = 1000), mean = .7, sd = sd(rho7[,2]))), color = "red") +  
  ylab("Simulated Density") + xlab("FE Estimator") + coord_cartesian(xlim = c(-0, 1))

ggplot(data = data.table(rho7), aes(x = V1)) + geom_histogram(aes(y = ..density..)) + geom_line(aes(x = seq(from = .5, to = 1.5, length.out = 1000), y = dnorm(x = seq(from = .5, to = 1.5, length.out = 1000), mean = .7, sd = sd(rho7[,1]))), color = "red") +  
  ylab("Simulated Density") + xlab("POLS Estimator") + coord_cartesian(xlim = c(.5, 1.5))

#Now we do 3 time periods
rhoSims <- lapply(seq(0, 1, by= .1), runRhoSim, N = 100, periods = 3, sims = 1000)


dt.bias <- data.table(trueRho = numeric(), rhoOLS = numeric(), rhoFE = numeric(), rhoFD = numeric(), rhoAH = numeric())
dt.se <- data.table(trueRho = numeric(), rhoOLS = numeric(), rhoFE = numeric(), rhoFD = numeric(), rhoAH = numeric())
dt.rmse <- data.table(trueRho = numeric(), rhoOLS = numeric(), rhoFE = numeric(), rhoFD = numeric(), rhoAH = numeric())
rhoSeq <- seq(0, 1, by = .1)

for(i in 1:length(seq(0, 1, by = .1))){
  dt.bias <- rbind(dt.bias, matrix(c(rhoSeq[i], rhoSims[[i]]$bias), nrow = 1), use.names = F)
  dt.se <- rbind(dt.se, matrix(c(rhoSeq[i], rhoSims[[i]]$se), nrow = 1), use.names = F)
  dt.rmse <- rbind(dt.rmse, matrix(c(rhoSeq[i], rhoSims[[i]]$rmse), nrow = 1), use.names = F)
  
}

ggplot(melt(dt.bias, id.vars = "trueRho"), aes(x = trueRho, y = value, color = variable)) + geom_line() + 
  ylab("Value of Bias as Function of Rho") + xlab("True Rho")

ggplot(melt(dt.se, id.vars = "trueRho"), aes(x = trueRho, y = value, color = variable)) + geom_line() + 
  ylab("Value of SE as Function of Rho") + xlab("True Rho")

ggplot(melt(dt.rmse, id.vars = "trueRho"), aes(x = trueRho, y = value, color = variable)) + geom_line() + 
  ylab("Value of RMSE as Function of Rho") + xlab("True Rho")

rho7 <- rhoSims[[8]]$estimators
ggplot(data = data.table(rho7), aes(x = V4)) + geom_histogram(aes(y = ..density..)) + geom_line(aes(x = seq(from = -5, to = 5, length.out = 1000), y = dnorm(x = seq(from = -5, to = 5, length.out = 1000), mean = .7, sd = sd(rho7[,4]))), color = "red") +  
  ylab("Simulated Density") + xlab("AH Estimator") + coord_cartesian(xlim = c(-5, 5))

ggplot(data = data.table(rho7), aes(x = V3)) + geom_histogram(aes(y = ..density..)) + geom_line(aes(x = seq(from = -1, to = 1, length.out = 1000), y = dnorm(x = seq(from = -1, to = 1, length.out = 1000), mean = .7, sd = sd(rho7[,3]))), color = "red") +  
  ylab("Simulated Density") + xlab("FD Estimator") + coord_cartesian(xlim = c(-1, 1))

ggplot(data = data.table(rho7), aes(x = V2)) + geom_histogram(aes(y = ..density..)) + geom_line(aes(x = seq(from = -0, to = 1, length.out = 1000), y = dnorm(x = seq(from = -0, to = 1, length.out = 1000), mean = .7, sd = sd(rho7[,2]))), color = "red") +  
  ylab("Simulated Density") + xlab("FE Estimator") + coord_cartesian(xlim = c(-0, 1))

ggplot(data = data.table(rho7), aes(x = V1)) + geom_histogram(aes(y = ..density..)) + geom_line(aes(x = seq(from = .5, to = 1.5, length.out = 1000), y = dnorm(x = seq(from = .5, to = 1.5, length.out = 1000), mean = .7, sd = sd(rho7[,1]))), color = "red") +  
  ylab("Simulated Density") + xlab("POLS Estimator") + coord_cartesian(xlim = c(.5, 1.5))

#Now we do 9 time periods
rhoSims <- lapply(seq(0, 1, by= .1), runRhoSim, N = 100, periods = 9, sims = 1000)


dt.bias <- data.table(trueRho = numeric(), rhoOLS = numeric(), rhoFE = numeric(), rhoFD = numeric(), rhoAH = numeric())
dt.se <- data.table(trueRho = numeric(), rhoOLS = numeric(), rhoFE = numeric(), rhoFD = numeric(), rhoAH = numeric())
dt.rmse <- data.table(trueRho = numeric(), rhoOLS = numeric(), rhoFE = numeric(), rhoFD = numeric(), rhoAH = numeric())
rhoSeq <- seq(0, 1, by = .1)

for(i in 1:length(seq(0, 1, by = .1))){
  dt.bias <- rbind(dt.bias, matrix(c(rhoSeq[i], rhoSims[[i]]$bias), nrow = 1), use.names = F)
  dt.se <- rbind(dt.se, matrix(c(rhoSeq[i], rhoSims[[i]]$se), nrow = 1), use.names = F)
  dt.rmse <- rbind(dt.rmse, matrix(c(rhoSeq[i], rhoSims[[i]]$rmse), nrow = 1), use.names = F)
  
}

ggplot(melt(dt.bias, id.vars = "trueRho"), aes(x = trueRho, y = value, color = variable)) + geom_line() + 
  ylab("Value of Bias as Function of Rho") + xlab("True Rho")

ggplot(melt(dt.se, id.vars = "trueRho"), aes(x = trueRho, y = value, color = variable)) + geom_line() + 
  ylab("Value of SE as Function of Rho") + xlab("True Rho")

ggplot(melt(dt.rmse, id.vars = "trueRho"), aes(x = trueRho, y = value, color = variable)) + geom_line() + 
  ylab("Value of RMSE as Function of Rho") + xlab("True Rho")

rho7 <- rhoSims[[8]]$estimators
ggplot(data = data.table(rho7), aes(x = V4)) + geom_histogram(aes(y = ..density..)) + geom_line(aes(x = seq(from = -5, to = 5, length.out = 1000), y = dnorm(x = seq(from = -5, to = 5, length.out = 1000), mean = .7, sd = sd(rho7[,4]))), color = "red") +  
  ylab("Simulated Density") + xlab("AH Estimator") + coord_cartesian(xlim = c(-5, 5))

ggplot(data = data.table(rho7), aes(x = V3)) + geom_histogram(aes(y = ..density..)) + geom_line(aes(x = seq(from = -1, to = 1, length.out = 1000), y = dnorm(x = seq(from = -1, to = 1, length.out = 1000), mean = .7, sd = sd(rho7[,3]))), color = "red") +  
  ylab("Simulated Density") + xlab("FD Estimator") + coord_cartesian(xlim = c(-1, 1))

ggplot(data = data.table(rho7), aes(x = V2)) + geom_histogram(aes(y = ..density..)) + geom_line(aes(x = seq(from = -0, to = 1, length.out = 1000), y = dnorm(x = seq(from = -0, to = 1, length.out = 1000), mean = .7, sd = sd(rho7[,2]))), color = "red") +  
  ylab("Simulated Density") + xlab("FE Estimator") + coord_cartesian(xlim = c(-0, 1))

ggplot(data = data.table(rho7), aes(x = V1)) + geom_histogram(aes(y = ..density..)) + geom_line(aes(x = seq(from = .5, to = 1.5, length.out = 1000), y = dnorm(x = seq(from = .5, to = 1.5, length.out = 1000), mean = .7, sd = sd(rho7[,1]))), color = "red") +  
  ylab("Simulated Density") + xlab("POLS Estimator") + coord_cartesian(xlim = c(.5, 1.5))
