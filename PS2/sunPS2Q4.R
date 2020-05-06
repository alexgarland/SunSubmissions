library(readxl)
library(data.table)
library(plm)
library(Matrix)
library(lmtest)
library(dynpanel)

dt.data <- data.table(read_xls("cigar.xls"))
colnames(dt.data) <- c("state", "year", "lnc", "lnp", "lni", "lnm")
dt.data[, lagc := c(NA, dt.data[1:(NROW(dt.data)-1)]$lnc)]
dt.data

#Pooling data
plm.pooled <- plm(lnc ~ lagc + lnp + lni + lnm, model = "pooling", data = dt.data)
coeftest(plm.pooled, vcov = vcovHC(plm.pooled, "white1"))

plm.yearFE <- plm(lnc ~ lagc + lnp + lni + lnm + factor(year), model = "pooling", data = dt.data)
coeftest(plm.yearFE, vcov = vcovHC(plm.yearFE, "white1"))


plm.within <- plm(lnc ~ lagc + lnp + lni + lnm, model = "within", data = dt.data, index = "state")
coeftest(plm.within, vcov = vcovHC(plm.within, "white1"))


plm.withinYear <- plm(lnc ~ lagc + lnp + lni + lnm + factor(year), model = "within", data = dt.data, index = "state")
coeftest(plm.withinYear, vcov = vcovHC(plm.withinYear, "white1"))
waldtest(plm.withinYear, plm.within, vcov = vcovHC(plm.withinYear, "white1"))


#AH Estimator with pgmm
ahEstimator <- pgmm(lnc ~ lag(lnc, 1) + lnp + lni + lnm | lag(lnc, 2), data = dt.data, index = c("state", "year"), effect = c("twoways"), fsm = "I", collapse = T)
summary(ahEstimator, robust=T)

#AB Estimator with pgmm
abEstimatorAll <- pgmm(lnc ~ lag(lnc, 1) + lnp + lni + lnm | lag(lnc, 2:99), data = dt.data, index = c("state", "year"), effect = c("twoways"), model = c("onestep"), fsm = "G")
abEstimatorFour <- pgmm(lnc ~ lag(lnc, 1) + lnp + lni + lnm | lag(lnc, 2:4), data = dt.data, index = c("state", "year"), effect = c("twoways"), model = c("onestep"), fsm = "G")
abEstimator <- pgmm(lnc ~ lag(lnc, 1) + lnp + lni + lnm | lag(lnc, 2), data = dt.data, index = c("state", "year"), effect = c("twoways"), fsm = "G")
summary(abEstimatorAll, robust= T)
summary(abEstimatorFour, robust = T)
summary(abEstimator, robust= T)
