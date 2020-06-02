library(data.table)
library(mlogit)
library(survival)

dt.data <- setDT(read.csv("~/Downloads/ca_heating.csv"))
dt.data
dt.dataLong <- melt(dt.data, id.vars = c("depvar", "idcase"))

dt.dataLong <- dt.dataLong[variable %like% "ic" | variable %like% "oc"]
dt.dataLong[, `:=`(cost = substr(variable,1,2), heating = as.numeric(substr(variable, 3, 3)))]
dt.dataLong[, `:=`(variable = NULL, case = (heating == depvar), depvar = NULL)]
dt.dataLong <- dcast(dt.dataLong, idcase + case + heating ~ cost)
dt.dataLong[, heating := as.factor(heating)]

firstModel <- clogit(case ~ ic + oc + strata(idcase), dt.dataLong)
coefModel <- coef(firstModel)

dt.prediction <- copy(dt.dataLong)
dt.prediction[, prob := exp(coefModel[1] * ic + coefModel[2] * oc)]
dt.prediction[, prob := prob / sum(prob), by = idcase]
dt.prediction[, mean(prob), keyby = heating]
table(dt.data$depvar) / sum(table(dt.data$depvar))

wtp <- coefModel[2] / coefModel[1]
r <- 1/wtp

secondModel <- clogit(case ~ heating + ic + oc + strata(idcase), dt.dataLong)
secondModel
coefModel <- c(0,coef(secondModel))


dt.prediction <- copy(dt.dataLong)
dt.prediction[, prob := exp(coefModel[6] * ic + coefModel[7] * oc + coefModel[as.numeric(heating)])]
dt.prediction[, prob := prob / sum(prob), by = idcase]
dt.prediction[, mean(prob), keyby = heating]
table(dt.data$depvar) / sum(table(dt.data$depvar))

wtp2 <- coefModel[7] / coefModel[6]
r2 <- 1/wtp2

dt.prediction[, probnic := exp(coefModel[6] * ic * ifelse(as.numeric(heating) == 5, .9, 1) + coefModel[7] * oc + coefModel[as.numeric(heating)])]
dt.prediction[, probnic := probnic / sum(probnic), by = idcase]
dt.prediction[, .(mean(prob), mean(probnic)), keyby = heating]

coefModel <- c(0,coef(secondModel)[1:4], coef(secondModel)[2], coef(secondModel)[5:6])


dt.prediction <- copy(dt.dataLong)
dt.new <- dt.prediction[heating == 3]
dt.new[, `:=`(heating = NULL, ic = ic + 200, oc = .75 * oc)]
dt.new[, heating := 6]
dt.prediction[, heating := as.numeric(heating)]
dt.prediction <- rbind(dt.prediction, dt.new)
dt.prediction[, prob := exp(coefModel[7] * ic + coefModel[8] * oc + coefModel[as.numeric(heating)])]
dt.prediction[, prob := prob / sum(prob), by = idcase]
dt.prediction[, mean(prob), keyby = heating]
table(dt.data$depvar) / sum(table(dt.data$depvar))

