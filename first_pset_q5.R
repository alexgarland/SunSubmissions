library(foreign)
library(data.table)
library(sandwich)
library(lmtest)
library(plm)
dt <- data.table(read.dta("handguns.dta"))

#Part I
lm.simpleVio <- lm(log(vio) ~ shall, data = dt)
coeftest(lm.simpleVio, vcov = vcovHC(lm.simpleVio, "HC1"))
summary(lm.simpleVio)$r.squared

lm.simpleMur <- lm(log(mur) ~ shall, data = dt)
coeftest(lm.simpleMur, vcov = vcovHC(lm.simpleMur, "HC1"))
summary(lm.simpleMur)$r.squared

lm.simpleRob <- lm(log(rob) ~ shall, data = dt)
coeftest(lm.simpleRob, vcov = vcovHC(lm.simpleRob, "HC1"))
summary(lm.simpleRob)$r.squared

#Part II
lm.Vio <- lm(log(vio) ~ shall + incarc_rate + density + pop + pm1029 + avginc, data = dt)
coeftest(lm.Vio, vcov = vcovHC(lm.Vio, "HC1"))
summary(lm.Vio)$r.squared

lm.Mur <- lm(log(mur) ~ shall + incarc_rate + density + pop + pm1029 + avginc, data = dt)
coeftest(lm.Mur, vcov = vcovHC(lm.Mur, "HC1"))
summary(lm.Mur)$r.squared

lm.Rob <- lm(log(rob) ~ shall + incarc_rate + density + pop + pm1029 + avginc, data = dt)
coeftest(lm.Rob, vcov = vcovHC(lm.Rob, "HC1"))
summary(lm.Rob)$r.squared

#Part IV
var <- 'rob'

for (var in c("vio", "mur", "rob")){
  print(var)
  lm.basic <- lm(log(dt[[eval(var)]]) ~ shall + incarc_rate + density + pop + pm1029 + avginc, data = dt)
  lm.state <- lm(log(dt[[eval(var)]]) ~ shall + incarc_rate + density + pop + pm1029 + avginc + factor(state), data = dt)
  lm.year <-  lm(log(dt[[eval(var)]]) ~ shall + incarc_rate + density + pop + pm1029 + avginc + factor(year), data = dt)
  lm.all <- lm(log(dt[[eval(var)]]) ~ shall + incarc_rate + density + pop + pm1029 + avginc + factor(year) + factor(state), data = dt)
  plm.all <-  plm(log(dt[[eval(var)]]) ~ shall + incarc_rate + density + pop + pm1029 + avginc + factor(year) + 
                    factor(state), data = dt, model = "pooling")
  coeftest(lm.state, vcov = vcovHC(lm.state, "HC1"))
  waldtest(lm.basic, lm.state, vcov = vcovHC(lm.state, "HC1"))
  
  coeftest(lm.all, vcov = vcovHC(lm.all, "HC1"))
  waldtest(lm.all, lm.state, vcov = vcovHC(lm.all, "HC1"))
  waldtest(lm.all, lm.year, vcov = vcovHC(lm.all, "HC1"))
  
  coeftest(plm.all, vcov = vcovHC(plm.all, "white1", cluster = c("group")))
  
}