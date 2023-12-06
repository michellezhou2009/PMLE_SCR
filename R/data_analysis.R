rm(list=ls())
source("R/fitSPT.R")
source("R/helpers.R")
source("R/SemiCompCopSTPFuns.R")
library(dplyr)
library(VineCopula)
library(trust)
library(SemiCompRisks)

# Read and process data ----
data("BMT")
workdata = BMT %>% 
  select_("disease" = "g", "t_relapse" = "T2", "delta_relapse" = "delta2",
          "t_death"= "T1", "delta_death" = "delta1") %>% 
  mutate(disease = factor(disease,levels = c(2, 3, 1), 
                          labels = c("AML-low", "AML-high", "ALL"))) 

Xi = workdata$t_relapse
Ci = workdata$t_death
deltaT = workdata$delta_relapse
deltaD = workdata$delta_death
Z1 = workdata$disease
data = data.frame(Xi, Ci, deltaT, deltaD, Z1)
var1 = var2 = "Z1"; 
time1 = "Xi"; time2 = "Ci"; status1 = "deltaT"; status2 = "deltaD"
N = nrow(data); col.nm = colnames(data)
formula.marginal = list(
  T = paste0(" ~ ", paste0(var1, collapse = "+")), 
  D = paste0(" ~ ", paste0(var2, collapse = "+"))
)
T.fmla = as.formula(formula.marginal$T)
D.fmla = as.formula(formula.marginal$D)
Xi = data[, col.nm == time1]; Ci = data[, col.nm == time2]; 
deltaT = data[, col.nm == status1]; deltaD = data[, col.nm == status2]
Zmat.T = model.matrix(T.fmla, data)[ , -1, drop = F]
Zmat.D = model.matrix(D.fmla, data)[ , -1, drop = F]
tk = sort(Xi[deltaT == 1]); n.tk = length(unique(tk))
dk = sort(Ci[deltaD == 1]); n.dk = length(unique(dk))

# Stage I of PMLE: estimate the marginal for death ----
fitD = fitSPT(data, time = time2, status = status2, 
              formula = ~ Z1, Gfun = "PH")
betaD = as.vector(fitD$beta$est)
dLambdaD = as.vector(fitD$dLambda$est)
Psi.thetaD = do.call(cbind, fitD$Psi.theta)
Psi.uD = predict.fitSPT(fitD, data)$Psi.surv
thetaD.se = sqrt(diag(fitD$varcov$robust))

# Initial value for the marginal of relapse for stage II of PMLE ----
fitT = fitSPT(data, time = time1, status = status1, 
              formula = ~ Z1, Gfun = "PH")
betaT = as.vector(fitT$beta$est)
dLambdaT  = as.vector(fitT$dLambda$est)


# Create "exc.fun" ----
## obtain MLE and PMLE for a given copula family
exc.fun = function(copula.fam){
  switch(copula.fam,
         "Clayton" = {
           copula.index = 3; 
           copula.lower = 0.1; copula.upper = 28
           tau.alpha = function(alpha) alpha / (alpha + 2)
           Dtau.alpha = function(alpha){}
           body(Dtau.alpha) = D(expression(alpha/(alpha + 2)), "alpha")
         },
         "Frank" = {
           copula.index = 5; 
           copula.lower = 0.1; copula.upper = 50
           tau.alpha = function(alpha){
             fun0 = function(t) t / (exp(t) - 1)
             1 - 4 / alpha + (4 / alpha ^ 2) * 
               integrate(fun0, lower = 0, upper = alpha)$value
           }
           Dtau.alpha = function(alpha){
             fun0 = function(t) t / (exp(t) - 1)
             4 / alpha ^ 2 + 4 / (alpha * (exp(alpha) - 1)) - 
               (8 / alpha ^ 3) * integrate(fun0, lower = 0, upper = alpha)$value
           }
         },
         "Gumbel" = {
           copula.index = 4; 
           copula.lower = 1; copula.upper = 17
           tau.alpha = function(alpha) 1 - 1 / alpha
           Dtau.alpha = function(alpha) 1 / (alpha ^ 2)
           }
  )
  ## initial value of copula parameter
  Gfun = list(T = "PH", D = "PH")
  obj.fun <- function(alpha, betaT, dLambdaT, betaD, dLambdaD, 
                      Xi, Ci, deltaT, deltaD, 
                      Zmat.T, Zmat.D, copula.index, Gfun){
    theta = c( betaT, dLambdaT, betaD, dLambdaD)
    lln.fun1(copula.para=alpha,theta=theta, thetaD = NULL, Xi = Xi, Ci = Ci, 
             deltaT = deltaT, deltaD = deltaD, Zmat.T = Zmat.T, Zmat.D = Zmat.D, 
             copula.index = copula.index, Gfun = Gfun)
  }
  alpha = optimize(obj.fun, lower = copula.lower, upper = copula.upper, 
                   betaT = betaT, dLambdaT = dLambdaT, 
                   betaD = betaD, dLambdaD = dLambdaD, 
                   Xi = Xi, Ci = Ci, deltaT = deltaT, deltaD = deltaD,
                   Zmat.T = Zmat.T, Zmat.D = Zmat.D, 
                   copula.index = copula.index, Gfun = Gfun,
                   maximum = T)$maximum
  
  theta.ini = c(alpha, betaT, dLambdaT, betaT, dLambdaD)
  theta.ini.0 = c(alpha, betaT, dLambdaT)
  thetaD.est = c(betaD, dLambdaD)
  
  ## ===== PMLE ==== ##
  objfun<- function(x){
    f = lln.fun(theta = x, thetaD = thetaD.est, Xi, Ci, deltaT, deltaD, 
                Zmat.T, Zmat.D, copula.index, Gfun)
    g = dlln.fun(theta = x, thetaD = thetaD.est, Xi, Ci, deltaT, deltaD, 
                 Zmat.T, Zmat.D, copula.index, Gfun)
    B = ddlln.fun(theta = x, thetaD = thetaD.est, Xi, Ci, deltaT, deltaD, 
                  Zmat.T, Zmat.D, copula.index, Gfun)
    list(value = f, gradient = g, hessian = B)
  }
  est.res <- trust(
    objfun, theta.ini.0, 5, 100, iterlim = 300, 
    minimize= FALSE, blather = T)
  theta.est = est.res$argument
  names(theta.est) = c(
    "copula", paste0("T.",colnames(Zmat.T)), paste0("dLambdaT.", 1 : n.tk))
  Imat = - ddlln.fun(theta = theta.est, thetaD = thetaD.est, 
                     Xi, Ci, deltaT, deltaD, 
                     Zmat.T, Zmat.D, copula.index, Gfun)
  dll.k = ddll.fun2(theta = theta.est, thetaD = thetaD.est, 
                    Xi, Ci, deltaT, deltaD, 
                    Zmat.T, Zmat.D, copula.index, Gfun) 
  dll.k[is.na(dll.k)] = 0
  Psi =  dll.fun(theta = theta.est, thetaD = thetaD.est, 
                 Xi, Ci, deltaT, deltaD, 
                 Zmat.T, Zmat.D, copula.index, Gfun, keep = T) 
  Psi[is.na(Psi)] = 0
  Psi = Psi + Psi.uD %*% dll.k / N
  Vmat = crossprod(Psi, Psi) / N
  theta1.cov = solve(Imat) %*% Vmat %*% solve(Imat) / N
  colnames(theta1.cov) = rownames(theta1.cov) = names(theta.est)
  theta1.se = sqrt(diag(theta1.cov))
  theta1.est = theta.est
  Psi.theta1 = Psi %*% solve(Imat)
  theta.est = c(theta1.est, thetaD.est)
  names(theta.est) = c(names(theta1.est), 
                       paste0("D.",colnames(Zmat.D)), paste0("dLambdaD.", 1 : n.dk))
  tt = c(NA, rep(NA, length(betaT)), tk, rep(NA, length(betaD)), dk)
  theta.se = c(theta1.se, thetaD.se)
  Psi = cbind(Psi.theta1, Psi.thetaD)
  index = grep("copula", names(theta.est))
  alpha.est = theta.est[index]
  alpha.se = theta.se[index]
  index = grep("T.Z", names(theta.est))
  betaT.est = theta.est[index]
  betaT.se = theta.se[index]
  index = grep("D.Z", names(theta.est))
  betaD.est = theta.est[index]
  betaD.se = theta.se[index]
  tau.est = tau.alpha(alpha.est)
  tau.se = abs(Dtau.alpha(alpha.est)) * alpha.se
  out.PMLE = data.frame(
      family = copula.fam, method = "PMLE",
      para = c("copula", "tau", 
               "betaT_High", "betaT_All", "betaD_High", "betaD_All", 
               "iterations"),
      est = c(alpha.est, tau.est, betaT.est, betaD.est, est.res$iterations),
      se = c(alpha.se, tau.se, betaT.se, betaD.se, NA)
    )
  
  ## ============= MLE =========== ##
  objfun<- function(x){
    f = lln.fun(theta = x,thetaD = NULL, Xi, Ci, deltaT, deltaD, 
                Zmat.T, Zmat.D, copula.index, Gfun)
    g = dlln.fun(theta = x,thetaD = NULL, Xi, Ci, deltaT, deltaD, 
                 Zmat.T, Zmat.D, copula.index, Gfun)
    B = ddlln.fun(theta = x,thetaD = NULL, Xi, Ci, deltaT, deltaD, 
                  Zmat.T, Zmat.D, copula.index, Gfun)
    list(value = f, gradient = g, hessian = B)
  }
  est.res <- trust(
    objfun, theta.ini, 5, 100, iterlim = 300, 
    minimize= FALSE, blather = T)
  theta.est = est.res$argument
  names(theta.est) = c(
    "copula", paste0("T.",colnames(Zmat.T)), paste0("dLambdaT.", 1 : n.tk),
    paste0("D.",colnames(Zmat.D)), paste0("dLambdaD.", 1 : n.dk))
  
  tt = c(NA, rep(NA, length(betaT)), tk, rep(NA, length(betaD)), dk)
  Imat = - ddlln.fun(theta = theta.est, thetaD = NULL, Xi, Ci, deltaT, deltaD, 
                     Zmat.T, Zmat.D, copula.index)
  Psi = dll.fun(theta = theta.est, thetaD = NULL, Xi, Ci, deltaT, deltaD, 
                Zmat.T, Zmat.D, copula.index)
  Vmat = crossprod(Psi, Psi) / N
  est.cov = solve(Imat) %*% Vmat %*% solve(Imat) / N
  colnames(est.cov) = rownames(est.cov) = names(theta.est)
  est.se = sqrt(diag(est.cov))
  Psi = Psi %*% solve(Imat)
  index = grep("copula", names(theta.est))
  alpha.est = theta.est[index]
  alpha.se = est.se[index]
  index = grep("T.Z", names(theta.est))
  betaT.est = theta.est[index]
  betaT.se = est.se[index]
  index = grep("D.Z", names(theta.est))
  betaD.est = theta.est[index]
  betaD.se = est.se[index]
  tau.est = tau.alpha(alpha.est)
  tau.se = abs(Dtau.alpha(alpha.est)) * alpha.se
  out.MLE = data.frame(
        family = copula.fam, method = "MLE",
        para = c("copula", "tau", 
                 "betaT_High", "betaT_All", "betaD_High", "betaD_All", 
                 "iterations"),
        est = c(alpha.est, tau.est, betaT.est, betaD.est, est.res$iterations),
        se = c(alpha.se, tau.se, betaT.se, betaD.se, NA)
      )
  
  rbind(out.MLE, out.PMLE)
}

# Fit Gumbel, Clayton, and Frank copula ----

fam.all = c("Gumbel", "Clayton", "Frank")
res.all = lapply(fam.all, function(fam){
  exc.fun(fam)
}) %>% do.call(rbind, .)

# Create Table 6 ----

tab = res.all %>% filter(para != "iterations") %>%
  mutate(res = paste0(format(round(est, 3), nsmall = 3, digits = 4), " (", 
                      format(round(se, 3), nsmall = 3, digits = 4), ")")) %>%
 dplyr::select(family, para, method, res) %>%
  tidyr::spread(key = "method", value = "res") %>%
  mutate(
    family = factor(family, levels = c("Gumbel", "Clayton", "Frank")),
    para = factor(para, levels = c("copula", "tau", "betaT_High", "betaT_All",
                                   "betaD_High", "betaD_All"))
  ) %>%
  arrange(family, para)
write.csv(tab, "Tables/tab6_BMT.csv", row.names = F)
