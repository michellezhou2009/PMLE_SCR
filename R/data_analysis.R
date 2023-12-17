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
data = BMT %>% 
  mutate(g = factor(g, levels = c(2, 3, 1), 
                    labels = c("AML-low", "AML-high", "ALL")))
PMLE_SCR = function(data, time, death, status_time, status_death,
                    T.fmla = ~ 1, D.fmla = ~ 1, 
                    Gfun, copula.family, 
                    copula.control = 
                      list(link = NULL, formula = ~ 1, initial = NULL)){
  
  link.fun = copula.control$link; copula.fmla = copula.control$formula
  copula.initial = copula.control$initial
  Zmat.nm = unique(c(all.vars(T.fmla), all.vars(D.fmla), 
                     all.vars(copula.fmla)))
  dat = data[, c(time, death, status_time, status_death, Zmat.nm)];
  colnames(dat) = c("time", "death", "status_time", "status_death", Zmat.nm)
  N = nrow(dat); col.nm = colnames(dat)
  Xi = dat[, "time"]; Ci = dat[, "death"]; 
  deltaT = dat[, "status_time"]; deltaD = dat[, "status_death"]
  Zmat.T = model.matrix(T.fmla, dat)[ , -1, drop = F]; n.bT = ncol(Zmat.T)
  Zmat.D = model.matrix(D.fmla, dat)[ , -1, drop = F]; n.bD = ncol(Zmat.D)
  Wmat = model.matrix(copula.fmla, dat); n.gamma = ncol(Wmat)
  tk = sort(unique(Xi[deltaT == 1])); n.tk = length(tk)
  dk = sort(unique(Ci[deltaD == 1])); n.dk = length(dk)
  
  switch(copula.family,
         "Clayton" = {
           copula.index = 3; copula.lwr = 0; copula.upr = 28
           tau.alpha = function(alpha) alpha / (alpha + 2)
           Dtau.alpha = function(alpha){}
           body(Dtau.alpha) = D(expression(alpha/(alpha + 2)), "alpha")
           link.default = "log"
         },
         "Frank" = {
           copula.index = 5; 
           copula.lwr = 0; copula.upr = 50
           tau.alpha = function(alpha){
             fun0 = function(t) t / (exp(t) - 1)
             sapply(alpha, function(a){
               1 - 4 / a + (4 / a ^ 2) * 
                 integrate(fun0, lower = 0, upper = a)$value
             })
           }
           Dtau.alpha = function(alpha){
             fun0 = function(t) t / (exp(t) - 1)
             sapply(alpha, function(a){
               4 / a ^ 2 + 4 / (a * (exp(a) - 1)) - 
                 (8 / a ^ 3) * integrate(fun0, lower = 0, upper = a)$value
             })
           }
           link.default = "identity"
         },
         "Gumbel" = {
           copula.index = 4; 
           copula.lwr = 1; copula.upr = 17
           tau.alpha = function(alpha) 1 - 1 / alpha
           Dtau.alpha = function(alpha) 1 / (alpha ^ 2)
           link.default = "log-1"
         }
  )
  if (is.null(link.fun)) link.fun = link.default
  switch(link.fun,
         "identity" = {
           copula.link = list(h.fun = function(x) {x}, 
                              dot.h.fun = function(x) {rep(1, length(x))},
                              ddot.h.fun = function(x) {rep(0, length(x))})
         },
         "log" = {
           copula.link = list(h.fun = function(x) {exp(x)}, 
                              dot.h.fun = function(x) {exp(x)},
                              ddot.h.fun = function(x) {exp(x)})
         },
         "log-1" = {
           copula.link = list(h.fun = function(x) {exp(x) + 1}, 
                              dot.h.fun = function(x) {exp(x)},
                              ddot.h.fun = function(x) {exp(x)})
         }
  )
  control = list(copula.lwr = copula.lwr, copula.upr = copula.upr)
  
  # Stage I of PMLE: estimate the marginal for death ----
  fitD = fitSPT(dat, time = "death", status = "status_death", 
                formula = D.fmla, Gfun = "PH")
  betaD = as.vector(fitD$beta$est)
  dLambdaD = as.vector(fitD$dLambda$est)
  Psi.thetaD = do.call(cbind, fitD$Psi.theta)
  Psi.uD = predict.fitSPT(fitD, dat)$Psi.surv
  thetaD.est = c(betaD, dLambdaD)
  thetaD.cov = fitD$varcov$robust
  thetaD.se = sqrt(diag(thetaD.cov))
  
  # Initial value for the marginal of relapse for stage II of PMLE ----
  fitT = fitSPT(dat, time = "time", status = "status_time", 
                formula = T.fmla, Gfun = "PH")
  betaT = as.vector(fitT$beta$est)
  dLambdaT  = as.vector(fitT$dLambda$est)
  thetaT.naive = c(betaT, dLambdaT)
  
  # Naive estimate
  objfun<- function(x){
    f = lln.fun(theta = x, thetaT = thetaT.naive, thetaD = thetaD.est, 
                Xi, Ci, deltaT, deltaD, 
                Zmat.T, Zmat.D, copula.index, Gfun, 
                copula.link, Wmat, control)
    g = dlln.fun(theta = x, thetaT = thetaT.naive, thetaD = thetaD.est, 
                 Xi, Ci, deltaT, deltaD, 
                 Zmat.T, Zmat.D, copula.index, Gfun, 
                 copula.link, Wmat, control)
    B = ddlln.fun(theta = x, thetaT = thetaT.naive, thetaD = thetaD.est, 
                  Xi, Ci, deltaT, deltaD, 
                  Zmat.T, Zmat.D, copula.index, Gfun, 
                  copula.link, Wmat, control) 
    list(value = f, gradient = g, hessian = B)
  }
  est.res <- trust(
    objfun, copula.initial, 5, 100, iterlim = 300, 
    minimize= FALSE, blather = T)
  gamma.naive = est.res$argument
  para.naive = data.frame(type = "copula", para = colnames(Wmat), time = NA,
                          est = gamma.naive)
  
  para.naive = rbind(
    para.naive,
    data.frame(type = "betaT", para = colnames(Zmat.T), time = NA,
               est = thetaT.naive[c(1 : n.bT)])
  )
  para.naive = rbind(
    para.naive,
    data.frame(type = "dLambdaT", para = "dLambdaT", time = tk,
               est = thetaT.naive[n.bT + c(1 : n.tk)])
  )
  
  ## ===== PMLE ==== ##
  PMLE.ini = c(gamma.naive, thetaT.naive)
  objfun<- function(x){
    f = lln.fun(theta = x, thetaT = NULL, thetaD = thetaD.est, 
                Xi, Ci, deltaT, deltaD, 
                Zmat.T, Zmat.D, copula.index, Gfun, 
                copula.link, Wmat, control)
    g = dlln.fun(theta = x, thetaT = NULL, thetaD = thetaD.est, 
                 Xi, Ci, deltaT, deltaD, 
                 Zmat.T, Zmat.D, copula.index, Gfun, 
                 copula.link, Wmat, control)
    B = ddlln.fun(theta = x, thetaT = NULL, thetaD = thetaD.est, 
                  Xi, Ci, deltaT, deltaD, 
                  Zmat.T, Zmat.D, copula.index, Gfun, 
                  copula.link, Wmat, control) 
    list(value = f, gradient = g, hessian = B)
  }
  
  est.res <- trust(
    objfun, PMLE.ini, 5, 100, iterlim = 300, 
    minimize= FALSE, blather = T)
  theta.est = est.res$argument
  Imat = - ddlln.fun(theta = theta.est, thetaT = NULL, thetaD = thetaD.est, 
                     Xi, Ci, deltaT, deltaD, 
                     Zmat.T, Zmat.D, copula.index, Gfun, 
                     copula.link, Wmat, control) 
  dll.k = ddll.fun2(theta = theta.est, thetaD = thetaD.est, 
                    Xi, Ci, deltaT, deltaD, 
                    Zmat.T, Zmat.D, copula.index, Gfun, 
                    copula.link, Wmat, control) 
  dll.k[is.na(dll.k)] = 0
  Psi =  dll.fun(theta = theta.est, thetaT = NULL, thetaD = thetaD.est, 
                 Xi, Ci, deltaT, deltaD, 
                 Zmat.T, Zmat.D, copula.index, Gfun, 
                 copula.link, Wmat, control) 
  Psi[is.na(Psi)] = 0
  Psi.theta = Psi + Psi.uD %*% dll.k / N
  Vmat = crossprod(Psi.theta, Psi.theta) / N
  theta.cov = solve(Imat) %*% Vmat %*% solve(Imat) / N
  
  
  para.est = data.frame(type = "copula", para = colnames(Wmat), time = NA,
                        est = theta.est[1 : n.gamma], 
                        se = sqrt(diag(theta.cov))[1 : n.gamma])
  
  para.est = rbind(
    para.est,
    data.frame(type = "betaT", para = colnames(Zmat.T), time = NA,
               est = theta.est[n.gamma + c(1 : n.bT)], 
               se = sqrt(diag(theta.cov))[n.gamma + c(1 : n.bT)])
  )
  para.est = rbind(
    para.est,
    data.frame(type = "dLambdaT", para = "dLambdaT", time = tk,
               est = theta.est[n.gamma + n.bT + c(1 : n.tk)],
               se = sqrt(diag(theta.cov))[n.gamma + n.bT + c(1 : n.tk)])
  )
  para.est = rbind(
    para.est,
    data.frame(type = "betaD", para = colnames(Zmat.D), time = NA,
               est = thetaD.est[c(1 : n.bD)], 
               se = sqrt(diag(thetaD.cov))[c(1 : n.bD)])
  )
  para.est = rbind(
    para.est,
    data.frame(type = "dLambdaD", para = "dLambdaD", time = dk,
               est = thetaD.est[n.bD + c(1 : n.dk)],
               se = sqrt(diag(thetaD.cov))[n.bD + c(1 : n.dk)])
  )
  
  theta.cov = list(
    copula = theta.cov[1 : n.gamma, 1 : n.gamma, drop = F],
    thetaT = theta.cov[n.gamma + c(1 : (n.bT + n.tk)), 
                       n.gamma + c(1 : (n.bT + n.tk))],
    thetaD = thetaD.cov)
  
  out.PMLE = list(est = para.est, cov = theta.cov)
  
  ## ===== MLE ==== ##
  MLE.ini = c(gamma.naive, thetaT.naive, thetaD.est)
  objfun<- function(x){
    f = lln.fun(theta = x, thetaT = NULL, thetaD = NULL, 
                Xi, Ci, deltaT, deltaD, 
                Zmat.T, Zmat.D, copula.index, Gfun, 
                copula.link, Wmat, control)
    g = dlln.fun(theta = x, thetaT = NULL, thetaD = NULL, 
                 Xi, Ci, deltaT, deltaD, 
                 Zmat.T, Zmat.D, copula.index, Gfun, 
                 copula.link, Wmat, control)
    B = ddlln.fun(theta = x, thetaT = NULL, thetaD = NULL, 
                  Xi, Ci, deltaT, deltaD, 
                  Zmat.T, Zmat.D, copula.index, Gfun, 
                  copula.link, Wmat, control) 
    list(value = f, gradient = g, hessian = B)
  }
  
  est.res <- trust(
    objfun, MLE.ini, 5, 100, iterlim = 300, 
    minimize= FALSE, blather = T)
  theta.est = est.res$argument
  Imat = - ddlln.fun(theta = theta.est, thetaT = NULL, thetaD = NULL, 
                     Xi, Ci, deltaT, deltaD, 
                     Zmat.T, Zmat.D, copula.index, Gfun, 
                     copula.link, Wmat, control) 
  Psi =  dll.fun(theta = theta.est, thetaT = NULL, thetaD = NULL, 
                 Xi, Ci, deltaT, deltaD, 
                 Zmat.T, Zmat.D, copula.index, Gfun, 
                 copula.link, Wmat, control) 
  Vmat = crossprod(Psi, Psi) / N
  theta.cov = solve(Imat) %*% Vmat %*% solve(Imat) / N
  theta.cov.model = solve(Imat) / N
  
  para.est = data.frame(
    type = "copula", para = colnames(Wmat), time = NA,
    est = theta.est[1 : n.gamma], 
    se = sqrt(diag(theta.cov))[1 : n.gamma],
    se.model = sqrt(diag(theta.cov.model))[1 : n.gamma])
  
  para.est = rbind(
    para.est,
    data.frame(
      type = "betaT", para = colnames(Zmat.T), time = NA,
      est = theta.est[n.gamma + c(1 : n.bT)], 
      se = sqrt(diag(theta.cov))[n.gamma + c(1 : n.bT)],
      se.model = sqrt(diag(theta.cov.model))[n.gamma + c(1 : n.bT)])
  )
  para.est = rbind(
    para.est,
    data.frame(
      type = "dLambdaT", para = "dLambdaT", time = tk,
      est = theta.est[n.gamma + n.bT + c(1 : n.tk)],
      se = sqrt(diag(theta.cov))[n.gamma + n.bT + c(1 : n.tk)],
      se.model = sqrt(diag(theta.cov.model))[n.gamma + n.bT + c(1 : n.tk)])
  )
  para.est = rbind(
    para.est,
    data.frame(
      type = "betaD", para = colnames(Zmat.D), time = NA,
      est = theta.est[n.gamma + n.bT + n.tk + c(1 : n.bD)], 
      se = sqrt(diag(theta.cov))[n.gamma + n.bT + n.tk + c(1 : n.bD)],
      se.model = 
        sqrt(diag(theta.cov.model))[n.gamma + n.bT + n.tk + c(1 : n.bD)])
  )
  para.est = rbind(
    para.est,
    data.frame(
      type = "dLambdaD", para = "dLambdaD", time = dk,
      est = theta.est[n.gamma + n.bT + n.tk + n.bD + c(1 : n.dk)],
      se = sqrt(diag(theta.cov))[n.gamma + n.bT + n.tk + n.bD + c(1 : n.dk)],
      se.model = 
        sqrt(diag(theta.cov.model))[n.gamma + n.bT + n.tk + n.bD + c(1 : n.dk)])
  )
  theta.cov = list(
    copula = theta.cov[1 : n.gamma, 1 : n.gamma, drop = F],
    thetaT = theta.cov[n.gamma + c(1 : (n.bT + n.tk)), 
                       n.gamma + c(1 : (n.bT + n.tk))],
    thetaD = theta.cov[n.gamma + n.bT + n.tk + c(1 : (n.bD + n.dk)), 
                       n.gamma + n.bT + n.tk + c(1 : (n.bD + n.dk))])
  
  out.MLE = list(est = para.est, cov = theta.cov)
  
  list(PMLE = out.PMLE, MLE = out.MLE, naive = para.naive,
       copula.family = copula.family, copula.link = copula.link,
       Tau2Par = list(tau.alpha = tau.alpha, Dtau.alpha = Dtau.alpha))
}

# Anaysis with different copula families ----
copula.control = list(link = "identity", formula = ~ g, initial = c(2, 0, 0))
copula.fam.all = c("Gumbel", "Clayton", "Frank")
res.all = lapply(copula.fam.all, function(copula.fam){
  PMLE_SCR(data = data, time = "T2", death = "T1", 
           status_time = "delta2", status_death = "delta1", 
           T.fmla = ~ g, D.fmla = ~ g, Gfun = list(T = "PH", D = "PH"),
           copula.family = copula.fam, 
           copula.control = copula.control)
})

# Create Table 6 ----
tab.all = lapply(res.all, function(res){
  out = res$PMLE
  res.betaT = out$est %>% filter(type == "betaT")
  out.betaT = data.frame(
    type = "betaT",
    para = c("AML high", "ALL"),
    est = paste0(
      format(round(res.betaT$est, 3), nsmall = 3, digits = 4),
      " (", format(round(res.betaT$se, 3), nsmall = 3, digits = 4), ")")
  )
  res.betaD = out$est %>% filter(type == "betaD")
  out.betaD = data.frame(
    type = "betaD",
    para = c("AML high", "ALL"),
    est = paste0(
      format(round(res.betaD$est, 3), nsmall = 3, digits = 4),
      " (", format(round(res.betaD$se, 3), nsmall = 3, digits = 4), ")")
  )
  res.gamma = out$est %>% filter(type == "copula")
  gamma.cov =  out$cov$copula
  newdata = data.frame(g = levels(data$g)) %>%
    mutate(g = factor(g, levels = levels(data$g)))
  Wmat = model.matrix(copula.control$formula, newdata)
  gamma.est = matrix(res.gamma$est[match(colnames(Wmat), res.gamma$para)], 
                     byrow = F, ncol = 1)
  lp.est = Wmat %*% gamma.est
  lp.se = sqrt(diag(Wmat %*% gamma.cov %*% t(Wmat)))
  alpha.est = res$copula.link$h.fun(lp.est)
  alpha.se = abs(res$copula.link$dot.h.fun(lp.est)) * lp.se
  tau.est = res$Tau2Par$tau.alpha(alpha.est)
  tau.se = abs(res$Tau2Par$Dtau.alpha(alpha.est)) * alpha.se
  out.tau = data.frame(
    type = "tau", para = c("AML low", "AML high", "ALL"),
    est = paste0(format(round(tau.est, 3), nsmall = 3, digits = 4),
                 " (", format(round(tau.se, 3), nsmall = 3, digits = 4), ")")
  )
  out.PMLE = rbind(out.betaT, out.betaD, out.tau) 
  colnames(out.PMLE) = c("type", "para", paste0("PMLE_", res$copula.family))  
  
  out = res$MLE
  res.betaT = out$est %>% filter(type == "betaT")
  out.betaT = data.frame(
    type = "betaT",
    para = c("AML high", "ALL"),
    est = paste0(
      format(round(res.betaT$est, 3), nsmall = 3, digits = 4),
      " (", format(round(res.betaT$se, 3), nsmall = 3, digits = 4), ")")
  )
  res.betaD = out$est %>% filter(type == "betaD")
  out.betaD = data.frame(
    type = "betaD",
    para = c("AML high", "ALL"),
    est = paste0(
      format(round(res.betaD$est, 3), nsmall = 3, digits = 4),
      " (", format(round(res.betaD$se, 3), nsmall = 3, digits = 4), ")")
  )
  res.gamma = out$est %>% filter(type == "copula")
  gamma.cov =  out$cov$copula
  newdata = data.frame(g = levels(data$g)) %>%
    mutate(g = factor(g, levels = levels(data$g)))
  Wmat = model.matrix(copula.control$formula, newdata)
  gamma.est = matrix(res.gamma$est[match(colnames(Wmat), res.gamma$para)], 
                     byrow = F, ncol = 1)
  lp.est = Wmat %*% gamma.est
  lp.se = sqrt(diag(Wmat %*% gamma.cov %*% t(Wmat)))
  alpha.est = res$copula.link$h.fun(lp.est)
  alpha.se = abs(res$copula.link$dot.h.fun(lp.est)) * lp.se
  tau.est = res$Tau2Par$tau.alpha(alpha.est)
  tau.se = abs(res$Tau2Par$Dtau.alpha(alpha.est)) * alpha.se
  out.tau = data.frame(
    type = "tau", para = c("AML low", "AML high", "ALL"),
    est = paste0(format(round(tau.est, 3), nsmall = 3, digits = 4),
                 " (", format(round(tau.se, 3), nsmall = 3, digits = 4), ")")
  )
  out.MLE = rbind(out.betaT, out.betaD, out.tau) 
  colnames(out.MLE) = c("type", "para", paste0("MLE_", res$copula.family))  
  
  left_join(out.PMLE, out.MLE)
})
tab = tab.all[[1]] %>%
  left_join(., tab.all[[2]]) %>% left_join(., tab.all[[3]])
write.csv(tab, "tab6_BMT.csv")

# Create Figure 1
Xi = data$T2; Ci = data$T1;
deltaT = data$delta2; deltaD = data$delta1
tk = sort(unique(Xi[deltaT == 1]))
dk = sort(unique(Ci[deltaD == 1]))

newdatT = purrr::cross(list(time = tk, g = c("AML-low", "AML-high", "ALL")))
newdatD = purrr::cross(list(time = dk, g = c("AML-low", "AML-high", "ALL")))

newdatT = tidyr::crossing(time = tk, g = c("AML-low", "AML-high", "ALL"))
newdatD = tidyr::crossing(time = dk, g = c("AML-low", "AML-high", "ALL"))

S.fun = function(b, dLambda, tt0, zz0){
  b = matrix(b, byrow = F, ncol = 1)
  tt = dLambda$time
  dLambda = dLambda$est
  Lambda.tt0 = sapply(tt0, function(t0) sum(dLambda[tt <= t0]))
  exp(- Lambda.tt0 *  exp(zz0 %*% b))
}

dot.S.fun = function(b, dLambda, tt0, zz0){
  b = matrix(b, byrow = F, ncol = 1)
  tt = dLambda$time
  dLambda = dLambda$est
  Lambda.tt0 = sapply(tt0, function(t0) sum(dLambda[tt <= t0]))
  dot.b = - as.vector(exp(- Lambda.tt0 *  exp(zz0 %*% b)) * 
                        Lambda.tt0 * exp(zz0 %*% b)) * zz0
  dot.dLambda = - as.vector(
    exp(- Lambda.tt0 * exp(zz0 %*% b)) * exp(zz0 %*% b))* 
    (matrix(tt, byrow = T, nrow = length(tt0), ncol = length(tt)) <= tt0)
  cbind(dot.b, dot.dLambda)
}

out.betaT = lapply(res.all, function(res){
  out = res$PMLE
  res.theta = out$est %>% filter(type %in% c("betaT", "dLambdaT"))
  theta.cov = out$cov$thetaT
  newdat = newdatT
  b = res.theta$est[res.theta$type == "betaT"]
  dLambda = res.theta[res.theta$type == "dLambdaT", c("time", "est")]
  newdat = newdat %>% 
    mutate(g = factor(g, levels = c("AML-low", "AML-high", "ALL")))
  zz0 = model.matrix(~ g, newdat)[, -1, drop = F]
  newdat$est = as.vector(S.fun(b, dLambda, newdat$time, zz0))
  dd = dot.S.fun(b, dLambda, newdat$time, zz0)
  cov.mat = dd %*% theta.cov %*% t(dd)
  newdat$se = sqrt(diag(cov.mat))
  newdat.PMLE = newdat %>% mutate(method = "PMLE")
  
  out = res$MLE
  res.theta = out$est %>% filter(type %in% c("betaT", "dLambdaT"))
  theta.cov = out$cov$thetaT
  newdat = newdatT
  b = res.theta$est[res.theta$type == "betaT"]
  dLambda = res.theta[res.theta$type == "dLambdaT", c("time", "est")]
  newdat = newdat %>% 
    mutate(g = factor(g, levels = c("AML-low", "AML-high", "ALL")))
  zz0 = model.matrix(~ g, newdat)[, -1, drop = F]
  newdat$est = as.vector(S.fun(b, dLambda, newdat$time, zz0))
  dd = dot.S.fun(b, dLambda, newdat$time, zz0)
  cov.mat = dd %*% theta.cov %*% t(dd)
  newdat$se = sqrt(diag(cov.mat))
  newdat.MLE = newdat %>% mutate(method = "MLE") 
  
  rbind(newdat.PMLE, newdat.MLE) %>% mutate(family = res$copula.family)
}) %>% do.call(rbind, .)

library(ggplot2)
p1 = out.betaT %>% filter(method == "PMLE") %>%
  mutate(lwr = est - 1.96 * se, upr = est + 1.96 * se,
         family = factor(family, levels = copula.fam.all)) %>%
  ggplot(.) + 
  geom_step(aes(x = time, y = est, linewidth = g), linetype = "solid") +
  geom_step(aes(x = time, y = lwr, linewidth = g), linetype = "dashed") +
  geom_step(aes(x = time, y = upr, linewidth = g), linetype = "dashed") +
  
  scale_linewidth_manual(
    values = c("AML-low" = 0.2, "AML-high" = 0.5, "ALL" = 1.2), name = "") +
  facet_grid(~ family) +
  labs(x = "time", y = "replapse-free survival rate") 

pdf("fig1_BMT.pdf", height = 5, width = 10)  
p1 + theme(text = element_text(size = 14))
dev.off()
  

