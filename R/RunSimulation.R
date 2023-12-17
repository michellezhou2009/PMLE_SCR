rm(list = ls())
library(parallel)
library(doSNOW)
library(foreach)
library(VineCopula)

execute.fun = function(nrep, nsamp, tau, copula.true, copula.fam, theta.true, 
                       gfun=list(Gfun1 = "PH", Gfun2 = "PH"), seed1=1234, 
                       res.dir){
  eta.true = c(0.2, 0)
  switch(copula.true,
         "Clayton" = {index.true = 3},
         "Frank" = {index.true = 5},
         "Gumbel" = {index.true = 4}
  )
  alpha.true = BiCopTau2Par(family = index.true, tau = tau)
  
  fam.nm = paste0(substr(copula.true, 1, 1), substr(copula.fam, 1, 1))
  
  set.seed(seed1)
  seed.all = round(runif(nrep, 10, 100000000))
  
  foreach(i = 1:nrep, 
          .packages=c("survival", "VineCopula", "trust", "dplyr", "tidyr")
  ) %dopar% {
    file.nm =  paste0(res.dir, fam.nm, "_nsamp", nsamp, "_", "tau", tau, 
                      "_", "theta", theta.true[1], ",", theta.true[2], "_rep", 
                      i, ".RData")
    if (!file.exists(file.nm)){
      cat("rep", i)
      source("R/helpers.R")
      source("R/fitSPT.R")
      source("R/SemiCompCopSTPFuns.R")
      set.seed(seed.all[i])
      
      # data generation 
      zz1 = rnorm(nsamp, 1,0.5); zz1[zz1>=2] = 2; zz1[zz1<=0] = 0 
      zz2 = rbinom(nsamp, 1,0.8)
      uu = BiCopSim(N=nsamp, family = index.true, par = alpha.true)
      ee = log(-log(uu))
      Ai = runif(nsamp, min = 1, max = 10)
      
      T1i = 3* exp(-theta.true[1]*zz1 - theta.true[2]*zz2 + ee[,1])
      T2i = 3* exp(-eta.true[1]*zz1 - eta.true[2]*zz2 + ee[,2])
      Zi = pmin(T1i, T2i)
      Xi = pmin(Zi, Ai)
      Ci = pmin(T2i, Ai)
      deltaT = 1 * (T1i <= Ci)
      deltaD = 1 * (T2i <= Ai)
      delta0 = 1 * (Zi <= Ai)
      
      mydata =data.frame(Xi = Xi, Ci = Ci,  deltaT = deltaT, deltaD = deltaD, 
                         Z1 = zz1, Z2 = zz2)
      
      data = mydata; var1 = var2 = c("Z1", "Z2"); 
      time1 = "Xi"; time2 = "Ci"; status1 = "deltaT"; status2 = "deltaD"
      switch(copula.fam,
             "Clayton" = {copula.index = 3; copula.lower = 0.1; copula.upper = 28},
             "Frank" = {copula.index = 5; copula.lower = 0.1; copula.upper = 50},
             "Gumbel" = {copula.index = 4; copula.lower = 1; copula.upper = 17}
      )
      
      N = nrow(data); col.nm = colnames(data)
      formula.marginal = list(T = paste0(" ~ ", paste0(var1, collapse = "+")), 
                              D = paste0(" ~ ", paste0(var2, collapse = "+")))
      T.fmla = as.formula(formula.marginal$T)
      D.fmla = as.formula(formula.marginal$D)
      Xi = data[, col.nm == time1]; Ci = data[, col.nm == time2]; 
      deltaT = data[, col.nm == status1]; deltaD = data[, col.nm == status2]
      Zmat.T = model.matrix(T.fmla, data)[ , -1, drop = F]
      Zmat.D = model.matrix(D.fmla, data)[ , -1, drop = F]
      tk = sort(Xi[deltaT == 1]); n.tk = length(tk)
      dk = sort(Ci[deltaD == 1]); n.dk = length(dk)
      
      fitD = fitSPT(data, time = time2, status = status2, formula = D.fmla, 
                    Gfun = gfun$Gfun2)
      betaD = as.vector(fitD$beta$est)
      dLambdaD = as.vector(fitD$dLambda$est)
      Psi.thetaD = do.call(cbind, fitD$Psi.theta)
      Psi.uD = predict.fitSPT(fitD, data)$Psi.surv
      thetaD.se = sqrt(diag(fitD$varcov$robust))
      
      ## ==== initial values ======= ##
      fitT = fitSPT(data, time = time1, status = status1, 
                    formula = T.fmla, Gfun = gfun$Gfun1)
      betaT = as.vector(fitT$beta$est)
      dLambdaT  = as.vector(fitT$dLambda$est)
      
      ## initial value of copula parameter
      Gfun = list(T = gfun$Gfun1, D = gfun$Gfun2)
      obj.fun <- function(alpha, betaT, dLambdaT, betaD, dLambdaD, 
                          Xi, Ci, deltaT, deltaD, 
                          Zmat.T, Zmat.D, copula.index, Gfun){
        theta = c(alpha, betaT, dLambdaT, betaD, dLambdaD)
        lln.fun(theta, thetaD = NULL, Xi = Xi, Ci = Ci, deltaT = deltaT, deltaD = deltaD,
                Zmat.T = Zmat.T, Zmat.D = Zmat.D, 
                copula.index = copula.index, Gfun = Gfun)
      }
      alpha = optimize(obj.fun, lower = copula.lower, upper = copula.upper, 
                       betaT = betaT, dLambdaT = dLambdaT, 
                       betaD = betaD, dLambdaD = dLambdaD, 
                       Xi = Xi, Ci = Ci, deltaT = deltaT, deltaD = deltaD,
                       Zmat.T = Zmat.T, Zmat.D = Zmat.D, 
                       copula.index = copula.index, Gfun = Gfun,
                       maximum = T)$maximum
      
      theta.ini.0 = c(alpha, betaT, dLambdaT)
      thetaD.est = c(betaD, dLambdaD)
      theta.ini = c(alpha, betaT, dLambdaT, betaD, dLambdaD)
      
      control = list(copula.lwr = copula.lower, copula.upr = copula.upper)
      ## ====== MLE ====== ##
      objfun<- function(x){
        f = lln.fun(theta = x, thetaT = NULL, thetaD = NULL, 
                    Xi, Ci, deltaT, deltaD, Zmat.T, Zmat.D, copula.index, Gfun, 
                    control = control)
        g = dlln.fun(theta = x, thetaT = NULL, thetaD = NULL, 
                     Xi, Ci, deltaT, deltaD, Zmat.T, Zmat.D, copula.index, Gfun, 
                     control = control)
        B = ddlln.fun(theta = x, thetaT = NULL, thetaD = NULL, 
                      Xi, Ci, deltaT, deltaD, Zmat.T, Zmat.D, copula.index, Gfun, 
                      control = control)
        list(value = f, gradient = g, hessian = B)
      }
      
      usedtime <- system.time({
        est.res <- trust(
          objfun, theta.ini, 5, 100, iterlim = 300, 
          minimize= FALSE, blather = T)})
      theta.est = est.res$argument; 
      names(theta.est) = c(
        "copula", paste0("T.",var1), paste0("dLambdaT.", 1 : n.tk), 
        paste0("D.",var2), paste0("dLambdaD.", 1 : n.dk))
      tt = c(NA, rep(NA, length(betaT)), tk, rep(NA, length(betaD)), dk)
      Imat = - ddlln.fun(theta = theta.est, thetaT = NULL, thetaD = NULL, Xi, Ci, 
                         deltaT, deltaD, Zmat.T, Zmat.D, copula.index, 
                         Gfun = Gfun, control = control)
      Psi = dll.fun(theta = theta.est, thetaT = NULL, thetaD = NULL, Xi, Ci, 
                    deltaT, deltaD, Zmat.T, Zmat.D, copula.index, 
                    Gfun = Gfun, control = control)
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
      index = grep("dLambdaT", names(theta.est))
      dLambdaT.est = theta.est[index]
      Psi.dLambdaT = Psi[, index]
      LambdaT.est = sum.I(tk, ">=", tk, dLambdaT.est)
      Psi.LambdaT = tcrossprod(
        Psi.dLambdaT, 
        1 * (tk >= matrix(tk, byrow = T, 
                          nrow = length(tk), ncol = length(tk))))
      LambdaT.se = sqrt(diag(crossprod(Psi.LambdaT, Psi.LambdaT) / N) / N)  
      index = grep("dLambdaD", names(theta.est))
      dLambdaD.est = theta.est[index]
      Psi.dLambdaD = Psi[, index]
      LambdaD.est = sum.I(dk, ">=", dk, dLambdaD.est)
      Psi.LambdaD = tcrossprod(
        Psi.dLambdaD, 
        1 * (dk >= matrix(dk, byrow = T, 
                          nrow = length(dk), ncol = length(dk))))
      LambdaD.se = sqrt(diag(crossprod(Psi.LambdaD, Psi.LambdaD) / N) / N)  
      
      out.MLE =list(convergence = est.res$converged, 
                    time = usedtime,
                    iterations = est.res$iterations,
                    theta.ini = theta.ini, 
                    alpha = data.frame(est = alpha.est, se = alpha.se), 
                    betaT = data.frame(est = betaT.est, se = betaT.se),
                    betaD = data.frame(est = betaD.est, se = betaD.se),
                    LambdaT = data.frame(time = tk, est = LambdaT.est, 
                                         se = LambdaT.se), 
                    LambdaD = data.frame(time = dk, est = LambdaD.est,
                                         se = LambdaD.se))
      
      ## ===== PMLE ==== ##
      
      objfun<- function(x){
        f = lln.fun(theta = x, thetaT = NULL, thetaD = thetaD.est, 
                    Xi, Ci, deltaT, deltaD, Zmat.T, Zmat.D, copula.index, Gfun, 
                    control = control)
        g = dlln.fun(theta = x, thetaT = NULL, thetaD = thetaD.est, Xi, Ci, 
                     deltaT, deltaD, Zmat.T, Zmat.D, copula.index, Gfun, 
                     control = control)
        B = ddlln.fun(theta = x, thetaT = NULL, thetaD = thetaD.est, 
                      Xi, Ci, deltaT, deltaD, Zmat.T, Zmat.D, copula.index, Gfun, 
                      control = control)
        list(value = f, gradient = g, hessian = B)
      }
      usedtime <- system.time({
        est.res <- trust(
          objfun, theta.ini.0, 5, 100, iterlim = 300, 
          minimize= FALSE, blather = T)})
      theta.est = est.res$argument
      names(theta.est) = c(
        "copula", paste0("T.",var1), paste0("dLambdaT.", 1 : n.tk))
      Imat = - ddlln.fun(theta = theta.est, thetaT = NULL, thetaD = thetaD.est, 
                         Xi, Ci, deltaT, deltaD, 
                         Zmat.T, Zmat.D, copula.index, Gfun, control = control)
      
      dll.k = ddll.fun2(theta = theta.est, thetaD = thetaD.est, 
                        Xi, Ci, deltaT, deltaD, 
                        Zmat.T, Zmat.D, copula.index, Gfun, control = control) 
      dll.k[is.na(dll.k)] = 0
      
      Psi =  dll.fun(theta = theta.est, thetaT = NULL, thetaD = thetaD.est, 
                     Xi, Ci, deltaT, deltaD, 
                     Zmat.T, Zmat.D, copula.index, Gfun, control = control) 
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
                           paste0("D.",var2), paste0("dLambdaD.", 1 : n.dk))
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
      
      index = grep("dLambdaT", names(theta.est))
      dLambdaT.est = theta.est[index]
      Psi.dLambdaT = Psi[, index]
      LambdaT.est = sum.I(tk, ">=", tk, dLambdaT.est)
      Psi.LambdaT = tcrossprod(
        Psi.dLambdaT, 
        1 * (tk >= matrix(tk, byrow = T, 
                          nrow = length(tk), ncol = length(tk))))
      LambdaT.se = sqrt(diag(crossprod(Psi.LambdaT, Psi.LambdaT) / N) / N)  
      index = grep("dLambdaD", names(theta.est))
      dLambdaD.est = theta.est[index]
      Psi.dLambdaD = Psi[, index]
      LambdaD.est = sum.I(dk, ">=", dk, dLambdaD.est)
      Psi.LambdaD = tcrossprod(
        Psi.dLambdaD, 
        1 * (dk >= matrix(dk, byrow = T, 
                          nrow = length(dk), ncol = length(dk))))
      LambdaD.se = sqrt(diag(crossprod(Psi.LambdaD, Psi.LambdaD) / N) / N)  
      
      out.PMLE =list(convergence = est.res$converged, 
                     time = usedtime,
                     iterations = est.res$iterations,
                     theta.ini = theta.ini, 
                     alpha = data.frame(est = alpha.est, se = alpha.se), 
                     betaT = data.frame(est = betaT.est, se = betaT.se),
                     betaD = data.frame(est = betaD.est, se = betaD.se),
                     LambdaT = data.frame(time = tk, est = LambdaT.est, 
                                          se = LambdaT.se), 
                     LambdaD = data.frame(time = dk, est = LambdaD.est,
                                          se = LambdaD.se))
      save(out.MLE, out.PMLE, file = file.nm)
    }
  }
}

cl = makeCluster(4)
registerDoSNOW(cl)



# Simulation Study I ----
res.dir = "sim_results/"
family = "Gumbel";  nsamp = 200; tau = 0.4; theta.true = c(1, 0.5)
execute.fun(nrep = 2, nsamp = nsamp, tau = tau, 
            copula.true = family, copula.fam = family, 
            theta.true = theta.true, 
            seed1 = 20230709, res.dir=res.dir)

for (family in c("Clayton", "Gumbel")){
  for (nsamp in c(200, 400)){
    for (tau in c(0.4, 0.6, 0.8)){
      for (theta.true in list(c(1, 0.5), c(1, 1))){
        execute.fun(nrep = 2, nsamp = nsamp, tau = tau, 
                    copula.true = family, copula.fam = family, 
                    theta.true = theta.true, 
                    seed1 = 20230709, res.dir=res.dir)
      }
    }
  }
}

# Simulation Study II
execute.fun(nrep = 2, nsamp = 400, tau = 0.6, 
            copula.true = "Gumbel", copula.fam = "Clayton", 
            theta.true = c(1, 0.5), 
            seed1 = 20230709, res.dir = res.dir)
execute.fun(nrep = 2, nsamp = 400, tau = 0.6, 
            copula.true = "Clayton", copula.fam = "Gumbel", 
            theta.true = c(1, 0.5), 
            seed1 = 20230709, res.dir = res.dir)
stopCluster(cl)

