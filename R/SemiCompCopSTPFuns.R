# theta does not include copula parameters
lln.fun1 <- function(copula.para, theta, thetaD = NULL, Xi, Ci, deltaT, deltaD, 
                     Zmat.T, Zmat.D, copula.index, 
                     Gfun = list(T = "PH", D = "PH")){
  N = length(Xi); n.theta = length(theta)
  n.bT = ncol(Zmat.T); n.bD = ncol(Zmat.D)
  tk = sort(unique(Xi[deltaT == 1])); n.tk = length(unique(tk))
  dk = sort(unique(Ci[deltaD == 1])); n.dk = length(unique(dk)) 
  copula.para = copula.para;
  betaT = theta[ c(1 : n.bT)]; dLambdaT = theta[n.bT + c(1 : n.tk)]
  if (is.null(thetaD)){
    betaD = theta[n.bT + n.tk + c(1 : n.bD)]
    dLambdaD = theta[ n.bT + n.tk + n.bD + c(1 : n.dk)]
  } else {
    betaD = thetaD[c(1 : n.bD)]
    dLambdaD = thetaD[n.bD + c(1 : n.dk)]
  }
  if (any(dLambdaT < 0) | any(dLambdaD < 0)) ll = - Inf else {
    G.all = G.funs(Gfun$T)
    gT.fun = G.all$g.fun; dgT.fun = G.all$dg.fun; ddgT.fun = G.all$ddg.fun
    dddgT.fun = G.all$dddg.fun
    G.all = G.funs(Gfun$D)
    gD.fun = G.all$g.fun; dgD.fun = G.all$dg.fun; ddgD.fun = G.all$ddg.fun
    dddgD.fun = G.all$dddg.fun
    
    LambdaT.Xi = sum.I(Xi, ">=", tk, dLambdaT)
    LambdaD.Ci = sum.I(Ci, ">=", dk, dLambdaD)
    bZ.T = as.vector(Zmat.T %*% betaT); ebZ.T = exp(bZ.T)
    bZ.D = as.vector(Zmat.D %*% betaD); ebZ.D = exp(bZ.D) 
    S.Xi = exp(- gT.fun(LambdaT.Xi * ebZ.T))
    S.Ci = exp(- gD.fun(LambdaD.Ci * ebZ.D))
    if (copula.index == 4) { # gumbel
      S.Xi[S.Xi == 1] = exp(- 1 / N)
      S.Ci[S.Ci == 1] = exp(- 1 / N)
    }
    log.dLambda.Xi = log.dLambda.Ci = rep(0, N); 
    log.dLambda.Xi[deltaT == 1] = log(dLambdaT[match(Xi[deltaT == 1], tk)])
    log.dLambda.Ci[deltaD == 1] = log(dLambdaD[match(Ci[deltaD == 1], dk)])
    ll = funsData(ll_funs, d1 = deltaT, d2 = deltaD, u1 = S.Xi, u2 = S.Ci,
                  copula.index = copula.index, para = copula.para) + 
      deltaT * (- gT.fun(LambdaT.Xi * ebZ.T) + 
                  log(dgT.fun(LambdaT.Xi * ebZ.T)) + bZ.T + log.dLambda.Xi) +
      deltaD * (- gD.fun(LambdaD.Ci * ebZ.D) +
                  log(dgD.fun(LambdaD.Ci * ebZ.D))+ bZ.D + log.dLambda.Ci) 
   sum(ll, na.rm = T) / N 
  }
}
# theta includes copula parameter
lln.fun <- function(theta, thetaD = NULL, Xi, Ci, deltaT, deltaD, 
                    Zmat.T, Zmat.D, copula.index, 
                    Gfun = list(T = "PH", D = "PH")){
  N = length(Xi); n.theta = length(theta)
  n.bT = ncol(Zmat.T); n.bD = ncol(Zmat.D)
  tk = sort(unique(Xi[deltaT == 1])); n.tk = length(unique(tk))
  dk = sort(unique(Ci[deltaD == 1])); n.dk = length(unique(dk)) 
  copula.para = theta[1];
  betaT = theta[1 + c(1 : n.bT)]; dLambdaT = theta[1 + n.bT + c(1 : n.tk)]
  if (is.null(thetaD)){
    betaD = theta[1 + n.bT + n.tk + c(1 : n.bD)]
    dLambdaD = theta[1 + n.bT + n.tk + n.bD + c(1 : n.dk)]
  } else {
    betaD = thetaD[c(1 : n.bD)]
    dLambdaD = thetaD[n.bD + c(1 : n.dk)]
  }
  if (any(dLambdaT < 0) | any(dLambdaD < 0)) ll = - Inf else {
    G.all = G.funs(Gfun$T)
    gT.fun = G.all$g.fun; dgT.fun = G.all$dg.fun; ddgT.fun = G.all$ddg.fun
    dddgT.fun = G.all$dddg.fun
    G.all = G.funs(Gfun$D)
    gD.fun = G.all$g.fun; dgD.fun = G.all$dg.fun; ddgD.fun = G.all$ddg.fun
    dddgD.fun = G.all$dddg.fun
    
    LambdaT.Xi = sum.I(Xi, ">=", tk, dLambdaT)
    LambdaD.Ci = sum.I(Ci, ">=", dk, dLambdaD)
    bZ.T = as.vector(Zmat.T %*% betaT); ebZ.T = exp(bZ.T)
    bZ.D = as.vector(Zmat.D %*% betaD); ebZ.D = exp(bZ.D) 
    S.Xi = exp(- gT.fun(LambdaT.Xi * ebZ.T))
    S.Ci = exp(- gD.fun(LambdaD.Ci * ebZ.D))
    if (copula.index == 4) { # gumbel
      S.Xi[S.Xi == 1] = exp(- 1 / N)
      S.Ci[S.Ci == 1] = exp(- 1 / N)
    }
    log.dLambda.Xi = log.dLambda.Ci = rep(0, N); 
    log.dLambda.Xi[deltaT == 1] = log(dLambdaT[match(Xi[deltaT == 1], tk)])
    log.dLambda.Ci[deltaD == 1] = log(dLambdaD[match(Ci[deltaD == 1], dk)])
    ll = funsData(ll_funs, d1 = deltaT, d2 = deltaD, u1 = S.Xi, u2 = S.Ci,
                  copula.index = copula.index, para = copula.para) + 
      deltaT * (- gT.fun(LambdaT.Xi * ebZ.T) + 
                  log(dgT.fun(LambdaT.Xi * ebZ.T)) + bZ.T + log.dLambda.Xi) +
      deltaD * (- gD.fun(LambdaD.Ci * ebZ.D) +
                  log(dgD.fun(LambdaD.Ci * ebZ.D))+ bZ.D + log.dLambda.Ci) 
    sum(ll, na.rm = T)/N 
  }
}

dll.fun <- function(theta, thetaD = NULL, Xi, Ci, deltaT, deltaD, 
                    Zmat.T, Zmat.D, copula.index, 
                    Gfun = list(T = "PH", D = "PH"), keep = F){
  N = length(Xi); n.theta = length(theta)
  n.bT = ncol(Zmat.T); n.bD = ncol(Zmat.D)
  tk = sort(unique(Xi[deltaT == 1])); n.tk = length(unique(tk))
  dk = sort(unique(Ci[deltaD == 1])); n.dk = length(unique(dk)) # sakie
  copula.para = theta[1];
  betaT = theta[1 + c(1 : n.bT)]; dLambdaT = theta[1 + n.bT + c(1 : n.tk)]
  if (is.null(thetaD)){
    betaD = theta[1 + n.bT + n.tk + c(1 : n.bD)]
    dLambdaD = theta[1 + n.bT + n.tk + n.bD + c(1 : n.dk)]
  } else {
    betaD = thetaD[c(1 : n.bD)]
    dLambdaD = thetaD[n.bD + c(1 : n.dk)]
  }

  if (any(dLambdaT < 0) | any(dLambdaD < 0))
    return(NA) else{
      G.all = G.funs(Gfun$T)
      gT.fun = G.all$g.fun; dgT.fun = G.all$dg.fun; ddgT.fun = G.all$ddg.fun
      dddgT.fun = G.all$dddg.fun
      G.all = G.funs(Gfun$D)
      gD.fun = G.all$g.fun; dgD.fun = G.all$dg.fun; ddgD.fun = G.all$ddg.fun
      dddgD.fun = G.all$dddg.fun

      LambdaT.Xi = sum.I(Xi, ">=", tk, dLambdaT)
      LambdaD.Ci = sum.I(Ci, ">=", dk, dLambdaD)
      bZ.T = as.vector(Zmat.T %*% betaT); ebZ.T = exp(bZ.T)
      bZ.D = as.vector(Zmat.D %*% betaD); ebZ.D = exp(bZ.D)
      S.Xi = exp(- gT.fun(LambdaT.Xi * ebZ.T))
      S.Ci = exp(- gD.fun(LambdaD.Ci * ebZ.D))
      if (copula.index == 4) { # gumbel
        S.Xi[S.Xi == 1] = exp(- 1 / N)
        S.Ci[S.Ci == 1] = exp(- 1 / N)
      }
      Xi.g.tk = 1 * (Xi >= matrix(tk, byrow = T, nrow = N, ncol = n.tk))
      dN.tk = 1*(Xi == matrix(tk, byrow = T, nrow = N, ncol = n.tk)) * deltaT

      dll.u1 = funsData(
        Funs = dll.u1_funs, d1 = deltaT, d2 = deltaD, u1 = S.Xi, u2 = S.Ci,
        copula.index = copula.index, para = copula.para)
      dll.para = funsData(
        Funs = dll.para_funs, d1 = deltaT, d2 = deltaD, u1 = S.Xi, u2 = S.Ci,
        copula.index = copula.index, para = copula.para)
      dll.betaT =
        (- dll.u1 * S.Xi * LambdaT.Xi * ebZ.T * dgT.fun(LambdaT.Xi * ebZ.T) -
           deltaT * LambdaT.Xi * ebZ.T * dgT.fun(LambdaT.Xi * ebZ.T) +
           deltaT * LambdaT.Xi * ebZ.T * ddgT.fun(LambdaT.Xi * ebZ.T) /
           dgT.fun(LambdaT.Xi * ebZ.T) + deltaT
        ) * Zmat.T
      dll.dLambdaT =
        - dll.u1 * S.Xi  * ebZ.T * Xi.g.tk * dgT.fun(LambdaT.Xi * ebZ.T) -
        deltaT * ebZ.T * Xi.g.tk * dgT.fun(LambdaT.Xi * ebZ.T) +
        deltaT * ebZ.T * Xi.g.tk * ddgT.fun(LambdaT.Xi * ebZ.T) /
        dgT.fun(LambdaT.Xi * ebZ.T) +
        dN.tk / matrix(dLambdaT, byrow = T, nrow = N, ncol = n.tk)
      out = cbind(dll.para, dll.betaT, dll.dLambdaT)
      if (is.null(thetaD)){
        Ci.g.dk = 1 * (Ci >= matrix(dk, byrow = T, nrow = N, ncol = n.dk))
        dN.dk = 1*(Ci == matrix(dk, byrow = T, nrow = N, ncol = n.dk)) * deltaD
        dll.u2 = funsData(
          Funs = dll.u2_funs, d1 = deltaT, d2 = deltaD, u1 = S.Xi, u2 = S.Ci,
          copula.index = copula.index, para = copula.para)
        dll.dLambdaD =
          - dll.u2 * S.Ci * ebZ.D * Ci.g.dk * dgD.fun(LambdaD.Ci * ebZ.D) -
          deltaD * ebZ.D * Ci.g.dk * dgD.fun(LambdaD.Ci * ebZ.D) +
          deltaD * ebZ.D * Ci.g.dk * ddgD.fun(LambdaD.Ci * ebZ.D) /
          dgD.fun(LambdaD.Ci * ebZ.D) +
          dN.dk / matrix(dLambdaD, byrow = T, nrow = N, ncol = n.dk)
        dll.betaD =
          (- dll.u2 * S.Ci * LambdaD.Ci * ebZ.D * dgD.fun(LambdaD.Ci * ebZ.D) -
             deltaD * LambdaD.Ci * ebZ.D * dgD.fun(LambdaD.Ci * ebZ.D) +
             deltaD * LambdaD.Ci * ebZ.D * ddgD.fun(LambdaD.Ci * ebZ.D) /
             dgD.fun(LambdaD.Ci * ebZ.D) + deltaD) * Zmat.D
        out = cbind(out, dll.betaD, dll.dLambdaD)
      }
      if (!keep) out = out[complete.cases(out), ] else out 
      return(out)
    }
}

dlln.fun <- function(theta, thetaD = NULL, Xi, Ci, deltaT, deltaD, 
                     Zmat.T, Zmat.D, copula.index, 
                     Gfun = list(T = "PH", D = "PH")){
  dll = dll.fun(theta, thetaD, Xi, Ci, deltaT, deltaD, 
                Zmat.T, Zmat.D, copula.index, 
                Gfun = Gfun)
  if (any(is.na(dll))) return(rep(-Inf, length(theta))) else
    return(apply(dll, 2, mean, na.rm = T))
}

ddlln.fun <- function(theta, thetaD = NULL, Xi, Ci, deltaT, deltaD, 
                      Zmat.T, Zmat.D, copula.index, 
                      Gfun = list(T = "PH", D = "PH")){
  N = length(Xi); n.theta = length(theta)
  n.bT = ncol(Zmat.T); n.bD = ncol(Zmat.D)
  tk = sort(unique(Xi[deltaT == 1])); n.tk = length(unique(tk))
  dk = sort(unique(Ci[deltaD == 1])); n.dk = length(unique(dk)) # sakie
  copula.para = theta[1];
  betaT = theta[1 + c(1 : n.bT)]; dLambdaT = theta[1 + n.bT + c(1 : n.tk)]
  if (is.null(thetaD)){
    betaD = theta[1 + n.bT + n.tk + c(1 : n.bD)]
    dLambdaD = theta[1 + n.bT + n.tk + n.bD + c(1 : n.dk)]
  } else {
    betaD = thetaD[c(1 : n.bD)]
    dLambdaD = thetaD[n.bD + c(1 : n.dk)]
  }
  
  if (any(dLambdaT < 0) | any(dLambdaD < 0)) 
    return(matrix(-Inf, ncol = n.theta, nrow = n.theta)) else{
      G.all = G.funs(Gfun$T)
      gT.fun = G.all$g.fun; dgT.fun = G.all$dg.fun; ddgT.fun = G.all$ddg.fun
      dddgT.fun = G.all$dddg.fun
      G.all = G.funs(Gfun$D)
      gD.fun = G.all$g.fun; dgD.fun = G.all$dg.fun; ddgD.fun = G.all$ddg.fun
      dddgD.fun = G.all$dddg.fun
      
      LambdaT.Xi = sum.I(Xi, ">=", tk, dLambdaT)
      LambdaD.Ci = sum.I(Ci, ">=", dk, dLambdaD)
      bZ.T = as.vector(Zmat.T %*% betaT); ebZ.T = exp(bZ.T)
      bZ.D = as.vector(Zmat.D %*% betaD); ebZ.D = exp(bZ.D) 
      S.Xi = exp(- gT.fun(LambdaT.Xi * ebZ.T))
      S.Ci = exp(- gD.fun(LambdaD.Ci * ebZ.D))
      if (copula.index == 4) { # gumbel
        S.Xi[S.Xi == 1] = exp(- 1 / N)
        S.Ci[S.Ci == 1] = exp(- 1 / N)
      }
      Xi.g.tk = 1 * (Xi >= matrix(tk, byrow = T, nrow = N, ncol = n.tk)) 
      dN.tk = 1*(Xi == matrix(tk, byrow = T, nrow = N, ncol = n.tk)) * deltaT
      
      dll.u1 = funsData(
        Funs = dll.u1_funs, d1 = deltaT, d2 = deltaD, u1 = S.Xi, u2 = S.Ci, 
        copula.index = copula.index, para = copula.para) 
      dll.para = funsData(
        Funs = dll.para_funs, d1 = deltaT, d2 = deltaD, u1 = S.Xi, u2 = S.Ci, 
        copula.index = copula.index, para = copula.para)
      dll.u1u1 = funsData(
        Funs = dll.u1u1_funs, d1 = deltaT, d2 = deltaD, u1 = S.Xi, u2 = S.Ci, 
        copula.index = copula.index, para = copula.para) 
      dll.para2 = mean(funsData(
        Funs = dll.para2_funs, d1 = deltaT, d2 = deltaD, u1 = S.Xi, u2 = S.Ci, 
        copula.index = copula.index, para = copula.para), na.rm = T)
      dll.u1para = funsData(
        Funs = dll.u1para_funs, d1 = deltaT, d2 = deltaD, u1 = S.Xi, u2 = S.Ci, 
        copula.index = copula.index, para = copula.para) 
      
      dll.dLambdaTdLambdaT = crossprod(
        (dll.u1u1 * S.Xi ^ 2 * ebZ.T ^ 2 * dgT.fun(LambdaT.Xi * ebZ.T) ^ 2 + 
           dll.u1 * S.Xi * ebZ.T ^ 2 * dgT.fun(LambdaT.Xi * ebZ.T) ^ 2 -
           dll.u1 * S.Xi * ebZ.T ^ 2 * ddgT.fun(LambdaT.Xi * ebZ.T) -
           deltaT * ebZ.T ^ 2 * ddgT.fun(LambdaT.Xi * ebZ.T) +
           deltaT * ebZ.T ^ 2 * 
           (dddgT.fun(LambdaT.Xi * ebZ.T) / dgT.fun(LambdaT.Xi * ebZ.T) - 
              (ddgT.fun(LambdaT.Xi * ebZ.T) / dgT.fun(LambdaT.Xi * ebZ.T)) ^ 2
           )) * 
          Xi.g.tk,  Xi.g.tk) / N - diag(apply(dN.tk, 2, mean) / dLambdaT ^ 2)
      
      dll.betaTbetaT = crossprod(
        (dll.u1u1 * S.Xi ^ 2 * LambdaT.Xi ^ 2 * ebZ.T ^ 2 * 
           dgT.fun(LambdaT.Xi * ebZ.T) ^ 2 +
           dll.u1 * S.Xi * LambdaT.Xi ^ 2 * ebZ.T ^ 2 * 
           dgT.fun(LambdaT.Xi * ebZ.T) ^ 2  -
           dll.u1 * S.Xi * LambdaT.Xi * ebZ.T * dgT.fun(LambdaT.Xi * ebZ.T) - 
           dll.u1 * S.Xi * LambdaT.Xi ^2  * ebZ.T ^2 * 
           ddgT.fun(LambdaT.Xi * ebZ.T) -
           deltaT * LambdaT.Xi * ebZ.T * dgT.fun(LambdaT.Xi * ebZ.T) - 
           deltaT * LambdaT.Xi ^2  * ebZ.T ^2  * 
           ddgT.fun(LambdaT.Xi * ebZ.T) +
           deltaT * LambdaT.Xi * ebZ.T * ddgT.fun(LambdaT.Xi * ebZ.T) /
           dgT.fun(LambdaT.Xi * ebZ.T) +
           deltaT * LambdaT.Xi ^ 2 * ebZ.T ^ 2 * 
           (dddgT.fun(LambdaT.Xi * ebZ.T) / dgT.fun(LambdaT.Xi * ebZ.T) - 
              (ddgT.fun(LambdaT.Xi * ebZ.T) / dgT.fun(LambdaT.Xi * ebZ.T)) ^ 2
           )) * Zmat.T, Zmat.T
      ) / N
      
      dll.dLambdaTbetaT = crossprod(
        (dll.u1u1 * S.Xi ^ 2 * ebZ.T ^ 2 * LambdaT.Xi * 
           dgT.fun(LambdaT.Xi * ebZ.T) ^ 2 + 
           dll.u1 * S.Xi * ebZ.T ^ 2 * LambdaT.Xi * 
           dgT.fun(LambdaT.Xi * ebZ.T) ^ 2 - 
           dll.u1 * S.Xi * ebZ.T * dgT.fun(LambdaT.Xi * ebZ.T) -
           dll.u1 * S.Xi * ebZ.T ^ 2 * LambdaT.Xi * 
           ddgT.fun(LambdaT.Xi * ebZ.T)  - 
           deltaT * ebZ.T * dgT.fun(LambdaT.Xi * ebZ.T) - 
           deltaT * ebZ.T ^ 2 * LambdaT.Xi * ddgT.fun(LambdaT.Xi * ebZ.T) + 
           deltaT * ebZ.T * ddgT.fun(LambdaT.Xi * ebZ.T) /
           dgT.fun(LambdaT.Xi * ebZ.T) +
           deltaT * ebZ.T ^ 2 *  LambdaT.Xi * 
           (dddgT.fun(LambdaT.Xi * ebZ.T) / dgT.fun(LambdaT.Xi * ebZ.T) - 
              (ddgT.fun(LambdaT.Xi * ebZ.T) / dgT.fun(LambdaT.Xi * ebZ.T)) ^ 2
           )) * Xi.g.tk, Zmat.T) / N
      
      dll.dLambdaTpara = apply(- dll.u1para * S.Xi * ebZ.T * Xi.g.tk * 
                                 dgT.fun(LambdaT.Xi * ebZ.T), 2, mean, na.rm = T)
      
      dll.betaTpara = apply(
        - dll.u1para * S.Xi * ebZ.T * LambdaT.Xi * Zmat.T * 
          dgT.fun(LambdaT.Xi * ebZ.T), 2, mean, na.rm = T)
      
      out = rbind(
        c(dll.para2, dll.betaTpara, dll.dLambdaTpara),
        cbind(dll.betaTpara, dll.betaTbetaT, t(dll.dLambdaTbetaT)),
        cbind(dll.dLambdaTpara, dll.dLambdaTbetaT, dll.dLambdaTdLambdaT)
      )
      
      if (is.null(thetaD)){
        Ci.g.dk = 1 * (Ci >= matrix(dk, byrow = T, nrow = N, ncol = n.dk)) 
        dN.dk = 1*(Ci == matrix(dk, byrow = T, nrow = N, ncol = n.dk)) * deltaD
        dll.u2 = funsData(
          Funs = dll.u2_funs, d1 = deltaT, d2 = deltaD, u1 = S.Xi, u2 = S.Ci, 
          copula.index = copula.index, para = copula.para)
        dll.u2u2 = funsData(
          Funs = dll.u2u2_funs, d1 = deltaT, d2 = deltaD, u1 = S.Xi, u2 = S.Ci, 
          copula.index = copula.index, para = copula.para)
        dll.u1u2 = funsData(
          Funs = dll.u1u2_funs, d1 = deltaT, d2 = deltaD, u1 = S.Xi, u2 = S.Ci, 
          copula.index = copula.index, para = copula.para) 
        dll.u2para = funsData(
          Funs = dll.u2para_funs, d1 = deltaT, d2 = deltaD, u1 = S.Xi, u2 = S.Ci, 
          copula.index = copula.index, para = copula.para) 
        
        dll.dLambdaDdLambdaD = crossprod(
          (dll.u2u2 * S.Ci ^ 2 * ebZ.D ^ 2 * dgD.fun(LambdaD.Ci * ebZ.D) ^ 2 + 
             dll.u2 * S.Ci * ebZ.D ^ 2 * dgD.fun(LambdaD.Ci * ebZ.D) ^ 2 -
             dll.u2 * S.Ci * ebZ.D ^ 2 * ddgD.fun(LambdaD.Ci * ebZ.D) -
             deltaD * ebZ.D ^ 2 * ddgD.fun(LambdaD.Ci * ebZ.D) +
             deltaD * ebZ.D ^ 2 * 
             (dddgD.fun(LambdaD.Ci * ebZ.D) / dgD.fun(LambdaD.Ci * ebZ.D) - 
                (ddgD.fun(LambdaD.Ci * ebZ.D) / dgD.fun(LambdaD.Ci * ebZ.D)) ^ 2
             )) * 
            Ci.g.dk,  Ci.g.dk) / N - diag(apply(dN.dk, 2, mean) / dLambdaD ^ 2)
        
        dll.betaDbetaD = crossprod(
          (dll.u2u2 * S.Ci ^ 2 * LambdaD.Ci ^ 2 * ebZ.D ^ 2 * 
             dgD.fun(LambdaD.Ci * ebZ.D) ^ 2 +
             dll.u2 * S.Ci * LambdaD.Ci ^ 2 * ebZ.D ^ 2 * 
             dgD.fun(LambdaD.Ci * ebZ.D) ^ 2  -
             dll.u2 * S.Ci * LambdaD.Ci * ebZ.D * dgD.fun(LambdaD.Ci * ebZ.D) - 
             dll.u2 * S.Ci * LambdaD.Ci ^2  * ebZ.D ^2 * 
             ddgD.fun(LambdaD.Ci * ebZ.D) -
             deltaD * LambdaD.Ci * ebZ.D * dgD.fun(LambdaD.Ci * ebZ.D) - 
             deltaD * LambdaD.Ci ^2  * ebZ.D ^2  * 
             ddgD.fun(LambdaD.Ci * ebZ.D) +
             deltaD * LambdaD.Ci * ebZ.D * ddgD.fun(LambdaD.Ci * ebZ.D) /
             dgD.fun(LambdaD.Ci * ebZ.D) +
             deltaD * LambdaD.Ci ^ 2 * ebZ.D ^ 2 * 
             (dddgD.fun(LambdaD.Ci * ebZ.D) / dgD.fun(LambdaD.Ci * ebZ.D) - 
                (ddgD.fun(LambdaD.Ci * ebZ.D) / dgD.fun(LambdaD.Ci * ebZ.D)) ^ 2
             )) * Zmat.D, Zmat.D
        ) / N
        
        dll.dLambdaDbetaD = crossprod(
          (dll.u2u2 * S.Ci ^ 2 * ebZ.D ^ 2 * LambdaD.Ci * 
             dgD.fun(LambdaD.Ci * ebZ.D) ^ 2 + 
             dll.u2 * S.Ci * ebZ.D ^ 2 * LambdaD.Ci * 
             dgD.fun(LambdaD.Ci * ebZ.D) ^ 2 - 
             dll.u2 * S.Ci * ebZ.D * dgD.fun(LambdaD.Ci * ebZ.D) -
             dll.u2 * S.Ci * ebZ.D ^ 2 * LambdaD.Ci * 
             ddgD.fun(LambdaD.Ci * ebZ.D)  - 
             deltaD * ebZ.D * dgD.fun(LambdaD.Ci * ebZ.D) - 
             deltaD * ebZ.D ^ 2 * LambdaD.Ci * ddgD.fun(LambdaD.Ci * ebZ.D) + 
             deltaD * ebZ.D * ddgD.fun(LambdaD.Ci * ebZ.D) /
             dgD.fun(LambdaD.Ci * ebZ.D) +
             deltaD * ebZ.D ^ 2 *  LambdaD.Ci * 
             (dddgD.fun(LambdaD.Ci * ebZ.D) / dgD.fun(LambdaD.Ci * ebZ.D) - 
                (ddgD.fun(LambdaD.Ci * ebZ.D) / dgD.fun(LambdaD.Ci * ebZ.D)) ^ 2
             )) * Ci.g.dk, Zmat.D) / N
        
        dll.dLambdaDpara = apply(
          - dll.u2para * S.Ci * ebZ.D * Ci.g.dk * dgD.fun(LambdaD.Ci * ebZ.D), 
          2, mean, na.rm = T)
        
        dll.betaDpara = apply(
          - dll.u2para * S.Ci * ebZ.D * LambdaD.Ci * 
            dgD.fun(LambdaD.Ci * ebZ.D) * Zmat.D, 2, mean, na.rm = T) 
        
        dll.dLambdaTdLambdaD = crossprod(
          dll.u1u2 * S.Xi * S.Ci * ebZ.T * ebZ.D *
            dgT.fun(LambdaT.Xi * ebZ.T) * dgD.fun(LambdaD.Ci * ebZ.D) *
            Xi.g.tk, Ci.g.dk) / N
        
        dll.dLambdaDbetaT = crossprod(
          dll.u1u2 * S.Xi * S.Ci * ebZ.T * ebZ.D * LambdaT.Xi * 
            dgT.fun(LambdaT.Xi * ebZ.T) * dgD.fun(LambdaD.Ci * ebZ.D) *
            Ci.g.dk, Zmat.T) / N
        
        dll.dLambdaTbetaD = crossprod(
          dll.u1u2 * S.Xi * S.Ci * ebZ.T * ebZ.D * LambdaD.Ci * 
            dgT.fun(LambdaT.Xi * ebZ.T) * dgD.fun(LambdaD.Ci * ebZ.D) *
            Xi.g.tk, Zmat.D) / N
        
        dll.betaTbetaD = crossprod(
          dll.u1u2 * S.Xi * S.Ci * LambdaT.Xi * LambdaD.Ci * 
            ebZ.T * ebZ.D * dgT.fun(LambdaT.Xi * ebZ.T) * 
            dgD.fun(LambdaD.Ci * ebZ.D) * 
            Zmat.T, Zmat.D) / N
        
        out22 = rbind(
          cbind(dll.betaDbetaD, t(dll.dLambdaDbetaD)),
          cbind(dll.dLambdaDbetaD, dll.dLambdaDdLambdaD)
        )
        out12 = rbind(
          c(dll.betaDpara, dll.dLambdaDpara),
          cbind(dll.betaTbetaD, t(dll.dLambdaDbetaT)),
          cbind(dll.dLambdaTbetaD, dll.dLambdaTdLambdaD)
        )
        out = rbind(cbind(out, out12), cbind(t(out12), out22))
      }
      out = out[complete.cases(out),] # sakie
      return(out)
    }
}

ddll.fun2 <- function(theta, thetaD, Xi, Ci, deltaT, deltaD, 
                      Zmat.T, Zmat.D, copula.index, 
                      Gfun = list(T = "PH", D = "PH")){
  N = length(Xi); n.theta = length(theta)
  n.bT = ncol(Zmat.T); n.bD = ncol(Zmat.D)
  tk = sort(unique(Xi[deltaT == 1])); n.tk = length(unique(tk))
  dk = sort(unique(Ci[deltaD == 1])); n.dk = length(unique(dk)) # sakie
  copula.para = theta[1];
  betaT = theta[1 + c(1 : n.bT)]; dLambdaT = theta[1 + n.bT + c(1 : n.tk)]
  betaD = thetaD[c(1 : n.bD)]
  dLambdaD = thetaD[n.bD + c(1 : n.dk)]
  G.all = G.funs(Gfun$T)
  gT.fun = G.all$g.fun; dgT.fun = G.all$dg.fun; ddgT.fun = G.all$ddg.fun
  dddgT.fun = G.all$dddg.fun
  G.all = G.funs(Gfun$D)
  gD.fun = G.all$g.fun; dgD.fun = G.all$dg.fun; ddgD.fun = G.all$ddg.fun
  dddgD.fun = G.all$dddg.fun
  
  LambdaT.Xi = sum.I(Xi, ">=", tk, dLambdaT)
  LambdaD.Ci = sum.I(Ci, ">=", dk, dLambdaD)
  bZ.T = as.vector(Zmat.T %*% betaT); ebZ.T = exp(bZ.T)
  bZ.D = as.vector(Zmat.D %*% betaD); ebZ.D = exp(bZ.D) 
  S.Xi = exp(- gT.fun(LambdaT.Xi * ebZ.T))
  S.Ci = exp(- gD.fun(LambdaD.Ci * ebZ.D))
  if (copula.index == 4) { # gumbel
    S.Xi[S.Xi == 1] = exp(- 1 / N)
    S.Ci[S.Ci == 1] = exp(-1 / N)
  }
  Xi.g.tk = 1 * (Xi >= matrix(tk, byrow = T, nrow = N, ncol = n.tk))
  dN.tk = 1*(Xi == matrix(tk, byrow = T, nrow = N, ncol = n.tk)) * deltaT
  
  dll.u2para = funsData(
    Funs = dll.u2para_funs, d1 = deltaT, d2 = deltaD, u1 = S.Xi, u2 = S.Ci, 
    copula.index = copula.index, para = copula.para) 
  dll.u1u2 = funsData(
    Funs = dll.u1u2_funs, d1 = deltaT, d2 = deltaD, u1 = S.Xi, u2 = S.Ci, 
    copula.index = copula.index, para = copula.para) 
  dll.betaTu2 = - dll.u1u2 * S.Xi * LambdaT.Xi * ebZ.T  * 
    dgT.fun(LambdaT.Xi * ebZ.T) * Zmat.T
  dll.dLambdaTu2 = - dll.u1u2 * S.Xi  * ebZ.T * Xi.g.tk * 
    dgT.fun(LambdaT.Xi * ebZ.T) 
  out = cbind(dll.u2para, dll.betaTu2, dll.dLambdaTu2)
  return(out)
  
}