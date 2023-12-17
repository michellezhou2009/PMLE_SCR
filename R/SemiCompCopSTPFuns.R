lln.fun <- function(theta, thetaT = NULL, thetaD = NULL, 
                    Xi, Ci, deltaT, deltaD, 
                    Zmat.T, Zmat.D, copula.index, 
                    Gfun = list(T = "PH", D = "PH"),
                    copula.link = 
                      list(
                        h.fun = function(x) {x}, 
                        dot.h.fun = function(x) {rep(1, length(x))},
                        ddot.h.fun = function(x) {rep(0, length(x))}), 
                    Wmat = NULL, 
                    control = list(copula.lwr = 0, copula.upr = 28)){
  N = length(Xi); n.theta = length(theta)
  n.bT = ncol(Zmat.T); n.bD = ncol(Zmat.D)
  tk = sort(unique(Xi[deltaT == 1])); n.tk = length(unique(tk))
  dk = sort(unique(Ci[deltaD == 1])); n.dk = length(unique(dk)) 

  h.fun = copula.link$h.fun
  dot.h.fun = copula.link$dot.h.fun
  ddot.h.fun = copula.link$ddot.h.fun
  copula.lwr = control$copula.lwr; copula.upr = control$copula.upr
  if (is.null(Wmat)) Wmat = matrix(1, ncol = 1, nrow = N)
  n.para = ncol(Wmat); copula.para = theta[1 : n.para]
  copula.lp = as.vector(Wmat %*% matrix(copula.para, byrow = F, ncol = 1))
  copula.para = h.fun(copula.lp)
  theta = theta[- c(1 : n.para)]
  if (is.null(thetaT)){
    betaT = theta[c(1 : n.bT)]; 
    dLambdaT = theta[n.bT + c(1 : n.tk)]  
    theta = theta[- c(1 : (n.bT + n.tk))]
  } else {
    betaT = thetaT[c(1 : n.bT)]
    dLambdaT = thetaT[n.bT + c(1 : n.tk)]
  }
  if (is.null(thetaD)){
    betaD = theta[c(1 : n.bD)]
    dLambdaD = theta[n.bD + c(1 : n.dk)]
  } else {
    betaD = thetaD[c(1 : n.bD)]
    dLambdaD = thetaD[n.bD + c(1 : n.dk)]
  }
  
  yes.constraint = any(dLambdaT < 0) | 
    any(dLambdaD < 0) | 
    min(copula.para) < copula.lwr | 
    max(copula.para) > copula.upr

  if (yes.constraint) ll = - Inf else {
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

dll.fun <- function(theta, thetaT = NULL, thetaD = NULL, 
                    Xi, Ci, deltaT, deltaD, 
                    Zmat.T, Zmat.D, copula.index, 
                    Gfun = list(T = "PH", D = "PH"),
                    copula.link = 
                      list(
                        h.fun = function(x) {x}, 
                        dot.h.fun = function(x) {rep(1, length(x))},
                        ddot.h.fun = function(x) {rep(0, length(x))}), 
                    Wmat = NULL, 
                    control = list(copula.lwr = 0, copula.upr = 28)){
  N = length(Xi); n.theta = length(theta)
  n.bT = ncol(Zmat.T); n.bD = ncol(Zmat.D)
  tk = sort(unique(Xi[deltaT == 1])); n.tk = length(unique(tk))
  dk = sort(unique(Ci[deltaD == 1])); n.dk = length(unique(dk)) 
  
  h.fun = copula.link$h.fun
  dot.h.fun = copula.link$dot.h.fun
  ddot.h.fun = copula.link$ddot.h.fun
  copula.lwr = control$copula.lwr; copula.upr = control$copula.upr
  if (is.null(Wmat)) Wmat = matrix(1, ncol = 1, nrow = N)
  n.para = ncol(Wmat); copula.para = theta[1 : n.para]
  copula.lp = as.vector(Wmat %*% matrix(copula.para, byrow = F, ncol = 1))
  copula.para = h.fun(copula.lp)
  theta = theta[- c(1 : n.para)]
  if (is.null(thetaT)){
    betaT = theta[c(1 : n.bT)]; 
    dLambdaT = theta[n.bT + c(1 : n.tk)]  
    theta = theta[- c(1 : (n.bT + n.tk))]
  } else {
    betaT = thetaT[c(1 : n.bT)]
    dLambdaT = thetaT[n.bT + c(1 : n.tk)]
  }
  if (is.null(thetaD)){
    betaD = theta[c(1 : n.bD)]
    dLambdaD = theta[n.bD + c(1 : n.dk)]
  } else {
    betaD = thetaD[c(1 : n.bD)]
    dLambdaD = thetaD[n.bD + c(1 : n.dk)]
  }
  
  yes.constraint = any(dLambdaT < 0) | 
    any(dLambdaD < 0) | 
    min(copula.para) < copula.lwr | 
    max(copula.para) > copula.upr
  if (yes.constraint)
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
      
      dll.para = funsData(
        Funs = dll.para_funs, d1 = deltaT, d2 = deltaD, u1 = S.Xi, u2 = S.Ci,
        copula.index = copula.index, para = copula.para)
      dll.para = dll.para * dot.h.fun(copula.lp) * Wmat
      # if (n.para == 1) dll.para = matrix(dll.para, byrow = F, ncol = 1)
      out = dll.para  
      if (is.null(thetaT)){
        dll.u1 = funsData(
          Funs = dll.u1_funs, d1 = deltaT, d2 = deltaD, u1 = S.Xi, u2 = S.Ci,
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
        out = cbind(out, dll.betaT, dll.dLambdaT)
      }
      
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
      #if (!keep) out = out[complete.cases(out), ] else out 
      return(out)
    }
}

dlln.fun <- function(theta, thetaT = NULL, thetaD = NULL, 
                     Xi, Ci, deltaT, deltaD, 
                     Zmat.T, Zmat.D, copula.index, 
                     Gfun = list(T = "PH", D = "PH"),
                     copula.link = 
                       list(
                         h.fun = function(x) {x}, 
                         dot.h.fun = function(x) {rep(1, length(x))},
                         ddot.h.fun = function(x) {rep(0, length(x))}), 
                     Wmat = NULL, 
                     control = list(copula.lwr = 0, copula.upr = 28)){
  
  dll = dll.fun(theta, thetaT, thetaD, Xi, Ci, deltaT, deltaD, 
                Zmat.T, Zmat.D, copula.index, 
                Gfun = Gfun, copula.link = copula.link, 
                Wmat = Wmat, control = control)
  if (any(is.na(dll))) return(rep(-Inf, length(theta))) else
    return(apply(dll, 2, mean, na.rm = T))
}

ddlln.fun <- function(theta, thetaT = NULL, thetaD = NULL, 
                      Xi, Ci, deltaT, deltaD, 
                      Zmat.T, Zmat.D, copula.index, 
                      Gfun = list(T = "PH", D = "PH"),
                      copula.link = 
                        list(
                          h.fun = function(x) {x}, 
                          dot.h.fun = function(x) {rep(1, length(x))},
                          ddot.h.fun = function(x) {rep(0, length(x))}), 
                      Wmat = NULL, 
                      control = list(copula.lwr = 0, copula.upr = 28)){
  N = length(Xi); n.theta = length(theta)
  n.bT = ncol(Zmat.T); n.bD = ncol(Zmat.D)
  tk = sort(unique(Xi[deltaT == 1])); n.tk = length(unique(tk))
  dk = sort(unique(Ci[deltaD == 1])); n.dk = length(unique(dk)) 
  
  h.fun = copula.link$h.fun
  dot.h.fun = copula.link$dot.h.fun
  ddot.h.fun = copula.link$ddot.h.fun
  copula.lwr = control$copula.lwr; copula.upr = control$copula.upr
  if (is.null(Wmat)) Wmat = matrix(1, ncol = 1, nrow = N)
  n.para = ncol(Wmat); copula.para = theta[1 : n.para]
  copula.lp = as.vector(Wmat %*% matrix(copula.para, byrow = F, ncol = 1))
  copula.para = h.fun(copula.lp)
  theta = theta[- c(1 : n.para)]
  if (is.null(thetaT)){
    betaT = theta[c(1 : n.bT)]; 
    dLambdaT = theta[n.bT + c(1 : n.tk)]  
    theta = theta[- c(1 : (n.bT + n.tk))]
  } else {
    betaT = thetaT[c(1 : n.bT)]
    dLambdaT = thetaT[n.bT + c(1 : n.tk)]
  }
  if (is.null(thetaD)){
    betaD = theta[c(1 : n.bD)]
    dLambdaD = theta[n.bD + c(1 : n.dk)]
  } else {
    betaD = thetaD[c(1 : n.bD)]
    dLambdaD = thetaD[n.bD + c(1 : n.dk)]
  }
  
  
  yes.constraint = any(dLambdaT < 0) | 
    any(dLambdaD < 0) | 
    min(copula.para) < copula.lwr | 
    max(copula.para) > copula.upr
  
  
  if (yes.constraint)
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
      
      dll.para = funsData(
        Funs = dll.para_funs, d1 = deltaT, d2 = deltaD, u1 = S.Xi, u2 = S.Ci, 
        copula.index = copula.index, para = copula.para)
      dll.para2 = funsData(
        Funs = dll.para2_funs, d1 = deltaT, d2 = deltaD, u1 = S.Xi, u2 = S.Ci, 
        copula.index = copula.index, para = copula.para)
      dll.para2 = matrix(
          apply((dll.para2 * (dot.h.fun(copula.lp) ^ 2) + 
                  dll.para * ddot.h.fun(copula.lp)) * 
                  Wmat[, rep(1 : n.para, n.para), drop = F] * 
                  Wmat[, rep(1 : n.para, each = n.para), drop = F], 
                2, mean, na.rm = T),
          byrow = T, nrow = n.para, ncol = n.para)
      out = dll.para2
      if (is.null(thetaT)){
        dll.u1 = funsData(
          Funs = dll.u1_funs, d1 = deltaT, d2 = deltaD, u1 = S.Xi, u2 = S.Ci, 
          copula.index = copula.index, para = copula.para) 
        dll.u1u1 = funsData(
          Funs = dll.u1u1_funs, d1 = deltaT, d2 = deltaD, u1 = S.Xi, u2 = S.Ci, 
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
        dll.u1para = funsData(
          Funs = dll.u1para_funs, d1 = deltaT, d2 = deltaD, u1 = S.Xi, u2 = S.Ci, 
          copula.index = copula.index, para = copula.para) 
        dll.dLambdaTpara = - dll.u1para * S.Xi * ebZ.T * Xi.g.tk * 
          dgT.fun(LambdaT.Xi * ebZ.T)
        dll.betaTpara = - dll.u1para * S.Xi * ebZ.T * LambdaT.Xi * Zmat.T * 
          dgT.fun(LambdaT.Xi * ebZ.T)
        dll.dLambdaTpara = matrix(
            apply(dll.dLambdaTpara[, rep(1 : n.tk, n.para)] *
                    dot.h.fun(copula.lp) * Wmat[, rep(1 : n.para, each = n.tk)],
                  2, mean, na.rm = T),
            byrow = T, ncol = n.tk, nrow = n.para) 
          dll.betaTpara = matrix(
            apply(dll.betaTpara[, rep(1 : n.bT, n.para)] * 
                    dot.h.fun(copula.lp) * Wmat[, rep(1 : n.para, each = n.bT)],
                  2, mean, na.rm = T),
            byrow = T, ncol = n.bT, nrow = n.para)
        out = rbind(
          cbind(out, dll.betaTpara, dll.dLambdaTpara),
          cbind(t(dll.betaTpara), dll.betaTbetaT, t(dll.dLambdaTbetaT)),
          cbind(t(dll.dLambdaTpara), dll.dLambdaTbetaT, dll.dLambdaTdLambdaT)
        )
      }
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
        
        dll.dLambdaDpara = 
          - dll.u2para * S.Ci * ebZ.D * Ci.g.dk * dgD.fun(LambdaD.Ci * ebZ.D)
        
        dll.betaDpara = 
          - dll.u2para * S.Ci * ebZ.D * LambdaD.Ci * 
          dgD.fun(LambdaD.Ci * ebZ.D) * Zmat.D 
        
      
          dll.dLambdaDpara = matrix(
            apply(dll.dLambdaDpara[, rep(1 : n.dk, n.para)] * 
                    dot.h.fun(copula.lp) * Wmat[, rep(1 : n.para, each = n.dk)],
                  2, mean, na.rm = T), byrow = T, nrow = n.para, ncol = n.dk)
          dll.betaDpara = matrix(
            apply(dll.betaDpara[, rep(1 : n.bD, n.para)] * 
                    dot.h.fun(copula.lp) * Wmat[, rep(1 : n.para, each = n.bD)], 
                  2, mean, na.rm = T), byrow = T, nrow = n.para, ncol = n.bD)
        
        out22 = rbind(
          cbind(dll.betaDbetaD, t(dll.dLambdaDbetaD)),
          cbind(dll.dLambdaDbetaD, dll.dLambdaDdLambdaD)
        )
        out12 = rbind(
          cbind(dll.betaDpara, dll.dLambdaDpara),
          cbind(dll.betaTbetaD, t(dll.dLambdaDbetaT)),
          cbind(dll.dLambdaTbetaD, dll.dLambdaTdLambdaD)
        )
        out = rbind(cbind(out, out12), cbind(t(out12), out22))
      }
      return(out)
    }
}

ddll.fun2 <- function(theta, thetaD, Xi, Ci, deltaT, deltaD, 
                      Zmat.T, Zmat.D, copula.index, 
                      Gfun = list(T = "PH", D = "PH"),
                      copula.link = 
                        list(
                          h.fun = function(x) {x}, 
                          dot.h.fun = function(x) {rep(1, length(x))},
                          ddot.h.fun = function(x) {rep(0, length(x))}), 
                      Wmat = NULL, 
                      control = list(copula.lwr = 0, copula.upr = 28)){
  N = length(Xi); n.theta = length(theta)
  n.bT = ncol(Zmat.T); n.bD = ncol(Zmat.D)
  tk = sort(unique(Xi[deltaT == 1])); n.tk = length(unique(tk))
  dk = sort(unique(Ci[deltaD == 1])); n.dk = length(unique(dk)) 
  
  h.fun = copula.link$h.fun
  dot.h.fun = copula.link$dot.h.fun
  ddot.h.fun = copula.link$ddot.h.fun
  copula.lwr = control$copula.lwr; copula.upr = control$copula.upr
  if (is.null(Wmat)) Wmat = matrix(1, ncol = 1, nrow = N)
  n.para = ncol(Wmat); copula.para = theta[1 : n.para]
  copula.lp = as.vector(Wmat %*% matrix(copula.para, byrow = F, ncol = 1))
  copula.para = h.fun(copula.lp)
  theta = theta[- c(1 : n.para)]
  
  betaT = theta[c(1 : n.bT)]; 
  dLambdaT = theta[n.bT + c(1 : n.tk)]  
  theta = theta[- c(1 : (n.bT + n.tk))]

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
  
  dll.u2para = dll.u2para * dot.h.fun(copula.lp) * Wmat
  out = cbind(dll.u2para, dll.betaTu2, dll.dLambdaTu2)
  return(out)
  
}