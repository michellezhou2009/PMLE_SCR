fitSPT = function(data, time = "time", status = "status", formula = ~ 1, 
                  Gfun = "PH", str = NULL){
  
  if (!rlang::is_formula(formula)) 
    stop("Argument \"formula\" needs to be a formula.")
  allvars = list(time, status); allnm = c("time", "status")
  tmp.index = sapply(allvars, function(x) is.character(x))
  if (any(tmp.index == FALSE)) 
    stop(
      paste0("Values of arguments: \"", 
             paste0(allnm[which(tmp.index == FALSE)], collapse = ", "),
             "\" needs to be characters."
      )
    )
  var.nm = all.vars(formula)
  allvars = c(time, status, var.nm)
  exists.index = allvars %in% colnames(data)
  if (any(exists.index == FALSE)) 
    stop(paste0("Variables: ", 
                paste0(allvars[exists.index == FALSE], collapse = ", ")),
         " do not exists in the data."
    )
  
  if (!is.null(str)){
    if (! is.character(str)) 
      stop("Value of argument \"str\" need to be characters")
    if (! str %in% colnames(data))
      stop(str, "does not exist in the data")
    id = data[, colnames(data) == str]
  } else id = seq(nrow(data))
  N = length(unique(id))  
  
  lln.fun <- function(theta, xi, di, zi, G.all){
    g.fun = G.all$g.fun; dg.fun = G.all$dg.fun;
    # N = length(xi); 
    n.b = ncol(zi); 
    tt = sort(unique(xi[di == 1])); n.tt = length(tt)
    b = theta[1 : n.b]; dLambda =  theta[-c(1 : n.b)]
    if (any(dLambda <= 0)) return(-Inf) else {
      bzi = as.vector(zi %*% b)
      log.dLambda = rep(0, length(xi)); 
      log.dLambda[di == 1] = log(dLambda[match(xi[di == 1], tt)])
      Lambda.xi = sum.I(xi, ">=", tt, dLambda)
      
      sum(di * bzi + log.dLambda + di * log(dg.fun(Lambda.xi * exp(bzi))) - 
            g.fun(Lambda.xi * exp(bzi)))
    }
  }
  dll.fun <- function(theta, xi, di, zi, G.all){
    g.fun = G.all$g.fun; dg.fun = G.all$dg.fun; ddg.fun = G.all$ddg.fun 
    # N = length(xi); 
    n.b = ncol(zi); 
    tt = sort(unique(xi[di == 1])); n.tt = length(tt)
    b = theta[1 : n.b]; dLambda =  theta[-c(1 : n.b)]
    if (any(dLambda <=0)) return(NA) else{
      bzi = as.vector(zi %*% b)
      Lambda.xi = sum.I(xi, ">=", tt, dLambda)
      xi.g.tt = (xi >= matrix(tt, byrow = T, nrow = length(xi), ncol = n.tt))
      dN.tt = (xi == matrix(tt, byrow = T, nrow = length(xi), ncol = n.tt)) * di
      dll.b = zi * (
        di + di * Lambda.xi * exp(bzi) * ddg.fun(Lambda.xi * exp(bzi)) /
          dg.fun(Lambda.xi * exp(bzi)) - dg.fun(Lambda.xi * exp(bzi)) * 
          Lambda.xi * exp(bzi)
      )
      dll.dLambda = (
        di * exp(bzi) * ddg.fun(Lambda.xi * exp(bzi)) /
          dg.fun(Lambda.xi * exp(bzi)) - 
          dg.fun(Lambda.xi * exp(bzi)) * exp(bzi)
      ) * xi.g.tt + 
        dN.tt / matrix(dLambda, byrow = T, nrow = length(xi), ncol = n.tt)
      return(cbind(dll.b, dll.dLambda))
    }
  }
  dlln.fun <- function(theta, xi, di, zi, G.all){
    dll = dll.fun(theta, xi, di, zi, G.all)
    if (any(is.na(dll))) return(rep(-Inf, length(theta))) else {
      # return(apply(dll, 2, mean))
      return(apply(dll, 2, sum))
    }
  }
  ddlln.fun <- function(theta, xi, di, zi, G.all){
    g.fun = G.all$g.fun; dg.fun = G.all$dg.fun; ddg.fun = G.all$ddg.fun 
    dddg.fun = G.all$dddg.fun
    # N = length(xi); 
    n.b = ncol(zi); 
    tt = sort(unique(xi[di == 1])); n.tt = length(tt)
    b = theta[1 : n.b]; dLambda =  theta[-c(1 : n.b)]
    if (any(dLambda <=0)) 
      return(matrix(-Inf, length(theta), length(theta))) else{
        bzi = as.vector(zi %*% b)
        Lambda.xi = sum.I(xi, ">=", tt, dLambda)
        xi.g.tt = (xi >= matrix(tt, byrow = T, nrow = length(xi), ncol = n.tt))
        dN.tt = (xi == matrix(tt, byrow = T, nrow = length(xi), ncol = n.tt)) * di
        dll.bb = crossprod(
          (di * Lambda.xi * exp(bzi) * ddg.fun(Lambda.xi * exp(bzi)) /
             dg.fun(Lambda.xi * exp(bzi)) + 
             di * Lambda.xi ^ 2 * exp(bzi) ^ 2 * 
             (dddg.fun(Lambda.xi * exp(bzi)) / dg.fun(Lambda.xi * exp(bzi)) - 
                (ddg.fun(Lambda.xi * exp(bzi)) / dg.fun(Lambda.xi * exp(bzi))) ^ 2)
           - ddg.fun(Lambda.xi * exp(bzi)) * Lambda.xi ^ 2 * exp(bzi) ^ 2 - 
             dg.fun(Lambda.xi * exp(bzi)) * Lambda.xi * exp(bzi)) * 
            zi, zi) 
        dll.dLambdadLambda = crossprod(
          (di * exp(bzi) ^ 2 *  
             (dddg.fun(Lambda.xi * exp(bzi)) / dg.fun(Lambda.xi * exp(bzi)) - 
                (ddg.fun(Lambda.xi * exp(bzi)) / dg.fun(Lambda.xi * exp(bzi))) ^ 2) -
             ddg.fun(Lambda.xi * exp(bzi)) * exp(bzi) ^ 2
          ) * xi.g.tt, xi.g.tt)  - diag(apply(dN.tt, 2, sum) / dLambda ^ 2)
        dll.dLambdab = crossprod(
          (di * exp(bzi) * 
             ddg.fun(Lambda.xi * exp(bzi)) / dg.fun(Lambda.xi * exp(bzi)) + 
             di * exp(bzi) ^ 2 * Lambda.xi * 
             (dddg.fun(Lambda.xi * exp(bzi)) / dg.fun(Lambda.xi * exp(bzi)) - 
                (ddg.fun(Lambda.xi * exp(bzi)) / dg.fun(Lambda.xi * exp(bzi))) ^ 2) -
             exp(bzi) * dg.fun(Lambda.xi * exp(bzi)) - 
             exp(bzi) ^ 2 * Lambda.xi * ddg.fun(Lambda.xi * exp(bzi))) * xi.g.tt
          , zi) 
        return(rbind(
          cbind(dll.bb, t(dll.dLambdab)),
          cbind(dll.dLambdab, dll.dLambdadLambda)
        ))
      }
  }
  
  xi = data[, colnames(data) == time]
  di = data[, colnames(data) == status]
  zi = model.matrix(formula, data)[, - 1, drop = F]
  n.b = ncol(zi)
  tt = sort(unique(xi[di == 1])); n.tt = length(tt)
  b.ini = rep(0, n.b); dLambda.ini = rep(0.1, n.tt)
  G.all = G.funs(Gfun)
  theta.ini  = c(b.ini, dLambda.ini)
  
  objfun<- function(x){
    f = lln.fun(theta = x, xi, di, zi, G.all)
    g = dlln.fun(theta = x,xi, di, zi, G.all)
    B = ddlln.fun(theta = x, xi, di, zi, G.all)
    list(value = f, gradient = g, hessian = B)
  }
  usedtime <- system.time({
    res <- trust::trust(objfun, theta.ini, 5, 100, iterlim = 300, 
                        minimize= FALSE, blather = T)})
  theta.est = res$argument
  b.est = theta.est[c(1 : n.b)]; dLambda.est = theta.est[ - c(1 : n.b)]
  Imat = - ddlln.fun(theta.est, xi, di, zi, G.all) / N
  dll = dll.fun(theta.est, xi, di, zi, G.all)
  if (N == length(xi)) 
    Psi.theta = dll %*% solve(Imat) 
  else
    Psi.theta = aggregate(. ~ id, data.frame(id = id, dll), sum) %>%
    select(- id) %>% as.matrix() %*% solve(Imat)
  theta.cov.model = solve(Imat) / N
  theta.cov.robust = crossprod(Psi.theta, Psi.theta) / (N ^ 2)
  beta.summary = data.frame(
    est = b.est, se = sqrt(diag(theta.cov.model)[c(1 : n.b)]), 
    robust = sqrt(diag(theta.cov.robust)[c(1 : n.b)])
  ); rownames(beta.summary) = colnames(zi)
  dLambda.summary = data.frame(
    time = tt, est = dLambda.est, 
    se = sqrt(diag(theta.cov.model)[ - c(1 : n.b)]), 
    robust = sqrt(diag(theta.cov.robust)[ - c(1 : n.b)])
  );
  Psi.beta = Psi.theta[ , 1 : n.b, drop = F]
  Psi.dLambda = Psi.theta[, - c(1 : n.b), drop = F]
  Psi.Lambda = t(apply(Psi.dLambda, 1, cumsum))
  
  tmp = matrix(1, nrow = n.tt, ncol = n.tt)
  tmp[upper.tri(tmp, diag = FALSE)] = 0
  Lambda.summary = data.frame(
    time = tt, est = cumsum(dLambda.est),
    se = 
      sqrt(diag(tmp %*% theta.cov.model[- c(1 : n.b), - c(1 : n.b)] %*% t(tmp))),
    robust = sqrt(apply(Psi.Lambda ^ 2, 2, sum)) / N
  )
  return(list(
    call = list(time = time, status = status, G.all = G.all, formula = formula),
    beta = beta.summary, dLambda = dLambda.summary, Lambda = Lambda.summary,
    varcov = list(model = theta.cov.model, robust = theta.cov.robust),
    Psi.theta = list(beta = Psi.beta, dLambda = Psi.dLambda, 
                     Lambda = Psi.Lambda),
    convergence = res$converged, niter = res$iterations))
}

predict.fitSPT= function(obj, newdata){
  time = obj$call$time; status = obj$call$status
  formula = obj$call$formula
  g.fun = obj$call$G.all$g.fun; dg.fun = obj$call$G.all$dg.fun
  b = as.vector(obj$beta$est)
  tt = as.vector(obj$dLambda$time)
  dLambda = as.vector(obj$dLambda$est)
  xk = newdata[ , colnames(newdata) == time]
  dk = newdata[ , colnames(newdata) == status]
  zk = model.matrix(formula, newdata)[ , -1, drop = F]
  bzk = as.vector(zk %*% b)
  Lambda.xk = sum.I(xk, ">=", tt, dLambda)
  uk = exp( - g.fun(Lambda.xk * exp(bzk)))
  xk.g.tt = 1 * (xk >= matrix(tt, byrow = T, 
                              nrow = length(xk), ncol = length(tt)))
  Psi.theta = cbind(obj$Psi.theta$beta, obj$Psi.theta$dLambda)
  A = cbind(zk * Lambda.xk * exp(bzk), xk.g.tt * exp(bzk)) *
    uk * dg.fun(Lambda.xk * exp(bzk))
  se.uk = sqrt(diag(A %*% obj$varcov$model %*% t(A)))
  Psi.uk = - Psi.theta %*% t(A)
  robust.se.uk = sqrt(apply(Psi.uk ^ 2, 2, sum)) / nrow(Psi.uk)
  
  list(surv = data.frame(fit = uk, se = se.uk, robust = robust.se.uk), 
       Psi.surv = Psi.uk)
}
