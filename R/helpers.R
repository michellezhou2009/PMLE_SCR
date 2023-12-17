sum.I <- function(yy, FUN, Yi, Vi=NULL){
  if (FUN == "<" | FUN == ">=") {yy = -yy; Yi = -Yi}
  # for each distinct ordered failure time t[j], number of Xi < t[j]
  pos <- rank(c(yy, Yi),ties.method = 'f')[1 : length(yy)] - 
    rank(yy, ties.method = 'f')    
  if (substring(FUN, 2, 2) == "=") pos = length(Yi) - pos # number of Xi>= t[j]
  if (!is.null(Vi)) {
    ## if FUN contains '=', tmpind is the order of decending
    if(substring(FUN, 2, 2) == "=") tmpind = order(-Yi) else  tmpind = order(Yi)
    Vi <- apply(as.matrix(Vi)[tmpind, , drop=F], 2, cumsum)
    return(rbind(0, Vi)[pos + 1, ])
  } else return(pos)
}

MyTau2Par <- function(copula.fam, tau){
  switch(copula.fam,
         "Gumbel" = {simu.index = 14},
         "Clayton" = {simu.index = 13}
  )
  BiCopTau2Par(family = simu.index, tau = tau)
}
MyCopula = function(copula.index){
  switch(
    as.character(copula.index),
    "3" = { # Clayton
      C.exp = expression((u1^(-para) + u2^(-para) -1)^(-1/para))
    },
    "4" = { # Gumbel
      C.exp = expression(exp(-((-log(u1))^(para) + (-log(u2))^(para))^(1/para)))
    },
    "5" = { # Frank
      C.exp = 
        expression((-1/para) * log(1 + ((exp(-para*u1)-1)*(exp(-para*u2)-1))/(exp(-para)-1)))
    }
  )
  C.fun = function(u1, u2, para){};
  body(C.fun) = C.exp
  dC.para = function(u1, u2, para){};
  body(dC.para) = D(C.exp, "para")
  dC.para2 = function(u1, u2, para){};
  body(dC.para2) = D(D(C.exp, "para"), "para")
  dC.u1 = function(u1, u2, para){};
  body(dC.u1) = D(C.exp, "u1")
  dC.u1u1 = function(u1, u2, para){};
  body(dC.u1u1) = D(D(C.exp, "u1"),"u1")
  dC.u1para = function(u1, u2, para){};
  body(dC.u1para) = D(D(C.exp, "u1"), "para")
  dC.u1u1para = function(u1, u2, para){};
  body(dC.u1u1para) = D(D(D(C.exp, "u1"), "u1"), "para")
  dC.u1para2 = function(u1, u2, para){};
  body(dC.u1para2) = D(D(D(C.exp, "u1"), "para"), "para")
  dC.u2 = function(u1, u2, para){};
  body(dC.u2) = D(C.exp, "u2")
  dC.u2u2 = function(u1, u2, para){};
  body(dC.u2u2) = D(D(C.exp, "u2"), "u2")
  dC.u2para = function(u1, u2, para){};
  body(dC.u2para) = D(D(C.exp, "u2"), "para")
  dC.u2u2para = function(u1, u2, para){};
  body(dC.u2u2para) = D(D(D(C.exp, "u2"), "u2"), "para")
  dC.u2para2 = function(u1, u2, para){};
  body(dC.u2para2) = D(D(D(C.exp, "u2"), "para"), "para")
  dC.u1u2 = function(u1, u2, para){};
  body(dC.u1u2)= D(D(C.exp, "u1"), "u2")
  dC.u1u2para = function(u1, u2, para){};
  body(dC.u1u2para) = D(D(D(C.exp, "u1"), "u2"), "para")
  dC.u1u1u1 = function(u1, u2, para){};
  body(dC.u1u1u1) = D(D(D(C.exp, "u1"), "u1"), "u1")
  dC.u2u2u2 = function(u1, u2, para){};
  body(dC.u2u2u2) = D(D(D(C.exp, "u2"), "u2"), "u2")
  dC.u1u1u2 = function(u1, u2, para){};
  body(dC.u1u1u2) = D(D(D(C.exp, "u1"), "u1"), "u2")
  dC.u1u2u2 = function(u1, u2, para){};
  body(dC.u1u2u2) = D(D(D(C.exp, "u1"), "u2"), "u2")
  dC.u1u1u1u2 = function(u1, u2, para){};
  body(dC.u1u1u1u2) = D(D(D(D(C.exp, "u1"), "u1"), "u1"), "u2")
  dC.u1u2u2u2 = function(u1, u2, para){};
  body(dC.u1u2u2u2) = D(D(D(D(C.exp, "u1"), "u2"), "u2"), "u2")
  dC.u1u1u2u2 = function(u1, u2, para){};
  body(dC.u1u1u2u2) = D(D(D(D(C.exp, "u1"), "u1"), "u2"), "u2")
  dC.u1u1u2para = function(u1, u2, para){};
  body(dC.u1u1u2para) = D(D(D(D(C.exp, "u1"), "u1"), "u2"), "para")
  dC.u1u2u2para = function(u1, u2, para){};
  body(dC.u1u2u2para) = D(D(D(D(C.exp, "u1"), "u2"), "u2"), "para")
  dC.u1u2para2 = function(u1, u2, para){};
  body(dC.u1u2para2) = D(D(D(D(C.exp, "u1"), "u2"), "para"), "para")
  list(C.fun = C.fun, dC.para = dC.para, dC.para2 =dC.para2,
       dC.u1 = dC.u1, dC.u1u1 = dC.u1u1, dC.u1para = dC.u1para,
       dC.u1u1para = dC.u1u1para, dC.u1para2 = dC.u1para2,
       dC.u2 = dC.u2, dC.u2u2 = dC.u2u2, dC.u2para = dC.u2para,
       dC.u2u2para = dC.u2u2para, dC.u2para2 = dC.u2para2,
       dC.u1u2 = dC.u1u2, dC.u1u2para = dC.u1u2para,
       dC.u1u1u1 = dC.u1u1u1, dC.u2u2u2 = dC.u2u2u2,
       dC.u1u1u2 = dC.u1u1u2, dC.u1u2u2 = dC.u1u2u2,
       dC.u1u1u1u2 = dC.u1u1u1u2, dC.u1u2u2u2 = dC.u1u2u2u2,
       dC.u1u1u2u2 = dC.u1u1u2u2,
       dC.u1u1u2para = dC.u1u1u2para, dC.u1u2u2para = dC.u1u2u2para,
       dC.u1u2para2 = dC.u1u2para2)
  
}

ll_funs <- function(copula.index){
  ll = c();
  ll[[1]] = function(u1, u2, para) {
    log(BiCopCDF(u1, u2, family = copula.index, par = para)) # (0, 0)
  } # (0, 0)
  ll[[2]] = function(u1, u2, para) {
    log(BiCopHfunc1(u1, u2, family = copula.index, par = para)) # (1, 0)
  } # (1, 0)
  ll[[3]] = function(u1, u2, para) {
    log(BiCopHfunc2(u1, u2, family = copula.index, par = para)) # (0, 1)
  } # (0, 1)
  ll[[4]] = function(u1, u2, para) {
    log(BiCopPDF(u1, u2, family = copula.index, par = para))# (1, 1)
  } # (1, 1)
  return(ll)
}


dll.para_funs <- function(copula.index){
  dll.para = c()
  dll.para[[1]] = function(u1, u2, para) {
    MyCopula(copula.index)$dC.para(u1, u2, para) /
      BiCopCDF(u1, u2, family = copula.index, par = para)
  } # (0, 0)
  dll.para[[2]] = function(u1, u2, para) {
    BiCopHfuncDeriv(u2, u1, deriv = "par", family = copula.index, par = para) /
      BiCopHfunc1(u1, u2, family = copula.index, par = para)
  } # (1, 0)
  dll.para[[3]] = function(u1, u2, para) {
    BiCopHfuncDeriv(u1, u2, deriv = "par", family = copula.index, par = para) /
      BiCopHfunc2(u1, u2, family = copula.index, par = para)
  } # (1, 0)
  dll.para[[4]] = function(u1, u2, para) {
    BiCopDeriv(u1, u2,  deriv = "par", family = copula.index, par = para) /
      BiCopPDF(u1, u2, family = copula.index, par = para)
  } # (1, 1)
  return(dll.para)
}

dll.u1_funs <- function(copula.index){
  dll.para = c()
  dll.para[[1]] = function(u1, u2, para) {
    BiCopHfunc1(u1, u2, family = copula.index, par = para) /
      BiCopCDF(u1, u2, family = copula.index, par = para)
  } # (0, 0)
  dll.para[[2]] = function(u1, u2, para) {
    BiCopHfuncDeriv(u2, u1, deriv = "u2", family = copula.index, par = para) /
      BiCopHfunc1(u1, u2, family = copula.index, par = para)
  } # (1, 0)
  dll.para[[3]] = function(u1, u2, para) {
    BiCopPDF(u1, u2, family = copula.index, par = para) /
      BiCopHfunc2(u1, u2, family = copula.index, par = para)
  } # (0, 1)
  dll.para[[4]] = function(u1, u2, para) {
    BiCopDeriv(u1, u2,  deriv = "u1", family = copula.index, par = para) /
      BiCopPDF(u1, u2, family = copula.index, par = para)
  } # (1, 1)
  return(dll.para)
}

dll.u2_funs <- function(copula.index){
  dll.para = c()
  dll.para[[1]] = function(u1, u2, para) {
    BiCopHfunc2(u1, u2, family = copula.index, par = para) /
      BiCopCDF(u1, u2, family = copula.index, par = para)
  } # (0, 0)
  dll.para[[2]] = function(u1, u2, para) {
    BiCopPDF(u1, u2, family = copula.index, par = para) /
      BiCopHfunc1(u1, u2, family = copula.index, par = para)
  } # (1, 0)
  dll.para[[3]] = function(u1, u2, para) {
    BiCopHfuncDeriv(u1, u2, deriv = "u2", family = copula.index, par = para) /
      BiCopHfunc2(u1, u2, family = copula.index, par = para)
  } # (0, 1)
  dll.para[[4]] = function(u1, u2, para) {
    BiCopDeriv(u1, u2,  deriv = "u2", family = copula.index, par = para) /
      BiCopPDF(u1, u2, family = copula.index, par = para)
  } # (1, 1)
  return(dll.para)
}

dll.para2_funs <- function(copula.index){
  dll.para2 = c()
  dll.para = dll.para_funs(copula.index)
  dll.para2[[1]] = function(u1, u2, para) {
    MyCopula(copula.index)$dC.para2(u1, u2, para) /
      BiCopCDF(u1, u2, family = copula.index, par = para) -
      (dll.para[[1]](u1, u2, para))^2
  }
  dll.para2[[2]] = function(u1, u2, para) {
    BiCopHfuncDeriv2(u2, u1, deriv = "par", family = copula.index, par = para) /
      BiCopHfunc1(u1, u2, family = copula.index, par = para) -
      (dll.para[[2]](u1, u2, para))^2
  }
  dll.para2[[3]] = function(u1, u2, para) {
    BiCopHfuncDeriv2(u1, u2, deriv = "par", family = copula.index, par = para) /
      BiCopHfunc2(u1, u2, family = copula.index, par = para) -
      (dll.para[[3]](u1, u2, para))^2
  }
  dll.para2[[4]] = function(u1, u2, para) {
    BiCopDeriv2(u1, u2,  deriv = "par", family = copula.index, par = para) /
      BiCopPDF(u1, u2, family = copula.index, par = para) -
      (dll.para[[4]](u1, u2, para))^2
  }
  return(dll.para2)
}
dll.u1u1_funs <- function(copula.index){
  dll.para = c()
  dll.para[[1]] = function(u1, u2, para) {
    BiCopHfuncDeriv(u2, u1, deriv = "u2", 
                    family = copula.index, par = para) /
      BiCopCDF(u1, u2, family = copula.index, par = para) - 
      (BiCopHfunc1(u1, u2, family = copula.index, par = para) /
         BiCopCDF(u1, u2, family = copula.index, par = para)) ^ 2 
  } # (0, 0)
  dll.para[[2]] = function(u1, u2, para) {
    BiCopHfuncDeriv2(u2, u1, deriv = "u2", family = copula.index, par = para) /
      BiCopHfunc1(u1, u2, family = copula.index, par = para) -
      (BiCopHfuncDeriv(u2, u1, deriv = "u2", 
                       family = copula.index, par = para) / 
         BiCopHfunc1(u1, u2, family = copula.index, par = para)) ^ 2
  } # (1, 0)
  dll.para[[3]] = function(u1, u2, para) {
    BiCopDeriv(u1, u2, deriv = "u1", family = copula.index, par = para) /
      BiCopHfunc2(u1, u2, family = copula.index, par = para) - 
      (BiCopPDF(u1, u2, family = copula.index, par = para) / 
         BiCopHfunc2(u1, u2, family = copula.index, par = para)) ^ 2
  } # (0, 1)
  dll.para[[4]] = function(u1, u2, para) {
    BiCopDeriv2(u1, u2,  deriv = "u1", family = copula.index, par = para) /
      BiCopPDF(u1, u2, family = copula.index, par = para) -
      (BiCopDeriv(u1, u2, deriv = "u1", family = copula.index, par = para) /
         BiCopPDF(u1, u2, family = copula.index, par = para)) ^ 2
  } # (1, 1)
  return(dll.para)
}

dll.u2u2_funs <- function(copula.index){
  dll.para = c()
  dll.para[[1]] = function(u1, u2, para) {
    BiCopHfuncDeriv(u1, u2, deriv = "u2", 
                    family = copula.index, par = para) /
      BiCopCDF(u1, u2, family = copula.index, par = para) - 
      (BiCopHfunc2(u1, u2, family = copula.index, par = para) /
         BiCopCDF(u1, u2, family = copula.index, par = para)) ^ 2 
  } # (0, 0)
  dll.para[[2]] = function(u1, u2, para) {
    BiCopDeriv(u1, u2, deriv = "u2", family = copula.index, par = para) /
      BiCopHfunc1(u1, u2, family = copula.index, par = para) -
      (BiCopPDF(u1, u2, family = copula.index, par = para) / 
         BiCopHfunc1(u1, u2, family = copula.index, par = para)) ^ 2
  } # (1, 0)
  dll.para[[3]] = function(u1, u2, para) {
    BiCopHfuncDeriv2(u1, u2, deriv = "u2", family = copula.index, par = para)/
      BiCopHfunc2(u1, u2, family = copula.index, par = para) - 
      (BiCopHfuncDeriv(u1, u2, deriv = "u2", 
                       family = copula.index, par = para) / 
         BiCopHfunc2(u1, u2, family = copula.index, par = para)) ^ 2
  } # (0, 1)
  dll.para[[4]] = function(u1, u2, para) {
    BiCopDeriv2(u1, u2,  deriv = "u2", family = copula.index, par = para) /
      BiCopPDF(u1, u2, family = copula.index, par = para) -
      (BiCopDeriv(u1, u2, deriv = "u2", family = copula.index, par = para) /
         BiCopPDF(u1, u2, family = copula.index, par = para)) ^ 2
  } # (1, 1)
  return(dll.para)
}

dll.u1u2_funs <- function(copula.index){
  dll.para = c()
  dll.para[[1]] = function(u1, u2, para) {
    BiCopPDF(u1, u2, family = copula.index, par = para) /
      BiCopCDF(u1, u2, family = copula.index, par = para) - 
      BiCopHfunc1(u1, u2, family = copula.index, par = para) * 
      BiCopHfunc2(u1, u2, family = copula.index, par = para) /
      BiCopCDF(u1, u2, family = copula.index, par = para) ^ 2
  } # (0, 0)
  dll.para[[2]] = function(u1, u2, para) {
    BiCopDeriv(u1, u2, deriv = "u1", family = copula.index, par = para) /
      BiCopHfunc1(u1, u2, family = copula.index, par = para) -
      BiCopPDF(u1, u2, family = copula.index, par = para) * 
      BiCopHfuncDeriv(u2, u1, deriv = "u2",
                      family = copula.index, par = para)/ 
      BiCopHfunc1(u1, u2, family = copula.index, par = para) ^ 2
  } # (1, 0)
  dll.para[[3]] = function(u1, u2, para) {
    BiCopDeriv(u1, u2, deriv = "u2", family = copula.index, par = para)/
      BiCopHfunc2(u1, u2, family = copula.index, par = para) - 
      BiCopHfuncDeriv(u1, u2, deriv = "u2", 
                      family = copula.index, par = para) * 
      BiCopPDF(u1, u2, family = copula.index, par = para) / 
      BiCopHfunc2(u1, u2, family = copula.index, par = para) ^ 2
  } # (0, 1)
  dll.para[[4]] = function(u1, u2, para) {
    MyCopula(copula.index)$dC.u1u1u2u2(u1, u2, para)/
      BiCopPDF(u1, u2, family = copula.index, par = para) -
      BiCopDeriv(u1, u2, deriv = "u1", family = copula.index, par = para) *
      BiCopDeriv(u1, u2, deriv = "u2", family = copula.index, par = para)/
      BiCopPDF(u1, u2, family = copula.index, par = para) ^ 2
  } # (1, 1)
  return(dll.para)
}

dll.u1para_funs <- function(copula.index){
  dll.u1para = c()
  dll.u1para[[1]] = function(u1, u2, para) {
    BiCopHfuncDeriv(u2, u1, deriv = "par", family = copula.index, par = para) /
      BiCopCDF(u1, u2, family = copula.index, par = para) -
      MyCopula(copula.index)$dC.para(u1, u2, para) *
      BiCopHfunc1(u1, u2, family = copula.index, par = para) /
      BiCopCDF(u1, u2, family = copula.index, par = para) ^ 2
  }
  dll.u1para[[2]] = function(u1, u2, para) {
    BiCopHfuncDeriv2(u2, u1, deriv = "par1u2",
                     family = copula.index, par = para) /
      BiCopHfunc1(u1, u2, family = copula.index, par = para) -
      BiCopHfuncDeriv(u2, u1, deriv = "par",
                      family = copula.index, par = para) *
      BiCopHfuncDeriv(u2, u1, deriv = "u2",
                      family = copula.index, par = para)/
      BiCopHfunc1(u1, u2, family = copula.index, par = para) ^ 2
  }
  dll.u1para[[3]] = function(u1, u2, para) {
    BiCopDeriv(u1, u2, deriv = "par", family = copula.index, par = para) /
      BiCopHfunc2(u1, u2, family = copula.index, par = para) -
      BiCopHfuncDeriv(u1, u2, deriv = "par", family = copula.index, par = para) *
      BiCopPDF(u1, u2, family = copula.index, par = para)/
      BiCopHfunc2(u1, u2, family = copula.index, par = para) ^ 2
  }
  dll.u1para[[4]] = function(u1, u2, para) {
    BiCopDeriv2(u1, u2, deriv = "par1u1", family = copula.index, par = para) /
      BiCopPDF(u1, u2, family = copula.index, par = para) -
      BiCopDeriv(u1, u2,  deriv = "par", family = copula.index, par = para) *
      BiCopDeriv(u1, u2, deriv = "u1",
                 family = copula.index, par = para)/
      BiCopPDF(u1, u2, family = copula.index, par = para) ^ 2
  }
  return(dll.u1para)
}

dll.u2para_funs <- function(copula.index){
  dll.u2para = c()
  dll.u2para[[1]] = function(u1, u2, para) {
    BiCopHfuncDeriv(u1, u2, deriv = "par", family = copula.index, par = para) /
      BiCopCDF(u1, u2, family = copula.index, par = para) -
      MyCopula(copula.index)$dC.para(u1, u2, para) *
      BiCopHfunc2(u1, u2, family = copula.index, par = para) /
      BiCopCDF(u1, u2, family = copula.index, par = para) ^ 2
  }
  dll.u2para[[2]] = function(u1, u2, para) {
    BiCopDeriv(u1, u2, deriv = "par", family = copula.index, par = para) /
      BiCopHfunc1(u1, u2, family = copula.index, par = para) -
      BiCopHfuncDeriv(u2, u1, deriv = "par",
                      family = copula.index, par = para) *
      BiCopPDF(u1, u2, family = copula.index, par = para)/
      BiCopHfunc1(u1, u2, family = copula.index, par = para) ^ 2
  }
  dll.u2para[[3]] = function(u1, u2, para) {
    BiCopHfuncDeriv2(u1, u2, deriv = "par1u2",
                     family = copula.index, par = para) /
      BiCopHfunc2(u1, u2, family = copula.index, par = para) -
      BiCopHfuncDeriv(u1, u2, deriv = "par",
                      family = copula.index, par = para) *
      BiCopHfuncDeriv(u1, u2, deriv = "u2",
                      family = copula.index, par = para)/
      BiCopHfunc2(u1, u2, family = copula.index, par = para) ^ 2
  }
  dll.u2para[[4]] = function(u1, u2, para) {
    BiCopDeriv2(u1, u2, deriv = "par1u2", family = copula.index, par = para) /
      BiCopPDF(u1, u2, family = copula.index, par = para) -
      BiCopDeriv(u1, u2,  deriv = "par", family = copula.index, par = para) *
      BiCopDeriv(u1, u2, deriv = "u2",
                 family = copula.index, par = para)/
      BiCopPDF(u1, u2, family = copula.index, par = para) ^ 2
  }
  return(dll.u2para)
}



funsData <- function(Funs, d1, d2, u1, u2, copula.index, para){
  funs = Funs(copula.index)
  grp = d1 + 2 * d2 + 1
  out = lapply(1:4, function(k){
    index = which(grp == k) 
    if (sum(index) == 0) return(NULL) else {
      if (length(para) == 1) {
        ll = funs[[k]](u1 = u1[index], u2 = u2[index], para = para)
      } else {
        ll = funs[[k]](u1 = u1[index], u2 = u2[index], para = para[index])
      }
      data.frame(index = index, ll = ll)
    }
  }) %>% do.call(rbind, .)
  out = out %>% arrange(index, )
  out$ll
}
G.funs = function(Gfun){
  switch(Gfun, 
         "PH" = {
           g.fun = function(x){x}
           dg.fun = function(x){rep(1, length(x))}
           ddg.fun = function(x){rep(0, length(x))}
           dddg.fun = function(x){rep(0, length(x))}
         },
         "PO" = {
           g.fun = function(x){log(1+x)}
           dg.fun = function(x){1 / (1 + x)}
           ddg.fun = function(x){- 1 / (1 + x) ^ 2}
           dddg.fun = function(x){2 / (1 + x) ^ 3}
         })
  list(g.fun = g.fun, dg.fun = dg.fun, ddg.fun = ddg.fun, 
       dddg.fun = dddg.fun)
}

OtherFuns = function(copula.index){
  switch(
    as.character(copula.index),
    "3" = { # Clayton
      C.exp = expression((u1^(-para) + u2^(-para) -1)^(-1/para))
    },
    "4" = { # Gumbel
      C.exp = expression(exp(-((-log(u1))^(para) + (-log(u2))^(para))^(1/para)))
    },
    "5" = { # Frank
      C.exp = 
        expression((-1/para) * log(1 + ((exp(-para*u1)-1)*(exp(-para*u2)-1))/(exp(-para)-1)))
    }
  )
  C.exp = expression((u1^(-alpha) + u2^(-alpha) -1)^(-1/alpha))
  C.fun = function(u1, u2, alpha){}; body(C.fun) = C.exp
  dC.alpha = function(u1, u2, alpha){}; body(dC.alpha) = D(C.exp, "alpha")
  dC.u1 = function(u1, u2, alpha){}; body(dC.u1) = D(C.exp, "u1")
  dC.u1alpha = function(u1, u2, alpha){}; body(dC.u1alpha) = D(D(C.exp, "u1"), "alpha")
  dC.u2 = function(u1, u2, alpha){}; body(dC.u2) = D(C.exp, "u2")
  dC.u2alpha = function(u1, u2, alpha){}; body(dC.u2alpha) = D(D(C.exp, "u2"), "alpha")
  dC.u1u2 = function(u1, u2, alpha){}; body(dC.u1u2)= D(D(C.exp, "u1"), "u2")
  dC.u1u2alpha = function(u1, u2, alpha){}; body(dC.u1u2alpha) = D(D(D(C.exp, "u1"), "u2"), "alpha")
  list(C.fun = C.fun, dC.alpha = dC.alpha, 
       dC.u1 = dC.u1, dC.u1alpha = dC.u1alpha,
       dC.u2 = dC.u2, dC.u2alpha = dC.u2alpha,
       dC.u1u2 = dC.u1u2, dC.u1u2alpha = dC.u1u2alpha)
}
