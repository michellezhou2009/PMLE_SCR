rm(list = ls())
library(dplyr)
load("RData/SimRes.RData")
summary.all = res.all %>% 
  group_by(nsamp, tau, betaT, TrueFam, FitFam, method, para) %>%
  summarise(AvgEst = mean(est, na.rm = T),
            ESD = sd(est, na.rm = T), 
            ASE = mean(se, na.rm = T),
            true = unique(true)) %>% 
  ungroup() %>%
  mutate(BIAS = abs(AvgEst - true), 
         MSE = sqrt(BIAS ^ 2 + ESD ^ 2))

save(summary.all,  
     file = "RData/Summary.RData")
