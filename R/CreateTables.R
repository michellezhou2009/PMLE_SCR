library(dplyr)
library(tidyr)
load("RData/Summary.RData")

# Table 1 ----
tab = summary.all %>% 
  filter(TrueFam == "Gumbel", FitFam == "Gumbel", 
         para %in% c("betaT.z1", "betaT.z2", "copula")) %>%
  mutate(BIAS = format(round(BIAS / true, 3), nsmall = 3, digits =4), 
         ESD = format(round(ESD / true, 3), nsmall = 3, digits =4), 
         MSE = format(round(MSE / true, 3), nsmall = 3, digits =4)) %>%
  arrange(para, nsamp, tau, desc(betaT), method) %>%
  mutate(betaT = paste0("(1,", betaT, ")")) %>%
  select(para, tau, betaT, method, nsamp, BIAS, ESD, MSE) 
tab.200 = tab %>% filter(nsamp == 200) %>% select(-nsamp)
colnames(tab.200) = c("para", "tau", "betaT", "method", 
                      "BIAS_n200", "ESD_n200", "MSE_n200")
tab.400 = tab %>% filter(nsamp == 400) %>% select(-nsamp)
colnames(tab.400) = c("para", "tau", "betaT", "method", 
                      "BIAS_4200", "ESD_4200", "MSE_4200")
tab = left_join(tab.200, tab.400)
write.csv(tab, "Tables/tab1_Gumbel.csv", row.names = F)

# Table 2 ----
tab = summary.all %>% 
  filter(TrueFam == "Clayton", FitFam == "Clayton", 
         para %in% c("betaT.z1", "betaT.z2", "copula")) %>%
  mutate(BIAS = format(round(BIAS / true, 3), nsmall = 3, digits =4), 
         ESD = format(round(ESD / true, 3), nsmall = 3, digits =4), 
         MSE = format(round(MSE / true, 3), nsmall = 3, digits =4)) %>%
  arrange(para, nsamp, tau, desc(betaT), method) %>%
  mutate(betaT = paste0("(1,", betaT, ")")) %>%
  select(para, tau, betaT, method, nsamp, BIAS, ESD, MSE) 
tab.200 = tab %>% filter(nsamp == 200) %>% select(-nsamp)
colnames(tab.200) = c("para", "tau", "betaT", "method", 
                      "BIAS_n200", "ESD_n200", "MSE_n200")
tab.400 = tab %>% filter(nsamp == 400) %>% select(-nsamp)
colnames(tab.400) = c("para", "tau", "betaT", "method", 
                      "BIAS_4200", "ESD_4200", "MSE_4200")
tab = left_join(tab.200, tab.400)
write.csv(tab, "Tables/tab2_Clayton.csv", row.names = F)

# Table 3 ----
tab = summary.all %>% 
  filter(TrueFam == FitFam, para %in% c("betaT.z1", "betaT.z2", "copula")) %>%
  mutate(ratio = format(round(ESD / ASE, 3), nsmall = 3, digits = 4),
         TrueFam = factor(TrueFam, levels = c("Gumbel", "Clayton"))) %>%
  select(betaT, TrueFam, para, tau, nsamp, method, ratio) %>%
  spread(key = "method", value = "ratio") %>%
  arrange(betaT, TrueFam, para, tau, nsamp) %>% 
  mutate(betaT = paste0("(1,", betaT, ")"))
tab.MLE = tab %>% filter(betaT == "(1,0.5)") %>% 
  select(TrueFam, tau, nsamp, betaT, para, MLE) %>%
  spread(key = "para", value = "MLE")
colnames(tab.MLE) = c("TrueFam", "tau", "nsamp", "betaT", "MLE_betaT1",
                   "MLE_betaT2", "MLE_alpha")
tab.PMLE = tab %>% filter(betaT == "(1,0.5)") %>% 
  select(TrueFam, tau, nsamp, betaT, para, PMLE) %>%
  spread(key = "para", value = "PMLE") %>%
  select(- TrueFam, - tau, - nsamp, - betaT)
colnames(tab.PMLE) = c("PMLE_betaT1",
                   "PMLE_betaT2", "PMLE_alpha")
tab1 = cbind(tab.MLE, tab.PMLE) %>%
  select(TrueFam, tau, nsamp, betaT, MLE_betaT1, PMLE_betaT1,
         MLE_betaT2, PMLE_betaT2, MLE_alpha, PMLE_alpha) 
tab.MLE = tab %>% filter(betaT == "(1,1)") %>% 
  select(TrueFam, tau, nsamp, betaT, para, MLE) %>%
  spread(key = "para", value = "MLE")
colnames(tab.MLE) = c("TrueFam", "tau", "nsamp", "betaT", "MLE_betaT1",
                      "MLE_betaT2", "MLE_alpha")
tab.PMLE = tab %>% filter(betaT == "(1,1)") %>% 
  select(TrueFam, tau, nsamp, betaT, para, PMLE) %>%
  spread(key = "para", value = "PMLE") %>%
  select(- TrueFam, - tau, - nsamp, - betaT)
colnames(tab.PMLE) = c("PMLE_betaT1",
                       "PMLE_betaT2", "PMLE_alpha")
tab2 = cbind(tab.MLE, tab.PMLE) %>%
  select(betaT, MLE_betaT1, PMLE_betaT1,
         MLE_betaT2, PMLE_betaT2, MLE_alpha, PMLE_alpha) 
tab = cbind(tab1, tab2)
write.csv(tab, "Tables/tab3_EDS_ASE_ratio.csv", row.names = F)

# Table 4 ----
tab = summary.all %>%
  filter(
    nsamp == 400, tau == 0.6, betaT == 0.5, 
    para %in% c("betaD.z1", "betaD.z2", "betaT.z1", "betaT.z2", "Tau")) %>%
  select(TrueFam, FitFam, para, method, BIAS, ESD, MSE) %>%
  mutate(
    setting = paste0(substr(TrueFam, 1, 1), substr(FitFam, 1, 1)),
    BIAS = format(round(BIAS, 3), digits = 4, nsmall = 3),
    ESD = format(round(ESD, 3), digits = 4, nsmall = 3),
    MSE = format(round(MSE, 3), digits = 4, nsmall = 3)
  ) 
tab.gumbel = tab %>% filter(TrueFam == "Gumbel") %>%
  select(TrueFam, para, method, setting, BIAS, ESD, MSE) %>%
  gather(BIAS, ESD, MSE, key = "type", value = "value") %>%
  arrange(para, method, type, setting) %>%
  mutate(type = paste0(type, "_", setting)) %>% select(- setting) %>%
  spread(key = "type", value = "value") %>%
  select(TrueFam, para, method, BIAS_GC, BIAS_GG, ESD_GC, ESD_GG, MSE_GC, MSE_GG)
tab.clayton = tab %>% filter(TrueFam == "Clayton") %>%
  select(TrueFam, para, method, setting, BIAS, ESD, MSE) %>%
  gather(BIAS, ESD, MSE, key = "type", value = "value") %>%
  arrange(para, method, type, setting) %>%
  mutate(type = paste0(type, "_", setting)) %>% select(- setting) %>%
  spread(key = "type", value = "value") %>%
  select(TrueFam, para, method, BIAS_CG, BIAS_CC, ESD_CG, ESD_CC, MSE_CG, MSE_CC)
tab = cbind(tab.gumbel, tab.clayton)
write.csv(tab, "Tables/tab4_StudyII.csv", row.names = F)

# Table 5 ----
zhu_res = data.frame(
  alpha = c(0.2, 0.2, 0.4, 0.4, 0.6, 0.6),
  para = "copula",
  method = "Zhu et al. (2021)",
  betaT = c(0.5, 1),
  BIAS = c(0.068, 0.078, 0.051, 0.058, 0.032, 0.034),
  ESD = c(0.048, 0.056, 0.060, 0.064, 0.061, 0.059),
  tau = c(0.8, 0.8, 0.6, 0.6, 0.4, 0.4),
  nsamp = 200) %>% 
  bind_rows(data.frame(
    alpha = c(0.2, 0.2, 0.4, 0.4, 0.6, 0.6),
    para = "copula",
    method = "Zhu et al. (2021)",
    betaT = c(0.5, 1),
    BIAS = c(0.055, 0.057, 0.039, 0.041, 0.029, 0.028),
    ESD = c(0.042, 0.042, 0.056, 0.057, 0.059, 0.058),
    tau = c(0.8, 0.8, 0.6, 0.6, 0.4, 0.4),
    nsamp = 400)) %>%
  mutate(BIAS = BIAS / alpha,
         ESD = ESD / alpha,
         MSE= sqrt(BIAS ^ 2 + ESD ^ 2)) %>%
  select(tau, betaT, method, nsamp, BIAS, ESD, MSE) %>%
  mutate(BIAS = format(round(BIAS, 3), nsmall = 3, digits = 4), 
         ESD = format(round(ESD, 3), nsmall = 3, digits = 4),
         MSE = format(round(MSE, 3), nsmall = 3, digits = 4))

myres = summary.all %>% 
  filter(TrueFam == "Gumbel", FitFam == "Gumbel",
         para == "copula", method == "PMLE") %>%
  mutate(BIAS = format(round(BIAS / true, 3), nsmall = 3, digits =4), 
         ESD = format(round(ESD / true, 3), nsmall = 3, digits =4), 
         MSE = format(round(MSE / true, 3), nsmall = 3, digits =4)) %>%
  select(tau, betaT, method, nsamp, BIAS, ESD, MSE)

tab = rbind(zhu_res, myres) 
tab.200 = tab %>% filter(nsamp == 200) %>% select(- nsamp) 
colnames(tab.200) = c("tau", "betaT", "method", 
                      "BIAS_n200", "ESD_n200", "MSE_n200")
tab.400 = tab %>% filter(nsamp == 400) %>% select(- nsamp) 
colnames(tab.400) = c("tau", "betaT", "method", 
                      "BIAS_n400", "ESD_n400", "MSE_n400")
tab = left_join(tab.200, tab.400) %>% 
  arrange(tau, desc(betaT), method) %>%
  mutate(betaT = paste0("(1,", betaT, ")"))
write.csv(tab, "Tables/tab5_vs_Zhu2021.csv", row.names = F)



