rm(list = ls(all.names = TRUE))

library(tidyverse)
library(mgcv)
library(emmeans)
library(eegUtils)
library(patchwork)
library(e1071)
library(DHARMa)
library(ggeffects)
library(beepr)

df_combined <- readRDS("C:/df_restingstate.rds")
roi_order <- c("O2", "O1", "P8", "P4", "Pz", "P3", "P7", "T4", "T3", "C4", "Cz", "C3", "F8", "F4", "Fz", "F3", "F7", "Fp2", "Fp1")

threeway_model <- gam(
  Residual_Alpha ~ 
    Cortisol_Baseline * CD_RISC_Total * ROI + epoch + Age + Gender + Baseline_Min_from_Midnight + s(Subject, bs = "re"),  
  data = df_combined,
  family = scat(),
  method = "REML"
)
summary(threeway_model)
anova(threeway_model)
gam.check(threeway_model)

print(threeway_model$df.residual)
skewness(residuals(threeway_model)) 
kurtosis(residuals(threeway_model))
beep(1)

sim_res <- simulateResiduals(threeway_model, n = 10000)
plot(sim_res)

mod_mean <- mean(df_combined$CD_RISC_Total, na.rm = TRUE)
mod_sd   <- sd(df_combined$CD_RISC_Total, na.rm = TRUE)
mod_vals <- c(mod_mean - mod_sd, mod_mean + mod_sd)
slopes <- emtrends(threeway_model, ~ CD_RISC_Total | ROI, 
                   var = "Cortisol_Baseline", 
                   at = list(CD_RISC_Total = mod_vals))
diffs <- contrast(slopes, method = "pairwise", adjust = "none")
res_table <- as.data.frame(summary(diffs, infer = c(TRUE, TRUE), adjust = "none"),level = 0.95)
res_table$p_holm <- p.adjust(res_table$p.value, method = "holm")
res_table$significance <- ifelse(res_table$p_holm < 0.05, "*", "")
res_table_final <- res_table[order(res_table$p_holm), ]
print(res_table_final[, c("ROI", "estimate", "t.ratio", "p.value", "p_holm", "significance")])
res_filtered <- res_table_final[res_table_final$p_holm < 0.05, ]
print(res_filtered[, c("ROI", "estimate", "t.ratio", "p.value", "p_holm")])