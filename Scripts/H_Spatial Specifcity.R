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

check_model <- gam(
  Exponent ~ 
    Cortisol_Baseline * CD_RISC_Total * ROI + epoch + Age + Gender + Baseline_Min_from_Midnight + Residual_Delta + s(Subject, bs = "re"),  
  data = df_combined,
  family = scat(),
  method = "REML"
)
summary(check_model)
anova(check_model)
gam.check(check_model)

print(check_model$df.residual)
skewness(residuals(check_model)) 
kurtosis(residuals(check_model))
beep(1)

sim_res <- simulateResiduals(check_model, n = 10000)
plot(sim_res)

mod_mean <- mean(df_combined$CD_RISC_Total, na.rm = TRUE)
mod_sd   <- sd(df_combined$CD_RISC_Total, na.rm = TRUE)
mod_vals <- c(mod_mean - mod_sd, mod_mean + mod_sd)
slopes <- emtrends(check_model, ~ CD_RISC_Total | ROI, 
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

df_topo <- res_table_final %>%
  mutate(
    quantity = estimate,
    electrode = ROI
  )
p_coupling <- topoplot(df_topo,
                       quantity = "quantity",
                       chan_marker  = "name",
                       method = "Biharmonic",
                       head = TRUE,
                       palette = "YlOrRd",
                       fill_title="")
p_coupling <- p_coupling +
  scale_fill_gradient2(
    low = "red",
    mid = "white",
    high = "blue",
    midpoint = 0,
    limits = c(-0.2,0),
    oob = scales::squish,
    name = "β"
  )
p_coupling <- p_coupling + labs(title = "Exponent") +
  theme(plot.title = element_text(hjust = 0.5, size = 15, colour="black"))
p_coupling


check_model <- gam(
  Residual_Delta ~ 
    Cortisol_Baseline * CD_RISC_Total * ROI + epoch + Age + Gender + Baseline_Min_from_Midnight + Exponent + s(Subject, bs = "re"),  
  data = df_combined,
  family = scat(),
  method = "REML"
)
summary(check_model)
anova(check_model)
gam.check(check_model)

print(check_model$df.residual)
skewness(residuals(check_model)) 
kurtosis(residuals(check_model))
beep(1)

sim_res <- simulateResiduals(check_model, n = 10000)
plot(sim_res)

mod_mean <- mean(df_combined$CD_RISC_Total, na.rm = TRUE)
mod_sd   <- sd(df_combined$CD_RISC_Total, na.rm = TRUE)
mod_vals <- c(mod_mean - mod_sd, mod_mean + mod_sd)
slopes <- emtrends(check_model, ~ CD_RISC_Total | ROI, 
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

df_topo <- res_table_final %>%
  mutate(
    quantity = estimate,
    electrode = ROI
  )
p_coupling2 <- topoplot(df_topo,
                       quantity = "quantity",
                       chan_marker  = "name",
                       method = "Biharmonic",
                       head = TRUE,
                       palette = "YlOrRd",
                       fill_title="")
p_coupling2 <- p_coupling2 +
  scale_fill_gradient2(
    low = "red",
    mid = "white",
    high = "blue",
    midpoint = 0,
    limits = c(0,0.2),
    oob = scales::squish,
    name = "β"
  )
p_coupling2 <- p_coupling2 + labs(title = "Residual Delta") +
  theme(plot.title = element_text(hjust = 0.5, size = 15, colour="black"))
p_coupling2 

combined <- (p_coupling + p_coupling2) +
  plot_layout(widths = c(1, 1))

combined
