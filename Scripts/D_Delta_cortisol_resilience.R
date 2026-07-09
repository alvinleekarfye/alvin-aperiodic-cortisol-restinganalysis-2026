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
  Residual_Delta ~ 
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

keep_rois <-res_filtered$ROI
pred_data <- predict_response(
  threeway_model,
  terms = c(
    "Cortisol_Baseline",
    paste0("CD_RISC_Total [", paste(mod_vals, collapse = ","), "]"),
    "ROI"
  ),  ci = 0.95)
pred_data$facet <- factor(pred_data$facet, levels = rev(roi_order))
pred_data <- pred_data[pred_data$facet %in% keep_rois, ]

pred_data$group <- factor(
  pred_data$group,
  labels = c("Low Resilience", "High Resilience")
)
pred_data$group <- factor(pred_data$group, levels = rev(levels(pred_data$group)))
p1 <- plot(pred_data) +
  facet_wrap(~facet, ncol = 4, nrow = 2) +
  labs(title = NULL, x = "Cortisol Concentration", y = "Residual Delta") +
  theme_minimal()+
  theme(
    legend.title = element_blank(),
    legend.position = "top",      
    legend.direction = "horizontal",
    legend.justification = "center"   
  ) +
  scale_color_manual(values = c(
    "Low Resilience" = "dark red", 
    "High Resilience" = "dark green"
  )) +
  scale_fill_manual(values = c(
    "Low Resilience" = "dark red", 
    "High Resilience" = "dark green"
  ))
p1

posthoc_df <- res_table_final %>%
  rename(estimate = estimate) %>%   
  mutate(
    signif = case_when(
      p_holm < 0.001 ~ "***",
      p_holm < 0.01  ~ "**",
      p_holm < 0.05  ~ "*",
      TRUE ~ " "
    )
  )
posthoc_df <- posthoc_df %>%
  mutate(ROI = factor(ROI, levels = roi_order))

p_coupling_forest <- ggplot(posthoc_df, aes(x = estimate, y = ROI)) +
  geom_point(size = 3, aes(color = signif)) +
  geom_errorbarh(aes(xmin = lower.CL, xmax = upper.CL), height = 0.5) +
  geom_text(aes(x = upper.CL, label = signif), hjust = -0.3, size = 4) +
  scale_color_manual(values = c("*" = "green", "**" = "green", "***" = "green")) +
  labs(x = "b Contrast (Low vs High Resilience)", y = NULL, color = "Significance") +
  theme_minimal() +
  theme(legend.position = "none")
p_coupling_forest <- p_coupling_forest + coord_cartesian(xlim = c(-0.5, 0.5))
p_coupling_forest

emm_df <- as.data.frame(slopes)
df_topo <- data.frame(
  electrode = emm_df$ROI,
  CD_RISC_Total = emm_df$CD_RISC_Total,
  quantity = emm_df[["Cortisol_Baseline.trend"]]
)
low_df  <- subset(df_topo, df_topo$CD_RISC_Total == mod_vals[1])
high_df <- subset(df_topo, df_topo$CD_RISC_Total == mod_vals[2])

p_coupling <- topoplot(low_df,
                       quantity = "quantity",
                       chan_marker  = "name",
                       method = "Biharmonic",
                       head = TRUE,
                       palette = "YlOrRd",
                       fill_title="b")
p_coupling <- p_coupling +
  scale_fill_gradient2(
    low = "red",
    mid = "white",
    high = "blue",
    midpoint = 0,
    limits = c(-0.18,0.06),
    oob = scales::squish,
    name = "b"
  )
p_coupling <- p_coupling + labs(title = "Low Resilience") +
  theme(plot.title = element_text(hjust = 0.5, size = 10, colour="dark red")) +
  guides(fill = "none")
p_coupling

p_coupling2 <- topoplot(high_df,
                        quantity = "quantity",
                        chan_marker  = "name",
                        method = "Biharmonic",
                        head = TRUE,
                        palette = "YlOrRd",
                        fill_title="b")
p_coupling2 <- p_coupling2 +
  scale_fill_gradient2(
    low = "red",
    mid = "white",
    high = "blue",
    midpoint = 0,
    limits = c(-0.18,0.06),
    oob = scales::squish,
    name = "b"
  )
p_coupling2 <- p_coupling2 + labs(title = "High Resilience") +
  theme(plot.title = element_text(hjust = 0.5, size = 10, colour="dark green"))
p_coupling2

combined <- (p_coupling + p_coupling_forest + p_coupling2) /
  p1 +
  plot_layout(widths = c(1, 1, 1))

combined