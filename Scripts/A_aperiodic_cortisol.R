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

twoway_model <- gam(
  Exponent ~ 
    Cortisol_Baseline * ROI + epoch + Age + Gender + Baseline_Min_from_Midnight + s(Subject, bs = "re"),  
  data = df_combined,
  family = scat(),
  method = "REML"
)
summary(twoway_model)
anova(twoway_model)
gam.check(twoway_model)

print(twoway_model$df.residual)
skewness(residuals(twoway_model)) 
kurtosis(residuals(twoway_model))
beep(1)

sim_res <- simulateResiduals(twoway_model, n = 10000)
plot(sim_res)

slopes <- emtrends(twoway_model, ~ ROI, var = "Cortisol_Baseline")
diffs <- contrast(slopes, method = "pairwise", adjust = "none")
res_table <- as.data.frame(summary(diffs, infer = c(TRUE, TRUE), adjust = "none"),level = 0.95)
res_table$p_holm <- p.adjust(res_table$p.value, method = "holm")
res_table$significance <- ifelse(res_table$p_holm < 0.05, "*", "")
res_table_final <- res_table[order(res_table$p_holm), ]
print(res_table_final[, c("contrast", "estimate", "t.ratio", "p.value", "p_holm", "significance")])
res_filtered <- res_table_final[res_table_final$p_holm < 0.05, ]
print(res_filtered[, c("contrast", "estimate", "t.ratio", "p.value", "p_holm")])

emm_df <- as.data.frame(slopes)
df_topo <- data.frame(
  electrode = emm_df$ROI,
  quantity = emm_df[["Cortisol_Baseline.trend"]]
)

p_coupling <- topoplot(df_topo,
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
    limits = c(-0.08,0),
    oob = scales::squish,
    name = "b"
  )
p_coupling <- p_coupling + 
  theme(plot.title = element_text(hjust = 0.5, size = 10, colour="dark red"))
p_coupling

res_split <- res_table_final %>%
  separate(
    contrast,
    into = c("ROI_A", "ROI_B"),
    sep = " - "
  )
res_mirror <- res_split %>%
  mutate(
    ROI_A_old = ROI_A,
    ROI_B_old = ROI_B
  ) %>%
  transmute(
    ROI_A = ROI_B_old,
    ROI_B = ROI_A_old,
    estimate = -estimate,
    SE = SE,
    df = df,
    lower.CL = -upper.CL,
    upper.CL = -lower.CL,
    t.ratio = -t.ratio,
    p.value = p.value,
    p_holm = p_holm
  )
res_full <- bind_rows(res_split, res_mirror)
res_full$signif <- case_when(
  res_full$p_holm < 0.001 ~ "***",
  res_full$p_holm < 0.01  ~ "**",
  res_full$p_holm < 0.05  ~ "*",
  TRUE ~ ""
)
res_final <- res_full %>%
  select(ROI_A, ROI_B, estimate, signif)
mat_est <- res_full %>%
  select(ROI_A, ROI_B, estimate) %>%
  pivot_wider(
    names_from = ROI_B,
    values_from = estimate
  ) %>%
  column_to_rownames("ROI_A") %>%
  as.matrix()
mat_est <- mat_est[rev(roi_order), rev(roi_order)]

mat_est_sig <- res_full %>%
  select(ROI_A, ROI_B, signif) %>%
  pivot_wider(
    names_from = ROI_B,
    values_from = signif
  ) %>%
  column_to_rownames("ROI_A") %>%
  as.matrix()
mat_est_sig <- mat_est_sig[rev(roi_order), rev(roi_order)]
est_df <- as.data.frame.table(mat_est) %>% 
  rename(ROI_A = Var1, ROI_B = Var2, estimate = Freq)
sig_df <- as.data.frame.table(mat_est_sig) %>% 
  rename(ROI_A = Var1, ROI_B = Var2, signif = Freq)

plot_df <- merge(est_df, sig_df, by = c("ROI_A", "ROI_B"))
plot_df <- plot_df %>%
  mutate(
    row = match(ROI_A, rownames(mat_est)),
    col = match(ROI_B, colnames(mat_est))
  ) %>%
  filter(row >= col)
matrix_plot <-ggplot(plot_df, aes(x = ROI_B, y = ROI_A, fill = estimate)) +
  geom_tile() +
  geom_text(aes(label = signif), size = 2) +
  coord_fixed() +
  scale_y_discrete(limits = rev(rownames(mat_est))) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5)
  ) +
  labs(
    x = NULL,
    y = NULL
  ) +
  scale_fill_gradient2(
    low = "red",
    mid = "white",
    high = "blue",
    midpoint = 0,
    limits = c(-0.04,0.04),
    oob = scales::squish,
    name = "Constrast Estimate"
  ) +
  theme(
    legend.title.position = "right",
    legend.title = element_text(
      angle = 270,
      vjust = 0.5,
      hjust = 0.5
    )
  )
matrix_plot

combined <- p_coupling + matrix_plot + 
  plot_layout(widths = c(1, 1))
combined
