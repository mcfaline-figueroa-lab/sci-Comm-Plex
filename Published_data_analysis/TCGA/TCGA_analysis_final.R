# Clear environment

rm(list=ls())
gc()  # Garbage collection to free memory

library(tidyverse)
library(survival)
library(lubridate)
library(ggsurvfit)
library(gtsummary)
library(tidycmprsk)
if (requireNamespace("gghalves", quietly = TRUE)) {
  library(gghalves)
}
#library(condsurv)

.plot_font_family <- "sans"

#conda activate tcga_analysis_monocle

# download data from https://www.cbioportal.org/ and select TCGA PanCancer Atlas Studies
# query by gene to download alteration and expression data
# explore selected studies to download corresponding clinical information
setwd("/home/user/Documents/Kinase_project/TCGA_data/")

dir.create("figures")

clin_data <- read_tsv(file = "data/combined_study_clinical_data.tsv",
                      col_names = TRUE)

expression_data <- read_tsv(file = "data/mRNA Expression, RSEM (Batch normalized from Illumina HiSeq_RNASeqV2)_PDGFRA.txt",
                            col_names = T)

colnames(expression_data) <- paste0(colnames(expression_data), "_", "exprs")

cn_data <-  read_tsv(file = "data/copy_number_alteration.txt",
                     col_names = T)

colnames(cn_data) <- paste0(colnames(cn_data), "_", "cn")

mut_data <- read_tsv(file = "data/mutations.txt",
                     col_names = T)

colnames(mut_data) <- paste0(colnames(mut_data), "_", "mut")

all_data <- left_join(clin_data,
                      expression_data,
                      by = c("Sample ID" = "SAMPLE_ID_exprs")) %>%
  left_join(cn_data,
            by = c("Sample ID" = "SAMPLE_ID_cn")) %>%
  left_join(mut_data,
            by = c("Sample ID" = "SAMPLE_ID_mut")) %>%
  filter(!is.na(PDGFRA_cn)) %>% # filter samples with no copy number data
  filter(!is.na(PDGFRA_exprs)) %>% # filter samples with no expression data
  mutate(log2_PDGFRA_exprs = log2(PDGFRA_exprs + 1)) %>%
  mutate(PDGFRA_cn_alteration = case_when(
    PDGFRA_cn > 0 ~ "Amplified",
    PDGFRA_cn < 0 ~ "Deleted",
    TRUE ~ "WT"
  )) %>%
  mutate(PDGFRA_cn_alteration = factor(PDGFRA_cn_alteration, levels = c("WT", "Amplified"))) %>%
  filter(PDGFRA_cn_alteration %in% c("WT", "Amplified"))  # Only keep WT and Amplified

# Expression of T-cell dose-response gene score
dose_response_expression_data <- read_tsv(file = "data/mRNA Expression, RSEM (Batch normalized from Illumina HiSeq_RNASeqV2)_Tcell_dose_response.txt",
                            col_names = T) %>%
  pivot_longer(cols = 3:last_col(),
               names_to = "gene",
               values_to = "expression") %>%
  drop_na() %>%
  group_by(STUDY_ID, SAMPLE_ID) %>%
  summarise(count = n(),
            agg_expression = mean(log2(expression + 1))) %>%
  rename(log2_agg_expression = agg_expression)
    
all_data <- left_join(all_data,
                      dose_response_expression_data,
                      by = c("Sample ID" = "SAMPLE_ID",
                             "Study ID" = "STUDY_ID")) %>%
  filter(!is.na(log2_agg_expression))

cat("\n===== SAMPLE COUNTS AT EACH STEP =====\n")
cat("Total pancancer samples after all filters:", nrow(all_data), "\n")

gbm_before_filter <- all_data %>% filter(`Cancer Type` == "Glioblastoma")
cat("GBM samples before subtype filter:", nrow(gbm_before_filter), "\n")

median_exprs <- all_data %>%
  group_by(PDGFRA_cn_alteration) %>%
  summarise(count = n(),
            mean_PDGFRA_exprs = mean(log2_PDGFRA_exprs),
            median_PDGFRA_exprs = median(log2_PDGFRA_exprs),
            mean_log2_agg = mean(log2_agg_expression),
            median_log2_agg = median(log2_agg_expression))

split_violin_data <- all_data %>%
  mutate(x_axis = 1) %>%
  mutate(cn_split = PDGFRA_cn_alteration) %>%
  mutate(cn_split = factor(cn_split, levels = c("WT", "Amplified")))

# Wilcoxon test for pancancer
wt_data <- split_violin_data %>% filter(PDGFRA_cn_alteration == "WT") %>% pull(log2_agg_expression)
ampl_data <- split_violin_data %>% filter(PDGFRA_cn_alteration == "Amplified") %>% pull(log2_agg_expression)
wilcox_test_pancancer <- wilcox.test(wt_data, ampl_data)
p_val_pancancer <- wilcox_test_pancancer$p.value

# Add sample counts for x-axis labels
sample_counts_pancancer <- split_violin_data %>%
  group_by(PDGFRA_cn_alteration) %>%
  summarise(n = n(), .groups = 'drop')

split_violin_data <- split_violin_data %>%
  left_join(sample_counts_pancancer, by = "PDGFRA_cn_alteration") %>%
  mutate(label = paste0(PDGFRA_cn_alteration, "\n(n = ", n, ")"))

ggplot(data = split_violin_data,
       aes(x = factor(label, levels = unique(label[order(PDGFRA_cn_alteration)])), 
           y = log2_agg_expression,
           fill = PDGFRA_cn_alteration)) +
  geom_boxplot(width = 0.4, alpha = 0.9, outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.5, size = 1.5) +
  scale_fill_manual(values = c("WT" = "#edf2f4", "Amplified" = "#8d99ae"),
                    guide = "none") +
  xlab("PDGFRA Status") +
  ylab("Up-regulated Signature Score (log2)") +
  ggtitle("Pan Cancer: T-cell Dose Response Signature vs PDGFRA CN Status") +
  annotate(geom = "text",
           x = 1.5, y = Inf,
           label = paste0("p = ", format(p_val_pancancer, digits = 3)),
           vjust = 1.5, hjust = 0.5,
           size = 3,
         family = .plot_font_family) +
  theme_minimal() +
    theme(axis.text = element_text(size = 10, family = .plot_font_family),
      axis.title = element_text(size = 11, family = .plot_font_family),
    plot.title = element_text(hjust = 0.5, face = "bold", family = .plot_font_family),
        plot.margin = margin(t = 10, r = 10, b = 10, l = 10),
        panel.grid = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
        axis.ticks = element_line(colour = "black", linewidth = 0.3))
ggsave("figures/T-cell_dose_response_genes_vs_PDGFRA_alteration_split_pancancer.png",
       dpi = 300, width = 4, height = 4)

# Pancancer: PDGFRA expression
# Wilcoxon test for PDGFRA expression in pancancer
wt_pdgfra_pan <- all_data %>% filter(PDGFRA_cn_alteration == "WT") %>% pull(log2_PDGFRA_exprs)
ampl_pdgfra_pan <- all_data %>% filter(PDGFRA_cn_alteration == "Amplified") %>% pull(log2_PDGFRA_exprs)
wilcox_test_pdgfra_pan <- wilcox.test(wt_pdgfra_pan, ampl_pdgfra_pan)
p_val_pdgfra_pan <- wilcox_test_pdgfra_pan$p.value

# Add sample counts
sample_counts_pdgfra_pan <- all_data %>%
  group_by(PDGFRA_cn_alteration) %>%
  summarise(n = n(), .groups = 'drop')

all_data_pan <- all_data %>%
  left_join(sample_counts_pdgfra_pan, by = "PDGFRA_cn_alteration") %>%
  mutate(label = paste0(PDGFRA_cn_alteration, "\n(n = ", n, ")"))

# Calculate median values for labels
median_values_pan <- all_data_pan %>%
  group_by(PDGFRA_cn_alteration) %>%
  summarise(median = median(log2_PDGFRA_exprs), .groups = 'drop')

p_pdgfra_pan <- ggplot(all_data_pan, aes(x = factor(label, levels = unique(label[order(PDGFRA_cn_alteration)])), 
                                          y = log2_PDGFRA_exprs, 
                                          fill = PDGFRA_cn_alteration)) +
  geom_boxplot(width = 0.4, alpha = 0.9, outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.5, size = 1.5) +
  geom_text(data = median_values_pan, 
            aes(x = c(1, 2), y = median, label = sprintf("%.2f", median)),
            vjust = -0.5, size = 3, fontface = "bold") +
  scale_fill_manual(values = c("WT" = "#edf2f4", "Amplified" = "#8d99ae"),
                    guide = "none") +
  xlab("PDGFRA Status") +
  ylab(expression(PDGFRA~Expression~(log[2](TPM + 1)))) +
  ggtitle("Pan Cancer: PDGFRA Expression by Copy Number Status") +
  annotate(geom = "text",
           x = 1.5, y = Inf,
           label = paste0("p = ", format(p_val_pdgfra_pan, digits = 3)),
           vjust = 1.5, hjust = 0.5,
           size = 3) +
  annotate(geom = "text",
           x = -Inf, y = Inf,
           label = "Pan Cancer",
           vjust = 1.5, hjust = -0.5,
           size = 4, fontface = "bold") +
  annotate(geom = "text",
           x = Inf, y = -Inf,
           label = "TCGA",
           vjust = -0.5, hjust = 1.2,
           size = 3, style = "italic") +
  theme_minimal() +
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 11),
      plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.margin = margin(t = 10, r = 10, b = 10, l = 10),
        panel.grid = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
        axis.ticks = element_line(colour = "black", linewidth = 0.3))

ggsave("figures/PDGFRA_expression_pancancer.png", p_pdgfra_pan, dpi = 300, width = 6, height = 6)
print(p_pdgfra_pan)

gbm_data <- all_data %>%
  filter(`Cancer Type` == "Glioblastoma") %>%
  filter(Subtype == "GBM_IDHwt")

cat("GBM samples in all_data (after T-cell signature filter):", nrow(gbm_data), "\n")
cat("  - WT samples:", nrow(gbm_data %>% filter(PDGFRA_cn_alteration == "WT")), "\n")
cat("  - Amplified samples:", nrow(gbm_data %>% filter(PDGFRA_cn_alteration == "Amplified")), "\n\n")

# # Check GBM Subtype distribution
# cat("\n===== GBM SUBTYPE DISTRIBUTION (Pre-filter) =====\n")
# subtype_counts <- gbm_data %>%
#   group_by(Subtype, .groups = 'drop') %>%
#   summarise(count = n(), .groups = 'drop') %>%
#   arrange(desc(count))
# print(subtype_counts)

# # # Filter to only GBM_IDHwt
# # gbm_data <- gbm_data %>%
# #   filter(Subtype == "GBM_IDHwt")

# cat("After GBM_IDHwt filter: n =", nrow(gbm_data), "\n\n")

median_exprs <- gbm_data %>%
  group_by(PDGFRA_cn_alteration) %>%
  summarise(count = n(),
            mean_PDGFRA_exprs = mean(log2_PDGFRA_exprs),
            median_PDGFRA_exprs = median(log2_PDGFRA_exprs),
            mean_log2_agg = mean(log2_agg_expression),
            median_log2_agg = median(log2_agg_expression))


split_violin_data <- gbm_data %>%
  mutate(x_axis = 1) %>%
  mutate(cn_split = PDGFRA_cn_alteration) %>%
  mutate(cn_split = factor(cn_split, levels = c("WT", "Amplified")))

# Wilcoxon test for GBM
wt_data_gbm <- split_violin_data %>% filter(PDGFRA_cn_alteration == "WT") %>% pull(log2_agg_expression)
ampl_data_gbm <- split_violin_data %>% filter(PDGFRA_cn_alteration == "Amplified") %>% pull(log2_agg_expression)
wilcox_test_gbm <- wilcox.test(wt_data_gbm, ampl_data_gbm)
p_val_gbm <- wilcox_test_gbm$p.value

# Add sample counts for x-axis labels
sample_counts_gbm <- split_violin_data %>%
  group_by(PDGFRA_cn_alteration) %>%
  summarise(n = n(), .groups = 'drop')

split_violin_data <- split_violin_data %>%
  left_join(sample_counts_gbm, by = "PDGFRA_cn_alteration") %>%
  mutate(label = paste0(PDGFRA_cn_alteration, "\n(n = ", n, ")"))

ggplot(data = split_violin_data,
       aes(x = factor(label, levels = unique(label[order(PDGFRA_cn_alteration)])), 
           y = log2_agg_expression,
           fill = PDGFRA_cn_alteration)) +
  geom_boxplot(width = 0.4, alpha = 0.9, outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.5, size = 1.5) +
  scale_fill_manual(values = c("WT" = "#edf2f4", "Amplified" = "#8d99ae"),
                    guide = "none") +
  xlab("PDGFRA Status") +
  ylab("Up-regulated Signature Score (log2)") +
  ggtitle("GBM (IDHwt): T-cell Dose Response Signature vs PDGFRA CN Status") +
  annotate(geom = "text",
           x = 1.5, y = Inf,
           label = paste0("p = ", format(p_val_gbm, digits = 3)),
           vjust = 1.5, hjust = 0.5,
           size = 3,
         family = .plot_font_family) +
  theme_minimal() +
    theme(axis.text = element_text(size = 10, family = .plot_font_family),
      axis.title = element_text(size = 11, family = .plot_font_family),
    plot.title = element_text(hjust = 0.5, face = "bold", family = .plot_font_family),
        plot.margin = margin(t = 10, r = 10, b = 10, l = 10),
        panel.grid = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
        axis.ticks = element_line(colour = "black", linewidth = 0.3))
ggsave("figures/T-cell_dose_response_genes_vs_PDGFRA_alteration_split_gbm.png",
       dpi = 300, width = 6, height = 6)

corr_exprs <- cor.test(x = gbm_data$log2_PDGFRA_exprs, 
                  y = gbm_data$log2_agg_expression,
                  method = "pearson")


p <- ggplot(gbm_data,
       aes(x = log2_PDGFRA_exprs,
           y = log2_agg_expression)) +
  geom_point(aes(color = PDGFRA_cn_alteration),
             size = 0.5) +
  geom_smooth(method = "lm",
              linewidth = 0.5) +
  xlab("Log2 PDGFRA Expression") +
  ylab("Log2 Agg Expression\nT-cell Dose Response Genes") +
  labs(title = "GBM (IDHwt): PDGFRA vs T-cell Dose Response Signature") +
  annotate(geom = "text",
           x = 4.75,
           y = 5.45,
           label = paste0("r = ", round(corr_exprs$estimate, 3), "\np = ", round(corr_exprs$p.value, 4)),
           size = 2) +
  scale_color_manual(values = c("WT" = "#999999", "Amplified" = "#FF6633"),
                    name = "CN Status") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        text = element_text(size = 8),
    plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "bottom",
        legend.key.height = unit(0.3, "cm"),
        legend.key.width = unit(0.3, "cm"),
        legend.margin = margin(t = -5, l = -5)) +
  guides(color = guide_legend(title = "PDGFRA CN Status",
                             byrow = T, 
                             nrow = 1))

ggsave("figures/gbm_scatter_pdgfra_vs_tcell.png", p, dpi = 300, width = 6, height = 5)

p1 <- ggExtra::ggMarginal(p, 
                          type = "density",
                          groupColour = TRUE, 
                          groupFill = TRUE, 
                          alpha = 0.25, size = 5)
p1

p1 <- ggExtra::ggMarginal(p, 
                          type = "density",
                          groupColour = TRUE, 
                          groupFill = TRUE, 
                          alpha = 0.25, size = 5)
ggsave("figures/gbm_scatter_pdgfra_vs_tcell_marginal.png", p1, dpi = 300, width = 6, height = 5)
p1


# GBM: PDGFRA expression
# Wilcoxon test for PDGFRA expression in GBM
wt_pdgfra <- gbm_data %>% filter(PDGFRA_cn_alteration == "WT") %>% pull(log2_PDGFRA_exprs)
ampl_pdgfra <- gbm_data %>% filter(PDGFRA_cn_alteration == "Amplified") %>% pull(log2_PDGFRA_exprs)
wilcox_test_pdgfra <- wilcox.test(wt_pdgfra, ampl_pdgfra)
p_val_pdgfra <- wilcox_test_pdgfra$p.value

# Add sample counts
sample_counts_pdgfra_gbm <- gbm_data %>%
  group_by(PDGFRA_cn_alteration) %>%
  summarise(n = n(), .groups = 'drop')

gbm_data_labeled <- gbm_data %>%
  left_join(sample_counts_pdgfra_gbm, by = "PDGFRA_cn_alteration") %>%
  mutate(label = paste0(PDGFRA_cn_alteration, "\n(n = ", n, ")"))

# Calculate median values for labels
median_values_gbm <- gbm_data_labeled %>%
  group_by(PDGFRA_cn_alteration) %>%
  summarise(median = median(log2_PDGFRA_exprs), .groups = 'drop')

p5 <- ggplot(gbm_data_labeled, aes(x = factor(label, levels = unique(label[order(PDGFRA_cn_alteration)])), 
                                    y = log2_PDGFRA_exprs, 
                                    fill = PDGFRA_cn_alteration)) +
  geom_boxplot(width = 0.4, alpha = 0.9, outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.5, size = 1.5) +
  geom_text(data = median_values_gbm, 
            aes(x = c(1, 2), y = median, label = sprintf("%.2f", median)),
            vjust = -0.5, size = 3, fontface = "bold") +
  scale_fill_manual(values = c("WT" = "#edf2f4", "Amplified" = "#8d99ae"),
                    guide = "none") +
  xlab("PDGFRA Status") +
  ylab(expression(PDGFRA~Expression~(log[2](TPM + 1)))) +
  ggtitle("GBM (IDHwt): PDGFRA Expression by Copy Number Status") +
  annotate(geom = "text",
           x = 1.5, y = Inf,
           label = paste0("p = ", format(p_val_pdgfra, digits = 3)),
           vjust = 1.5, hjust = 0.5,
           size = 3) +
  annotate(geom = "text",
           x = -Inf, y = Inf,
           label = "GBM",
           vjust = 1.5, hjust = -0.5,
           size = 4, fontface = "bold") +
  annotate(geom = "text",
           x = Inf, y = -Inf,
           label = "TCGA",
           vjust = -0.5, hjust = 1.2,
           size = 3, style = "italic") +
  theme_minimal() +
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 11),
      plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.margin = margin(t = 10, r = 10, b = 10, l = 10),
        panel.grid = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
        axis.ticks = element_line(colour = "black", linewidth = 0.3))

ggsave("figures/PDGFRA_expression_gbm.png", p5, dpi = 300, width = 6, height = 6)
print(p5)

# ============================================================================
# GLIOBLASTOMA + GLIOMA COMBINED ANALYSIS
# ============================================================================

glioma_data <- all_data %>%
  filter(`Cancer Type` %in% c("Glioblastoma", "Glioma"))

median_exprs_glioma <- glioma_data %>%
  group_by(PDGFRA_cn_alteration) %>%
  summarise(count = n(),
            mean_PDGFRA_exprs = mean(log2_PDGFRA_exprs),
            median_PDGFRA_exprs = median(log2_PDGFRA_exprs),
            mean_log2_agg = mean(log2_agg_expression),
            median_log2_agg = median(log2_agg_expression))

split_violin_data_glioma <- glioma_data %>%
  mutate(x_axis = 1) %>%
  mutate(cn_split = PDGFRA_cn_alteration) %>%
  mutate(cn_split = factor(cn_split, levels = c("WT", "Amplified")))

# Wilcoxon test for Glioblastoma + Glioma
wt_data_glioma <- split_violin_data_glioma %>% filter(PDGFRA_cn_alteration == "WT") %>% pull(log2_agg_expression)
ampl_data_glioma <- split_violin_data_glioma %>% filter(PDGFRA_cn_alteration == "Amplified") %>% pull(log2_agg_expression)
wilcox_test_glioma <- wilcox.test(wt_data_glioma, ampl_data_glioma)
p_val_glioma <- wilcox_test_glioma$p.value

# Add sample counts for x-axis labels
sample_counts_glioma <- split_violin_data_glioma %>%
  group_by(PDGFRA_cn_alteration) %>%
  summarise(n = n(), .groups = 'drop')

split_violin_data_glioma <- split_violin_data_glioma %>%
  left_join(sample_counts_glioma, by = "PDGFRA_cn_alteration") %>%
  mutate(label = paste0(PDGFRA_cn_alteration, "\n(n = ", n, ")"))

ggplot(data = split_violin_data_glioma,
       aes(x = factor(label, levels = unique(label[order(PDGFRA_cn_alteration)])), 
           y = log2_agg_expression,
           fill = PDGFRA_cn_alteration)) +
  geom_boxplot(width = 0.4, alpha = 0.9, outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.5, size = 1.5) +
  scale_fill_manual(values = c("WT" = "#edf2f4", "Amplified" = "#8d99ae"),
                    guide = "none") +
  xlab("PDGFRA Status") +
  ylab("Up-regulated Signature Score (log2)") +
  ggtitle("Glioblastoma + Glioma: T-cell Dose Response Signature vs PDGFRA CN Status") +
  annotate(geom = "text",
           x = 1.5, y = Inf,
           label = paste0("p = ", format(p_val_glioma, digits = 3)),
           vjust = 1.5, hjust = 0.5,
           size = 3,
         family = .plot_font_family) +
  theme_minimal() +
    theme(axis.text = element_text(size = 10, family = .plot_font_family),
      axis.title = element_text(size = 11, family = .plot_font_family),
    plot.title = element_text(hjust = 0.5, face = "bold", family = .plot_font_family),
        plot.margin = margin(t = 10, r = 10, b = 10, l = 10),
        panel.grid = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
        axis.ticks = element_line(colour = "black", linewidth = 0.3))
ggsave("figures/T-cell_dose_response_genes_vs_PDGFRA_alteration_split_glioblastoma_glioma.png",
       dpi = 300, width = 4, height = 4)

corr_exprs_glioma <- cor.test(x = glioma_data$log2_PDGFRA_exprs, 
                         y = glioma_data$log2_agg_expression,
                         method = "pearson")

p_glioma <- ggplot(glioma_data,
              aes(x = log2_PDGFRA_exprs,
                  y = log2_agg_expression)) +
  geom_point(aes(color = PDGFRA_cn_alteration),
             size = 0.5) +
  geom_smooth(method = "lm",
              linewidth = 0.5) +
  xlab("Log2 PDGFRA Expression") +
  ylab("Log2 Agg Expression\nT-cell Dose Response Genes") +
  labs(title = "Glioblastoma + Glioma: PDGFRA vs T-cell Dose Response Signature") +
  annotate(geom = "text",
           x = 4.75,
           y = 5.45,
           label = paste0("r = ", round(corr_exprs_glioma$estimate, 3), "\np = ", round(corr_exprs_glioma$p.value, 4)),
           size = 2) +
  scale_color_manual(values = c("WT" = "#999999", "Amplified" = "#FF6633"),
                    name = "CN Status") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        text = element_text(size = 8),
    plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "bottom",
        legend.key.height = unit(0.3, "cm"),
        legend.key.width = unit(0.3, "cm"),
        legend.margin = margin(t = -5, l = -5)) +
  guides(color = guide_legend(title = "PDGFRA CN Status",
                             byrow = T, 
                             nrow = 1))

ggsave("figures/glioblastoma_glioma_scatter_pdgfra_vs_tcell.png", p_glioma, dpi = 300, width = 6, height = 5)

p1_glioma <- ggExtra::ggMarginal(p_glioma, 
                                 type = "density",
                                 groupColour = TRUE, 
                                 groupFill = TRUE, 
                                 alpha = 0.25, size = 5)
ggsave("figures/glioblastoma_glioma_scatter_pdgfra_vs_tcell_marginal.png", p1_glioma, dpi = 300, width = 6, height = 5)

# Glioblastoma + Glioma: PDGFRA expression
# Wilcoxon test for PDGFRA expression
wt_pdgfra_glioma <- glioma_data %>% filter(PDGFRA_cn_alteration == "WT") %>% pull(log2_PDGFRA_exprs)
ampl_pdgfra_glioma <- glioma_data %>% filter(PDGFRA_cn_alteration == "Amplified") %>% pull(log2_PDGFRA_exprs)
wilcox_test_pdgfra_glioma <- wilcox.test(wt_pdgfra_glioma, ampl_pdgfra_glioma)
p_val_pdgfra_glioma <- wilcox_test_pdgfra_glioma$p.value

# Add sample counts
sample_counts_pdgfra_glioma <- glioma_data %>%
  group_by(PDGFRA_cn_alteration) %>%
  summarise(n = n(), .groups = 'drop')

glioma_data_labeled <- glioma_data %>%
  left_join(sample_counts_pdgfra_glioma, by = "PDGFRA_cn_alteration") %>%
  mutate(label = paste0(PDGFRA_cn_alteration, "\n(n = ", n, ")"))

# Calculate median values for labels
median_values_glioma <- glioma_data_labeled %>%
  group_by(PDGFRA_cn_alteration) %>%
  summarise(median = median(log2_PDGFRA_exprs), .groups = 'drop')

p_glioma_pdgfra <- ggplot(glioma_data_labeled, aes(x = factor(label, levels = unique(label[order(PDGFRA_cn_alteration)])), 
                                                    y = log2_PDGFRA_exprs, 
                                                    fill = PDGFRA_cn_alteration)) +
  geom_boxplot(width = 0.4, alpha = 0.9, outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.5, size = 1.5) +
  geom_text(data = median_values_glioma, 
            aes(x = c(1, 2), y = median, label = sprintf("%.2f", median)),
            vjust = -0.5, size = 3, fontface = "bold") +
  scale_fill_manual(values = c("WT" = "#edf2f4", "Amplified" = "#8d99ae"),
                    guide = "none") +
  xlab("PDGFRA Status") +
  ylab(expression(PDGFRA~Expression~(log[2](TPM + 1)))) +
  ggtitle("Glioblastoma + Glioma: PDGFRA Expression by Copy Number Status") +
  annotate(geom = "text",
           x = 1.5, y = Inf,
           label = paste0("p = ", format(p_val_pdgfra_glioma, digits = 3)),
           vjust = 1.5, hjust = 0.5,
           size = 3,
           family = .plot_font_family) +
  annotate(geom = "text",
           x = -Inf, y = Inf,
           label = "GBM + Glioma",
           vjust = 1.5, hjust = -0.5,
           size = 4, fontface = "bold") +
  annotate(geom = "text",
           x = Inf, y = -Inf,
           label = "TCGA",
           vjust = -0.5, hjust = 1.2,
           size = 3, style = "italic") +
  theme_minimal() +
  theme(axis.text = element_text(size = 10, family = .plot_font_family),
        axis.title = element_text(size = 11, family = .plot_font_family),
      plot.title = element_text(hjust = 0.5, face = "bold", family = .plot_font_family),
        plot.margin = margin(t = 10, r = 10, b = 10, l = 10),
        panel.grid = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
        axis.ticks = element_line(colour = "black", linewidth = 0.3))

ggsave("figures/PDGFRA_expression_glioblastoma_glioma.png", p_glioma_pdgfra, dpi = 300, width = 6, height = 6)
print(p_glioma_pdgfra)

# ============================================================================
# INDIVIDUAL GENE ANALYSIS: SOD2, CD274, CSF3, IDO1
# ============================================================================

# Re-read dose response data to get individual gene values
dose_response_raw <- read_tsv(file = "data/mRNA Expression, RSEM (Batch normalized from Illumina HiSeq_RNASeqV2)_Tcell_dose_response.txt",
                              col_names = T) %>%
  pivot_longer(cols = 3:last_col(),
               names_to = "gene",
               values_to = "expression") %>%
  drop_na()

# List of genes to analyze individually
genes_to_plot <- c("SOD2", "CD274", "CSF3", "IDO1")

for (gene_name in genes_to_plot) {
  
  cat("\n===== Processing gene:", gene_name, "=====\n")
  
  # Get individual gene expression
  gene_expression_data <- dose_response_raw %>%
    filter(gene == gene_name) %>%
    select(STUDY_ID, SAMPLE_ID, expression) %>%
    rename(gene_expression = expression)
  
  # Merge with all_data
  all_data_gene <- all_data %>%
    left_join(gene_expression_data,
              by = c("Sample ID" = "SAMPLE_ID",
                     "Study ID" = "STUDY_ID")) %>%
    filter(!is.na(gene_expression)) %>%
    mutate(log2_gene_expr = log2(gene_expression + 1))
  
  # Pancancer violin plot
  # Wilcoxon test for pancancer individual gene
  wt_gene_data <- all_data_gene %>% filter(PDGFRA_cn_alteration == "WT") %>% pull(log2_gene_expr)
  ampl_gene_data <- all_data_gene %>% filter(PDGFRA_cn_alteration == "Amplified") %>% pull(log2_gene_expr)
  wilcox_test_gene_pancancer <- wilcox.test(wt_gene_data, ampl_gene_data)
  p_val_gene_pancancer <- wilcox_test_gene_pancancer$p.value
  
  # Add sample counts
  sample_counts_gene_pan <- all_data_gene %>%
    group_by(PDGFRA_cn_alteration) %>%
    summarise(n = n(), .groups = 'drop')
  
  all_data_gene_labeled <- all_data_gene %>%
    left_join(sample_counts_gene_pan, by = "PDGFRA_cn_alteration") %>%
    mutate(label = paste0(PDGFRA_cn_alteration, "\n(n = ", n, ")"))
  
  p_pancancer <- ggplot(all_data_gene_labeled, aes(x = factor(label, levels = unique(label[order(PDGFRA_cn_alteration)])), 
                                           y = log2_gene_expr,
                                           fill = PDGFRA_cn_alteration)) +
    geom_boxplot(width = 0.4, alpha = 0.9, outlier.shape = NA) +
    geom_jitter(width = 0.15, alpha = 0.5, size = 1.5) +
    scale_fill_manual(values = c("WT" = "#edf2f4", "Amplified" = "#8d99ae"),
                      guide = "none") +
    xlab("PDGFRA Status") +
    ylab(bquote(.(gene_name)~Expression~(log[2](TPM + 1)))) +
      ggtitle(paste0("Pan Cancer: ", gene_name, " Expression by PDGFRA CN Status")) +
    annotate(geom = "text",
             x = 1.5, y = Inf,
             label = paste0("p = ", format(p_val_gene_pancancer, digits = 3)),
             vjust = 1.5, hjust = 0.5,
             size = 3,
             family = .plot_font_family) +
    theme_minimal() +
    theme(axis.text = element_text(size = 10, family = .plot_font_family),
          axis.title = element_text(size = 11, family = .plot_font_family),
        plot.title = element_text(hjust = 0.5, face = "bold", family = .plot_font_family),
          plot.margin = margin(t = 10, r = 10, b = 10, l = 10),
          panel.grid = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
          axis.ticks = element_line(colour = "black", linewidth = 0.3))
  
  ggsave(paste0("figures/", gene_name, "_pancancer_violin.png"), 
         p_pancancer, dpi = 300, width = 6, height = 6)
  
  # GBM analysis
  gbm_data_gene <- gbm_data %>%
    left_join(gene_expression_data,
              by = c("Sample ID" = "SAMPLE_ID",
                     "Study ID" = "STUDY_ID")) %>%
    filter(!is.na(gene_expression)) %>%
    mutate(log2_gene_expr = log2(gene_expression + 1))
  
  # GBM violin plot
  # Wilcoxon test for GBM individual gene
  wt_gene_data_gbm <- gbm_data_gene %>% filter(PDGFRA_cn_alteration == "WT") %>% pull(log2_gene_expr)
  ampl_gene_data_gbm <- gbm_data_gene %>% filter(PDGFRA_cn_alteration == "Amplified") %>% pull(log2_gene_expr)
  wilcox_test_gene_gbm <- wilcox.test(wt_gene_data_gbm, ampl_gene_data_gbm)
  p_val_gene_gbm <- wilcox_test_gene_gbm$p.value
  
  # Add sample counts
  sample_counts_gene_gbm <- gbm_data_gene %>%
    group_by(PDGFRA_cn_alteration) %>%
    summarise(n = n(), .groups = 'drop')
  
  gbm_data_gene_labeled <- gbm_data_gene %>%
    left_join(sample_counts_gene_gbm, by = "PDGFRA_cn_alteration") %>%
    mutate(label = paste0(PDGFRA_cn_alteration, "\n(n = ", n, ")"))
  
  p_gbm <- ggplot(gbm_data_gene_labeled, aes(x = factor(label, levels = unique(label[order(PDGFRA_cn_alteration)])), 
                                     y = log2_gene_expr,
                                     fill = PDGFRA_cn_alteration)) +
    geom_boxplot(width = 0.4, alpha = 0.9, outlier.shape = NA) +
    geom_jitter(width = 0.15, alpha = 0.5, size = 1.5) +
    scale_fill_manual(values = c("WT" = "#edf2f4", "Amplified" = "#8d99ae"),
                      guide = "none") +
    xlab("PDGFRA Status") +
    ylab(bquote(.(gene_name)~Expression~(log[2](TPM + 1)))) +
      ggtitle(paste0("GBM (IDHwt): ", gene_name, " Expression by PDGFRA CN Status")) +
    annotate(geom = "text",
             x = 1.5, y = Inf,
             label = paste0("p = ", format(p_val_gene_gbm, digits = 3)),
             vjust = 1.5, hjust = 0.5,
             size = 3,
             family = .plot_font_family) +
    theme_minimal() +
    theme(axis.text = element_text(size = 10, family = .plot_font_family),
          axis.title = element_text(size = 11, family = .plot_font_family),
        plot.title = element_text(hjust = 0.5, face = "bold", family = .plot_font_family),
          plot.margin = margin(t = 10, r = 10, b = 10, l = 10),
          panel.grid = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
          axis.ticks = element_line(colour = "black", linewidth = 0.3))
  
  ggsave(paste0("figures/", gene_name, "_gbm_violin.png"), 
         p_gbm, dpi = 300, width = 6, height = 6)
  
  # GBM scatter plot with correlation
  if (nrow(gbm_data_gene) > 0) {
    corr_gene <- cor.test(gbm_data_gene$log2_PDGFRA_exprs, 
                          gbm_data_gene$log2_gene_expr,
                          method = "pearson")
    
    p_scatter <- ggplot(gbm_data_gene,
                        aes(x = log2_PDGFRA_exprs,
                            y = log2_gene_expr)) +
      geom_point(aes(color = PDGFRA_cn_alteration),
                 size = 2, alpha = 0.6) +
      geom_smooth(method = "lm",
                  se = TRUE, color = "black", alpha = 0.2) +
      xlab("Log2 PDGFRA Expression") +
      ylab(bquote(.(gene_name)~Expression~(log[2](TPM + 1)))) +
      annotate(geom = "text",
               x = Inf, y = -Inf,
               label = paste0("r = ", round(corr_gene$estimate, 3), 
                             "\np = ", format(corr_gene$p.value, digits = 3)),
               vjust = -0.5, hjust = 1.1,
               size = 3) +
      scale_color_manual(values = c("WT" = "#999999", "Amplified" = "#FF6633"),
                        name = "CN Status") +
      theme_classic() +
      theme(legend.position = "bottom",
            text = element_text(size = 10),
            plot.title = element_text(hjust = 0.5, face = "bold")) +
      labs(title = paste0(gene_name, " vs PDGFRA Expression (GBM)"))
    
    # Marginal plot
    p_scatter_marginal <- ggExtra::ggMarginal(p_scatter, 
                                              type = "density",
                                              groupColour = TRUE, 
                                              groupFill = TRUE, 
                                              alpha = 0.25, size = 4)
    
    ggsave(paste0("figures/", gene_name, "_gbm_scatter_marginal.png"), 
           p_scatter_marginal, dpi = 300, width = 6, height =5)
  }
  
  cat("Completed:", gene_name, "\n")
}

# ==========================================================================
# INDIVIDUAL GENE ANALYSIS: GLIOBLASTOMA + GLIOMA COMBINED
# ==========================================================================

for (gene_name in genes_to_plot) {
  
  cat("\n===== Processing gene (Glioblastoma + Glioma):", gene_name, "=====\n")
  
  # Get individual gene expression
  gene_expression_data <- dose_response_raw %>%
    filter(gene == gene_name) %>%
    select(STUDY_ID, SAMPLE_ID, expression) %>%
    rename(gene_expression = expression)
  
  # Glioblastoma + Glioma analysis
  glioma_data_gene <- glioma_data %>%
    left_join(gene_expression_data,
              by = c("Sample ID" = "SAMPLE_ID",
                     "Study ID" = "STUDY_ID")) %>%
    filter(!is.na(gene_expression)) %>%
    mutate(log2_gene_expr = log2(gene_expression + 1))
  
  # Glioblastoma + Glioma violin plot
  # Wilcoxon test for combined cohort
  wt_gene_data_glioma <- glioma_data_gene %>% filter(PDGFRA_cn_alteration == "WT") %>% pull(log2_gene_expr)
  ampl_gene_data_glioma <- glioma_data_gene %>% filter(PDGFRA_cn_alteration == "Amplified") %>% pull(log2_gene_expr)
  wilcox_test_gene_glioma <- wilcox.test(wt_gene_data_glioma, ampl_gene_data_glioma)
  p_val_gene_glioma <- wilcox_test_gene_glioma$p.value
  
  # Add sample counts
  sample_counts_gene_glioma <- glioma_data_gene %>%
    group_by(PDGFRA_cn_alteration) %>%
    summarise(n = n(), .groups = 'drop')
  
  glioma_data_gene_labeled <- glioma_data_gene %>%
    left_join(sample_counts_gene_glioma, by = "PDGFRA_cn_alteration") %>%
    mutate(label = paste0(PDGFRA_cn_alteration, "\n(n = ", n, ")"))
  
  p_glioma_gene <- ggplot(glioma_data_gene_labeled, aes(x = factor(label, levels = unique(label[order(PDGFRA_cn_alteration)])), 
                                                y = log2_gene_expr,
                                                fill = PDGFRA_cn_alteration)) +
    geom_boxplot(width = 0.4, alpha = 0.9, outlier.shape = NA) +
    geom_jitter(width = 0.15, alpha = 0.5, size = 1.5) +
    scale_fill_manual(values = c("WT" = "#edf2f4", "Amplified" = "#8d99ae"),
                      guide = "none") +
    xlab("PDGFRA Status") +
    ylab(bquote(.(gene_name)~Expression~(log[2](TPM + 1)))) +
      ggtitle(paste0("Glioblastoma + Glioma: ", gene_name, " Expression by PDGFRA CN Status")) +
    annotate(geom = "text",
             x = 1.5, y = Inf,
             label = paste0("p = ", format(p_val_gene_glioma, digits = 3)),
             vjust = 1.5, hjust = 0.5,
             size = 3,
             family = .plot_font_family) +
    theme_minimal() +
    theme(axis.text = element_text(size = 10, family = .plot_font_family),
          axis.title = element_text(size = 11, family = .plot_font_family),
        plot.title = element_text(hjust = 0.5, face = "bold", family = .plot_font_family),
          plot.margin = margin(t = 10, r = 10, b = 10, l = 10),
          panel.grid = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
          axis.ticks = element_line(colour = "black", linewidth = 0.3))
  
  ggsave(paste0("figures/", gene_name, "_glioblastoma_glioma_violin.png"), 
         p_glioma_gene, dpi = 300, width = 6, height = 6)
  
  # Glioblastoma + Glioma scatter plot with correlation
  if (nrow(glioma_data_gene) > 0) {
    corr_gene_glioma <- cor.test(glioma_data_gene$log2_PDGFRA_exprs, 
                                 glioma_data_gene$log2_gene_expr,
                                 method = "pearson")
    
    p_scatter_glioma <- ggplot(glioma_data_gene,
                               aes(x = log2_PDGFRA_exprs,
                                   y = log2_gene_expr)) +
      geom_point(aes(color = PDGFRA_cn_alteration),
                 size = 2, alpha = 0.6) +
      geom_smooth(method = "lm",
                  se = TRUE, color = "black", alpha = 0.2) +
      xlab("Log2 PDGFRA Expression") +
      ylab(bquote(.(gene_name)~Expression~(log[2](TPM + 1)))) +
      annotate(geom = "text",
               x = Inf, y = -Inf,
               label = paste0("r = ", round(corr_gene_glioma$estimate, 3), 
                             "\np = ", format(corr_gene_glioma$p.value, digits = 3)),
               vjust = -0.5, hjust = 1.1,
               size = 3) +
      scale_color_manual(values = c("WT" = "#999999", "Amplified" = "#FF6633"),
                        name = "CN Status") +
      theme_classic() +
      theme(legend.position = "bottom",
            text = element_text(size = 10),
            plot.title = element_text(hjust = 0.5, face = "bold")) +
      labs(title = paste0(gene_name, " vs PDGFRA Expression (Glioblastoma + Glioma)"))
    
    # Marginal plot
    p_scatter_glioma_marginal <- ggExtra::ggMarginal(p_scatter_glioma, 
                                                     type = "density",
                                                     groupColour = TRUE, 
                                                     groupFill = TRUE, 
                                                     alpha = 0.25, size = 4)
    
    ggsave(paste0("figures/", gene_name, "_glioblastoma_glioma_scatter_marginal.png"), 
           p_scatter_glioma_marginal, dpi = 300, width = 6, height = 5)
  }
  
  cat("Completed (Glioblastoma + Glioma):", gene_name, "\n")
}

cat("\n===== Individual Gene Analysis Complete =====\n")
