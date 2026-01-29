#pert_type <-"crispri" 
library(devtools)
library(plyr)
library(dplyr)
library(monocle3)
library(ggplot2)
library(ggridges)
library(ggrepel)
#library(ggpubr)
library(gridExtra)
library(pheatmap)
library(tibble)
library(tidyr)
library(parallel)
library(DelayedArray)
library(Matrix)
library(magrittr)
library(purrr)
library(furrr)
library(future)
library(progressr)
library(cli)

pert_type <- "crispra"

base_dir <- paste0("/home/user/Documents/Kinase_project/", pert_type, "_final/")
setwd(base_dir)

output_folder <- paste0(toupper(pert_type), "_Figures/")



# Output path for QC images
output_path <- file.path(base_dir, output_folder)
# make the dir if it doesn't exist
if (!dir.exists(output_path)) {
  dir.create(output_path, recursive = TRUE)
}

cds_folder <- "preprocessed_cds_without_Tcells.rds"
cds_path <- file.path(base_dir, cds_folder)
cds <- readRDS(cds_path)
cds <- cds[, !is.na(colData(cds)$top_sg)]
dim(cds) # 58347 166061


genotype_by_Tcell_dose_diff_test.list <- readRDS(file.path(output_path, "Diff_test/genotype_by_Tcell_dose_diff_test.list.rds"))


kinases <- (unique(colData(cds)$gene_id))
kinases <- kinases[kinases != "NTC"]
kinases <- na.omit(kinases)


for (dose in names(genotype_by_Tcell_dose_diff_test.list)) {
  for (target in kinases) {
    genotype_by_Tcell_dose_diff_test.list[[dose]][[target]]$gene_id <- rep(target, nrow(genotype_by_Tcell_dose_diff_test.list[[dose]][[target]]))
  }
}

for (dose in names(genotype_by_Tcell_dose_diff_test.list)) {
  genotype_by_Tcell_dose_diff_test.list[[dose]] <- do.call("rbind", genotype_by_Tcell_dose_diff_test.list[[dose]])
  genotype_by_Tcell_dose_diff_test.list[[dose]]$treatment <- rep(dose, nrow(genotype_by_Tcell_dose_diff_test.list[[dose]]))
}

genotype_by_Tcell_dose_diff_test_results <- do.call("rbind", genotype_by_Tcell_dose_diff_test.list)

genotype_by_Tcell_dose_diff_test_results$q_value <- p.adjust(genotype_by_Tcell_dose_diff_test_results$p_value, method = "BH")



Tcell_dose_diff_test <- readRDS(file.path("/home/user/Documents/Kinase_project/crispra_final/CRISPRA_Figures/NTC_analysis/DEG_testing/NTC_CRISPR_dose_response_diff_test_replicates.rds"))
Tcell_dose_diff_test <- Tcell_dose_diff_test %>%
  filter(!is.na(p_value)) %>%
  mutate(q_value = p.adjust(p_value, method = "BH"))

Tcell_dose_response_diff_test_top_deg <- Tcell_dose_diff_test[Tcell_dose_diff_test$q_value < 0.05 & Tcell_dose_diff_test$term == "dose", ]
Tcell_dose_response_diff_test_top_deg <- Tcell_dose_response_diff_test_top_deg[order(Tcell_dose_response_diff_test_top_deg$q_value), ]

top50 <- head(Tcell_dose_response_diff_test_top_deg, 100)
top_Tcell_dose_gene <- top50[top50$normalized_effect > 1.2, ]

sod2_csf3_modulating_kinases <- genotype_by_Tcell_dose_diff_test_results %>%
  filter(gene_short_name %in% c(top_Tcell_dose_gene$gene_short_name, "CD274", "IDO1"), term != "(Intercept)", q_value < 0.05 & abs(normalized_effect) > 1) %>%
  arrange(desc(normalized_effect)) %>%
  group_by(gene_short_name) %>%
  dplyr::slice(1:10) %>%
  dplyr::select(gene_short_name, normalized_effect, term)


sod2_csf3_modulating_kinases <- genotype_by_Tcell_dose_diff_test_results %>%
  filter(gene_short_name %in% c(top_Tcell_dose_gene$gene_short_name, "CD274", "IDO1"), term != "(Intercept)", q_value < 0.05 & abs(normalized_effect) > 1) %>%
  arrange(desc(normalized_effect)) %>%
  group_by(gene_short_name) %>%
  dplyr::slice(1:10) %>%
  dplyr::select(gene_short_name, normalized_effect, term)

sod2_csf3_modulating_kinases$gene_id <- sapply(sod2_csf3_modulating_kinases$term, function(x) {
  stringr::str_split(x, pattern = "_id")[[1]][2]
})

sod2_csf3_modulating_kinases <- sod2_csf3_modulating_kinases %>%
  dplyr::distinct(gene_short_name, gene_id, .keep_all = TRUE)
write.csv(sod2_csf3_modulating_kinases, file.path(output_path, "Heatmaps/SOD2_CSF3_modulators_top_Tcell_dose_CD274_IDO1_100.csv"), row.names = FALSE)


sod2_csf3_modulating_kinases_beta_heatmap <- sod2_csf3_modulating_kinases %>%
  dplyr::select(gene_short_name, normalized_effect, gene_id) %>%
  tidyr::spread(key = gene_id, value = normalized_effect) %>%
  as.data.frame()

row.names(sod2_csf3_modulating_kinases_beta_heatmap) <- sod2_csf3_modulating_kinases_beta_heatmap$gene_short_name
sod2_csf3_modulating_kinases_beta_heatmap$gene_short_name <- NULL

sod2_csf3_modulating_kinases_beta_heatmap[is.na(sod2_csf3_modulating_kinases_beta_heatmap)] <- 0

hmcols <- colorRampPalette(colors = c("blue", "white", "red"))(35)
paletteLength <- 35
myBreaks <- c(seq(min(sod2_csf3_modulating_kinases_beta_heatmap), 0, length.out = ceiling(paletteLength / 2) + 1), seq(max(sod2_csf3_modulating_kinases_beta_heatmap) / floor(paletteLength / 2), max(sod2_csf3_modulating_kinases_beta_heatmap), length.out = floor(paletteLength / 2)))

pheatmap::pheatmap(sod2_csf3_modulating_kinases_beta_heatmap,
  clustering_method = "ward.D2",
  cluster_rows = F,
  color = hmcols,
  breaks = myBreaks,
  file = file.path(output_path, "Heatmaps/SOD2_CSF3_modulators_top_Tcell_dose_CD274_IDO1.png"),
  fontsize = 9,
  treeheight_col = 10,
  height = 3,
  width = 10
)

pheatmap::pheatmap(sod2_csf3_modulating_kinases_beta_heatmap,
  clustering_method = "ward.D2",
  cluster_rows = F,
  color = hmcols,
  breaks = myBreaks,
  file = file.path(output_path, "Heatmaps/SOD2_CSF3_modulators_top_Tcell_dose_CD274_IDO1.pdf"),
  fontsize = 9,
  treeheight_col = 10,
  height = 3,
  width = 10
)
