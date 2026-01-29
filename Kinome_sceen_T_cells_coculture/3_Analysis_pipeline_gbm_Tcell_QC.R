# 3. Analysis_pipeline_crisprI_gbm_Tcell_QC_20250104

# 0. Setup and Library Imports
rm(list = ls())
gc()

library(devtools)
library(parallel)
library(ggplot2)
library(ggridges)
library(gridExtra)
library(plyr)
library(dplyr)
library(tidyr)
library(tibble)
library(pheatmap)
# library(piano)
library(DelayedArray)
library(monocle3)
library(UpSetR)

source("/home/user/Documents/github_repos/sci-plex/bin/cell_cycle.R")
source("/home/user/Documents/github_repos/sci-plex/bin/dispersions_functions.R")
cc.genes <- readRDS("/home/user/Documents/github_repos/sci-plex/bin/cc.genes.RDS")

protein_to_gene_map <- readRDS("/home/user/Documents/github_repos/sci-Plex-GxE/Kinome_GxE_screen/KinMap_kinase_protein_to_gene_map.rds")

# 1. Functions for comparing knockdown effect
calculate_mean_target_knockdown <- function(cds, gene_short_name_list, protein_to_gene_map) {
  cds_subset <- cds[protein_to_gene_map$id, ]
  cds_exprs <- SingleCellExperiment::counts(cds_subset)
  cds_exprs <- Matrix::t(Matrix::t(cds_exprs) / monocle3:::size_factors(cds_subset))
  cds_exprs <- reshape2::melt(as.matrix(cds_exprs))

  colnames(cds_exprs) <- c("id", "cell", "expression")
  cds_exprs$id <- as.character(cds_exprs$id)
  cds_exprs$cell <- as.character(cds_exprs$cell)

  # protein_to_gene_map contains ensembl id, HGNC gene_short_name and the protein name
  cds_exprs <- dplyr::left_join(cds_exprs, protein_to_gene_map, by = "id")

  Mean_expression_by_target <- list()

  # Clean this up to a map function
  targets <- colData(cds) %>%
    as.data.frame() %>%
    dplyr::filter(!(gene_id %in% c("NTC", "random"))) %>%
    dplyr::pull(gene_id) %>%
    unique() %>%
    sort()

  for (target in targets) {
    target_cells <- colData(cds) %>%
      as.data.frame() %>%
      dplyr::filter(gene_id == target) %>%
      dplyr::pull(cell)

    target_cells_2_cutoff <- colData(cds) %>%
      as.data.frame() %>%
      dplyr::filter(gene_id == target & total_sgrna_read_per_cell >= 2) %>%
      dplyr::pull(cell)

    target_cells_5_cutoff <- colData(cds) %>%
      as.data.frame() %>%
      dplyr::filter(gene_id == target & total_sgrna_read_per_cell >= 5) %>%
      dplyr::pull(cell)

    target_cells_10_cutoff <- colData(cds) %>%
      as.data.frame() %>%
      dplyr::filter(gene_id == target & total_sgrna_read_per_cell >= 10) %>%
      dplyr::pull(cell)

    target_exprs <- cds_exprs %>%
      dplyr::filter(gene_id == target &
        cell %in% target_cells) %>%
      dplyr::pull(expression) %>%
      mean()

    target_exprs_2_cutoff <- cds_exprs %>%
      dplyr::filter(gene_id == target &
        cell %in% target_cells_2_cutoff) %>%
      dplyr::pull(expression) %>%
      mean()

    target_exprs_5_cutoff <- cds_exprs %>%
      dplyr::filter(gene_id == target &
        cell %in% target_cells_5_cutoff) %>%
      dplyr::pull(expression) %>%
      mean()

    target_exprs_10_cutoff <- cds_exprs %>%
      dplyr::filter(gene_id == target &
        cell %in% target_cells_10_cutoff) %>%
      dplyr::pull(expression) %>%
      mean()

    NTC_cells <- colData(cds) %>%
      as.data.frame() %>%
      dplyr::filter(gene_id == "NTC") %>%
      dplyr::pull(cell)

    NTC_cells_2_cutoff <- colData(cds) %>%
      as.data.frame() %>%
      dplyr::filter(gene_id == "NTC" & total_sgrna_read_per_cell >= 2) %>%
      dplyr::pull(cell)

    NTC_cells_5_cutoff <- colData(cds) %>%
      as.data.frame() %>%
      dplyr::filter(gene_id == "NTC" & total_sgrna_read_per_cell >= 5) %>%
      dplyr::pull(cell)

    NTC_cells_10_cutoff <- colData(cds) %>%
      as.data.frame() %>%
      dplyr::filter(gene_id == "NTC" & total_sgrna_read_per_cell >= 10) %>%
      dplyr::pull(cell)

    NTC_exprs <- cds_exprs %>%
      dplyr::filter(gene_id == target &
        cell %in% NTC_cells) %>%
      dplyr::pull(expression) %>%
      mean()

    NTC_exprs_2_cutoff <- cds_exprs %>%
      dplyr::filter(gene_id == target &
        cell %in% NTC_cells_2_cutoff) %>%
      dplyr::pull(expression) %>%
      mean()

    NTC_exprs_5_cutoff <- cds_exprs %>%
      dplyr::filter(gene_id == target &
        cell %in% NTC_cells_5_cutoff) %>%
      dplyr::pull(expression) %>%
      mean()

    NTC_exprs_10_cutoff <- cds_exprs %>%
      dplyr::filter(gene_id == target &
        cell %in% NTC_cells_10_cutoff) %>%
      dplyr::pull(expression) %>%
      mean()

    knockdown <- target_exprs / NTC_exprs
    knockdown_2_cutoff <- target_exprs_2_cutoff / NTC_exprs_2_cutoff
    knockdown_5_cutoff <- target_exprs_5_cutoff / NTC_exprs_5_cutoff
    knockdown_10_cutoff <- target_exprs_10_cutoff / NTC_exprs_10_cutoff

    set.seed(2016L)
    target_random_cells <- colData(cds) %>%
      as.data.frame() %>%
      dplyr::sample_n(length(target_cells)) %>%
      dplyr::pull(cell) %>%
      as.character()

    set.seed(2016L)
    target_random_cells_2_cutoff <- colData(cds) %>%
      as.data.frame() %>%
      dplyr::filter(total_sgrna_read_per_cell >= 2) %>%
      dplyr::sample_n(length(target_cells_2_cutoff)) %>%
      dplyr::pull(cell) %>%
      as.character()

    set.seed(2016L)
    target_random_cells_5_cutoff <- colData(cds) %>%
      as.data.frame() %>%
      dplyr::filter(total_sgrna_read_per_cell >= 5) %>%
      dplyr::sample_n(length(target_cells_5_cutoff)) %>%
      dplyr::pull(cell) %>%
      as.character()

    set.seed(2016L)
    target_random_cells_10_cutoff <- colData(cds) %>%
      as.data.frame() %>%
      dplyr::filter(total_sgrna_read_per_cell >= 10) %>%
      dplyr::sample_n(length(target_cells_10_cutoff)) %>%
      dplyr::pull(cell) %>%
      as.character()

    target_random_exprs <- cds_exprs %>%
      dplyr::filter(gene_id == target &
        cell %in% target_random_cells) %>%
      dplyr::pull(expression) %>%
      mean()

    target_random_exprs_2_cutoff <- cds_exprs %>%
      dplyr::filter(gene_id == target &
        cell %in% target_random_cells_2_cutoff) %>%
      dplyr::pull(expression) %>%
      mean()

    target_random_exprs_5_cutoff <- cds_exprs %>%
      dplyr::filter(gene_id == target &
        cell %in% target_random_cells_5_cutoff) %>%
      dplyr::pull(expression) %>%
      mean()

    target_random_exprs_10_cutoff <- cds_exprs %>%
      dplyr::filter(gene_id == target &
        cell %in% target_random_cells_10_cutoff) %>%
      dplyr::pull(expression) %>%
      mean()

    set.seed(2016L)
    NTC_random_cells <- colData(cds) %>%
      as.data.frame() %>%
      dplyr::sample_n(length(NTC_cells)) %>%
      dplyr::pull(cell) %>%
      as.character()

    set.seed(2016L)
    NTC_random_cells_2_cutoff <- colData(cds) %>%
      as.data.frame() %>%
      dplyr::filter(total_sgrna_read_per_cell >= 2) %>%
      dplyr::sample_n(length(NTC_cells_2_cutoff)) %>%
      dplyr::pull(cell) %>%
      as.character()

    set.seed(2016L)
    NTC_random_cells_5_cutoff <- colData(cds) %>%
      as.data.frame() %>%
      dplyr::filter(total_sgrna_read_per_cell >= 5) %>%
      dplyr::sample_n(length(NTC_cells_5_cutoff)) %>%
      dplyr::pull(cell) %>%
      as.character()

    set.seed(2016L)
    NTC_random_cells_10_cutoff <- colData(cds) %>%
      as.data.frame() %>%
      dplyr::filter(total_sgrna_read_per_cell >= 10) %>%
      dplyr::sample_n(length(NTC_cells_10_cutoff)) %>%
      dplyr::pull(cell) %>%
      as.character()

    NTC_random_exprs <- cds_exprs %>%
      dplyr::filter(gene_id == target &
        cell %in% NTC_random_cells) %>%
      dplyr::pull(expression) %>%
      mean()

    NTC_random_exprs_2_cutoff <- cds_exprs %>%
      dplyr::filter(gene_id == target &
        cell %in% NTC_random_cells_2_cutoff) %>%
      dplyr::pull(expression) %>%
      mean()

    NTC_random_exprs_5_cutoff <- cds_exprs %>%
      dplyr::filter(gene_id == target &
        cell %in% NTC_random_cells_5_cutoff) %>%
      dplyr::pull(expression) %>%
      mean()

    NTC_random_exprs_10_cutoff <- cds_exprs %>%
      dplyr::filter(gene_id == target &
        cell %in% NTC_random_cells_10_cutoff) %>%
      dplyr::pull(expression) %>%
      mean()

    random_knockdown <- target_random_exprs / NTC_random_exprs
    random_knockdown_2_cutoff <- target_random_exprs_2_cutoff / NTC_random_exprs_2_cutoff
    random_knockdown_5_cutoff <- target_random_exprs_5_cutoff / NTC_random_exprs_5_cutoff
    random_knockdown_10_cutoff <- target_random_exprs_10_cutoff / NTC_random_exprs_10_cutoff

    Mean_expression_by_target[[target]] <- data.frame(
      gene_id = target,
      mean_knockdown = knockdown,
      mean_random_knockdown = random_knockdown,
      mean_knockdown_2_cutoff = knockdown_2_cutoff,
      mean_random_knockdown_2_cutoff = random_knockdown_2_cutoff,
      mean_knockdown_5_cutoff = knockdown_5_cutoff,
      mean_random_knockdown_5_cutoff = random_knockdown_5_cutoff,
      mean_knockdown_10_cutoff = knockdown_10_cutoff,
      mean_random_knockdown_10_cutoff = random_knockdown_10_cutoff
    )

    message("finished ", target)
  }

  Mean_expression <- do.call("rbind", Mean_expression_by_target)

  return(Mean_expression)
}

# 2. Load in the data
pert_type <- "crispri"

base_dir <- paste0("/home/user/Documents/Kinase_project/", pert_type, "_final/")
setwd(base_dir)

output_folder <- paste0(toupper(pert_type), "_Figures")

# Output path for QC images
output_path <- file.path(base_dir, output_folder)

cds_folder <- "11-final-output/CDS_crop_hash_FINAL_CRISPR_GBM_T_cells_500_cutoff.rds"
cds_path <- file.path(base_dir, cds_folder)
cds <- readRDS(cds_path)
print(dim(cds))  # 58347 200473

# Remove doublets
doublet_cells <- read.csv("11-final-output/doublet_info.csv")
doublet_cells <- doublet_cells %>%
  dplyr::filter(is_doublet == "True") %>%
  dplyr::pull(cell_ID)
cds <- cds[, !colData(cds)$cell_ID %in% doublet_cells]
print(dim(cds))  # 58347 200436

# Remove cells with NA sgRNA
cds <- cds[, !is.na(colData(cds)$sgRNA)]
print(dim(cds))  # 58347 191805
output_path <- file.path(base_dir, output_folder)

print(stats::median(colData(cds)$n.umi))  # 2681

print(head(colData(cds)))

cds_data <- colData(cds) %>% as.data.frame()

hash_cds <- cds[, !is.na(colData(cds)$treatment)]

# 3. Evaluate sgRNA efficiency
# 3.1 Clean up guide names
colData(cds)$sgRNA <- sapply(colData(cds)$sgRNA, function(x) {
  if (grepl("negative_control", x)) {
    return("NTC")
  }
  if (grepl("non-targeting", x)) {
    return("NTC")
  }
  if (grepl("originalIDTorder", x)) {
    return("random")
  }
  if (grepl("randomregion", x)) {
    return("random")
  }
  return(x)
})

# Clean up gene names
cds <- cds[, !is.na(colData(cds)$sgRNA)]
colData(cds)$gene <- sapply(colData(cds)$sgRNA, function(x) {
  stringr::str_split(x, pattern = "_")[[1]][1]
})

colData(cds)$second_sg <- sapply(colData(cds)$second_sg, function(x) {
  if (grepl("negative_control", x)) {
    return("NTC")
  }
  if (grepl("non-targeting", x)) {
    return("NTC")
  }
  if (grepl("originalIDTorder", x)) {
    return("random")
  }
  if (grepl("randomregion", x)) {
    return("random")
  }
  return(x)
})

colData(cds)$third_sg <- sapply(colData(cds)$third_sg, function(x) {
  if (grepl("negative_control", x)) {
    return("NTC")
  }
  if (grepl("non-targeting", x)) {
    return("NTC")
  }
  if (grepl("originalIDTorder", x)) {
    return("random")
  }
  if (grepl("randomregion", x)) {
    return("random")
  }
  return(x)
})

colData(cds)$second_sg_gene <- sapply(colData(cds)$second_sg, function(x) {
  stringr::str_split(x, pattern = "_")[[1]][1]
})

colData(cds)$third_sg_gene <- sapply(colData(cds)$third_sg, function(x) {
  stringr::str_split(x, pattern = "_")[[1]][1]
})

# 3.2. Compare sgRNA counts in untreated and T cell treated



# guide distribution relative to plasmid
library(stringr)
library(dplyr)
if (pert_type == 'crispri') {
  U6_H05H05_R2 <- read.table("/home/user/Documents/Kinase_project/2024_guide_sequencing_crisprianda//H05H05_R2_containsU6.txt", header = FALSE, sep = "\t")
  U6_H05H05_R2_lengths <- U6_H05H05_R2
  U6_H05H05_R2_lengths$length_of_sequence <- str_length(U6_H05H05_R2_lengths$V1)
  U6_H05H05_R2_sgRNA_extracted <- U6_H05H05_R2
  U6_H05H05_R2_sgRNA_extracted$extracted_sgRNA_sequence <- substr(U6_H05H05_R2_sgRNA_extracted$V1, 24, 42)
  U6_H05H05_R2_sgRNA_extracted_to_merge <- U6_H05H05_R2_sgRNA_extracted
  U6_H05H05_R2_sgRNA_extracted_to_merge$sgRNA.sequence <- U6_H05H05_R2_sgRNA_extracted_to_merge$extracted_sgRNA_sequence
  CRISPRI <-read.table("/home/user/Documents/Kinase_project/2024_guide_sequencing_crisprianda//sciPlexGxE_2_gRNASampleSheet_CRISPRI.txt")
  colnames(CRISPRI) <- c('guide_name','sgRNA.sequence','V3')
  CRISPRI <-CRISPRI %>% mutate(sgrna_length = nchar(sgRNA.sequence))
  sgRNAs_in_library_H05H05 <- inner_join(U6_H05H05_R2_sgRNA_extracted_to_merge, CRISPRI, by = "sgRNA.sequence")
  head(sgRNAs_in_library_H05H05)

  #sgRNAs_in_library_H05H05$gene <- sapply(sgRNAs_in_library_H05H05$guide_name, function(x){stringr::str_split(x,pattern = "_")[[1]][1]})
  library(stringr)

  # Using str_extract to directly extract the gene name
  sgRNAs_in_library_H05H05$gene <- str_extract(sgRNAs_in_library_H05H05$guide_name, "^[^_]+")

  plasmid_sgRNA_distribution_summary_df <- sgRNAs_in_library_H05H05 %>%
    group_by(gene) %>%
    dplyr::mutate(plasmid_count = n()) %>%
    ungroup() %>%
    dplyr::mutate(total_count = n()) %>%
    dplyr::select(gene,plasmid_count, total_count) %>%
    distinct() %>%
    mutate(plasmid_frequency = (plasmid_count/total_count)*100)

  CRISPRI_sgRNA_distribution_summary_df <- colData(cds) %>%
    as.data.frame() %>%
    group_by(gene) %>%
    dplyr::mutate(CRISPRI_count = n()) %>%
    ungroup() %>%
    dplyr::mutate(total_CRISPRI_count = n()) %>%
    dplyr::select(gene,CRISPRI_count, total_CRISPRI_count) %>%
    distinct() %>%
    mutate(CRISPRI_frequency = (CRISPRI_count/total_CRISPRI_count)*100)

  CRISPRI_sgRNA_proportions <- inner_join(CRISPRI_sgRNA_distribution_summary_df,plasmid_sgRNA_distribution_summary_df,by = "gene")
  CRISPRI_sgRNA_proportions$relative_proportion <- log(CRISPRI_sgRNA_proportions$CRISPRI_frequency/CRISPRI_sgRNA_proportions$plasmid_frequency) %>% scale()
  CRISPRI_rank <- CRISPRI_sgRNA_proportions %>% arrange(relative_proportion) %>% pull(gene)

  ggplot(CRISPRI_sgRNA_proportions, aes(x = factor(gene, levels = CRISPRI_rank), y = relative_proportion, label = gene, color =  abs(relative_proportion) > 2)) +
  geom_point(size = 1, stroke = 0) +
  ggrepel::geom_label_repel(data = subset(CRISPRI_sgRNA_proportions, abs(relative_proportion) > 2), 
                            size = 1, 
                            segment.size = 0.1,
                            label.size = 0.1,
                            box.padding = 0.1,
                            label.padding = 0.1,
                            color = "black") +
  monocle3:::monocle_theme_opts() +
  scale_color_manual(values = c("TRUE" = "firebrick3", "FALSE" = "black")) +
  theme(text = element_text(size = 6),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none") +
  ylab("Proportion relative to\n plasmid library") +
  xlab("Perturbation")
  ggsave(file.path(output_path, "/QC_plots/gRNA_distribution_relative_to_plasmid_CRISPRI_.png", height = 1.5, width = 2, dpi = 300))
  ggsave(file.path(output_path, "/QC_plots/gRNA_distribution_relative_to_plasmid_CRISPRI_.pdf", height = 1.5, width = 2, dpi = 300))
  ggsave(file.path(output_path, "/QC_plots/gRNA_distribution_relative_to_plasmid_CRISPRI_.svg", height = 1.5, width = 2, dpi = 300))
  saveRDS(CRISPRI_sgRNA_proportions, file = file.path(output_path, "/QC_plots/gRNA_distribution_relative_to_plasmid_CRISPRI_.rds"))
}  




Untreated_cds <- cds[, colData(cds)$treatment %in% c("Untreated")]
CRISPR_sgRNA_distribution_summary_df_untreated <- colData(Untreated_cds) %>%
  as.data.frame() %>%
  dplyr::group_by(gene) %>%
  dplyr::mutate(CRISPR_count = dplyr::n()) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(total_CRISPR_count = dplyr::n()) %>%
  dplyr::select(gene, CRISPR_count, total_CRISPR_count) %>%
  dplyr::distinct() %>%
  dplyr::mutate(CRISPR_frequency = (CRISPR_count / total_CRISPR_count) * 100)

# 3.2.1. 1:1
tcell_cds <- cds[, colData(cds)$treatment %in% c("1to1")]
CRISPR_sgRNA_distribution_summary_df_tcell <- colData(tcell_cds) %>%
  as.data.frame() %>%
  dplyr::group_by(gene) %>%
  dplyr::mutate(CRISPR_count = dplyr::n()) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(total_CRISPR_count = dplyr::n()) %>%
  dplyr::select(gene, CRISPR_count, total_CRISPR_count) %>%
  dplyr::distinct() %>%
  dplyr::mutate(CRISPR_frequency = (CRISPR_count / total_CRISPR_count) * 100)

CRISPR_sgRNA_proportions <- dplyr::inner_join(CRISPR_sgRNA_distribution_summary_df_untreated, CRISPR_sgRNA_distribution_summary_df_tcell, by = "gene")
CRISPR_sgRNA_proportions$relative_proportion <- log(CRISPR_sgRNA_proportions$CRISPR_frequency.x / CRISPR_sgRNA_proportions$CRISPR_frequency.y) %>% scale()
CRISPR_rank <- CRISPR_sgRNA_proportions %>%
  dplyr::arrange(relative_proportion) %>%
  dplyr::pull(gene)

ggplot2::ggplot(CRISPR_sgRNA_proportions, ggplot2::aes(x = factor(gene, levels = CRISPR_rank), y = relative_proportion, label = gene, color = abs(relative_proportion) > 2)) +
  ggplot2::geom_point(size = 1, stroke = 0) +
  ggrepel::geom_label_repel(
    data = subset(CRISPR_sgRNA_proportions, abs(relative_proportion) > 2),
    size = 1,
    segment.size = 0.1,
    label.size = 0.1,
    box.padding = 0,
    label.padding = 0.1,
    color = "black"
  ) +
  monocle3:::monocle_theme_opts() +
  ggplot2::scale_color_manual(values = c("TRUE" = "firebrick3", "FALSE" = "black")) +
  ggplot2::theme(
    text = ggplot2::element_text(size = 6),
    axis.text.x = ggplot2::element_blank(),
    axis.ticks.x = ggplot2::element_blank(),
    legend.position = "none"
  ) +
  ggplot2::ylab("Proportion relative to\n T cell treatment") +
  ggplot2::xlab("Perturbation")

ggplot2::ggsave(file.path(output_path, "/QC_plots/gRNA_distribution_relative_to_1t1.png"), height = 2, width = 2.5, dpi = 300)

saveRDS(CRISPR_sgRNA_proportions, file = file.path(output_path, "/QC_plots/gRNA_distribution_relative_to_1t1.rds"))

# 3.2.2. 1:0.5
tcell_cds <- cds[, colData(cds)$treatment %in% c("1to0.5")]
CRISPR_sgRNA_distribution_summary_df_tcell <- colData(tcell_cds) %>%
  as.data.frame() %>%
  dplyr::group_by(gene) %>%
  dplyr::mutate(CRISPR_count = dplyr::n()) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(total_CRISPR_count = dplyr::n()) %>%
  dplyr::select(gene, CRISPR_count, total_CRISPR_count) %>%
  dplyr::distinct() %>%
  dplyr::mutate(CRISPR_frequency = (CRISPR_count / total_CRISPR_count) * 100)

CRISPR_sgRNA_proportions <- dplyr::inner_join(CRISPR_sgRNA_distribution_summary_df_untreated, CRISPR_sgRNA_distribution_summary_df_tcell, by = "gene")
CRISPR_sgRNA_proportions$relative_proportion <- log(CRISPR_sgRNA_proportions$CRISPR_frequency.x / CRISPR_sgRNA_proportions$CRISPR_frequency.y) %>% scale()
CRISPR_rank <- CRISPR_sgRNA_proportions %>%
  dplyr::arrange(relative_proportion) %>%
  dplyr::pull(gene)

ggplot2::ggplot(CRISPR_sgRNA_proportions, ggplot2::aes(x = factor(gene, levels = CRISPR_rank), y = relative_proportion, label = gene, color = abs(relative_proportion) > 2)) +
  ggplot2::geom_point(size = 1, stroke = 0) +
  ggrepel::geom_label_repel(
    data = subset(CRISPR_sgRNA_proportions, abs(relative_proportion) > 2),
    size = 1,
    segment.size = 0.1,
    label.size = 0.1,
    box.padding = 0,
    label.padding = 0.1,
    color = "black"
  ) +
  monocle3:::monocle_theme_opts() +
  ggplot2::scale_color_manual(values = c("TRUE" = "firebrick3", "FALSE" = "black")) +
  ggplot2::theme(
    text = ggplot2::element_text(size = 6),
    axis.text.x = ggplot2::element_blank(),
    axis.ticks.x = ggplot2::element_blank(),
    legend.position = "none"
  ) +
  ggplot2::ylab("Proportion relative to\n T cell treatment") +
  ggplot2::xlab("Perturbation")

ggplot2::ggsave(file.path(output_path, "/QC_plots/gRNA_distribution_relative_to_1t0.5.png"), height = 2, width = 2.5, dpi = 300)

saveRDS(CRISPR_sgRNA_proportions, file = file.path(output_path, "/QC_plots/gRNA_distribution_relative_to_1t0.5.rds"))

# 3.2.3. 1:0.25
tcell_cds <- cds[, colData(cds)$treatment %in% c("1to0.25")]
CRISPR_sgRNA_distribution_summary_df_tcell <- colData(tcell_cds) %>%
  as.data.frame() %>%
  dplyr::group_by(gene) %>%
  dplyr::mutate(CRISPR_count = dplyr::n()) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(total_CRISPR_count = dplyr::n()) %>%
  dplyr::select(gene, CRISPR_count, total_CRISPR_count) %>%
  dplyr::distinct() %>%
  dplyr::mutate(CRISPR_frequency = (CRISPR_count / total_CRISPR_count) * 100)

CRISPR_sgRNA_proportions <- dplyr::inner_join(CRISPR_sgRNA_distribution_summary_df_untreated, CRISPR_sgRNA_distribution_summary_df_tcell, by = "gene")
CRISPR_sgRNA_proportions$relative_proportion <- log(CRISPR_sgRNA_proportions$CRISPR_frequency.x / CRISPR_sgRNA_proportions$CRISPR_frequency.y) %>% scale()
CRISPR_rank <- CRISPR_sgRNA_proportions %>%
  dplyr::arrange(relative_proportion) %>%
  dplyr::pull(gene)

ggplot2::ggplot(CRISPR_sgRNA_proportions, ggplot2::aes(x = factor(gene, levels = CRISPR_rank), y = relative_proportion, label = gene, color = abs(relative_proportion) > 2)) +
  ggplot2::geom_point(size = 1, stroke = 0) +
  ggrepel::geom_label_repel(
    data = subset(CRISPR_sgRNA_proportions, abs(relative_proportion) > 2),
    size = 1,
    segment.size = 0.1,
    label.size = 0.1,
    box.padding = 0,
    label.padding = 0.1,
    color = "black"
  ) +
  monocle3:::monocle_theme_opts() +
  ggplot2::scale_color_manual(values = c("TRUE" = "firebrick3", "FALSE" = "black")) +
  ggplot2::theme(
    text = ggplot2::element_text(size = 6),
    axis.text.x = ggplot2::element_blank(),
    axis.ticks.x = ggplot2::element_blank(),
    legend.position = "none"
  ) +
  ggplot2::ylab("Proportion relative to\n T cell treatment") +
  ggplot2::xlab("Perturbation")

ggplot2::ggsave(file.path(output_path, "/QC_plots/gRNA_distribution_relative_to_1t0.25.png"), height = 2, width = 2.5, dpi = 300)

saveRDS(CRISPR_sgRNA_proportions, file = file.path(output_path, "/QC_plots/gRNA_distribution_relative_to_1t0.25.rds"))

# 3.3. sgRNA counts
df <- data.frame(colData(cds))
category_counts <- table(df$sgRNA)
df <- data.frame(category = names(category_counts), count = as.numeric(category_counts))

ggplot2::ggplot(df, ggplot2::aes(x = category, y = count)) +
  ggplot2::geom_bar(stat = "identity", fill = "skyblue") +
  ggplot2::labs(x = "Categories", y = "Count", title = "Histogram of Categories") +
  ggplot2::theme_minimal()

# Counts by gene
category_counts <- table(df$category)
df_gene <- data.frame(category = names(category_counts), count = as.numeric(category_counts))

ggplot2::ggplot(df_gene, ggplot2::aes(x = category, y = count)) +
  ggplot2::geom_bar(stat = "identity", fill = "skyblue") +
  ggplot2::labs(x = "Categories", y = "Count", title = "Histogram of Categories") +
  ggplot2::theme_minimal()

sorted_df <- df_gene %>%
  dplyr::arrange(dplyr::desc(count))
print(head(sorted_df))

# Filter out NTC and random
sorted_df <- sorted_df[sorted_df$count <= 1000, ]

ggplot2::ggplot(sorted_df, ggplot2::aes(x = category, y = count)) +
  ggplot2::geom_bar(stat = "identity", fill = "skyblue") +
  ggplot2::labs(x = "Categories", y = "Count", title = "Histogram of Categories") +
  ggplot2::theme_minimal()

# Histogram of counts
ggplot2::ggplot(sorted_df, ggplot2::aes(x = count)) +
  ggplot2::geom_histogram(binwidth = 1, fill = "gray60", color = "black") +
  ggplot2::geom_vline(xintercept = mean(sorted_df$count), linetype = "dashed", color = "red") +
  ggplot2::labs(
    title = "Histogram of sgRNA.sequence counts",
    x = "Count", y = "Frequency"
  ) +
  ggplot2::theme_minimal()

ggplot2::ggsave(file.path(output_path, "/QC_plots/sgRNA_sequence_counts_histogram.png"),
  width = 6, height = 4, dpi = 300
)

# 3.4 Treatment counts
category_counts <- table(colData(hash_cds)$treatment)
df_treat <- data.frame(category = names(category_counts), count = as.numeric(category_counts))

ggplot2::ggplot(df_treat, ggplot2::aes(x = category, y = count)) +
  ggplot2::geom_bar(stat = "identity", fill = "skyblue") +
  ggplot2::labs(x = "Categories", y = "Count", title = "Histogram of Categories") +
  ggplot2::theme_minimal()

ggplot2::ggsave(file.path(output_path, "/QC_plots/sgRNA_sequence_counts_histogram_hash.png"), width = 6, height = 4, dpi = 300)

print(stats::median(colData(cds)[colData(cds)$n.umi >= 500, ]$n.umi))  # 2681

# Create more columns for input
colData(cds)$cell <- row.names(colData(cds))
colData(cds)$RT_lig <- substr(colData(cds)$cell, 9, nchar(colData(cds)$cell))

# 4. Evaluate sgRNA efficiency
# 4.1. Generate size factor and gene statistics
print(length(unique(colData(cds)$gene)))  # Number of guides

cds <- detect_genes(cds)

rowData(cds) %>%
  as.data.frame() %>%
  dplyr::filter(gene_short_name %in% unique(colData(cds)$gene)) %>%
  dplyr::arrange(num_cells_expressed) %>%
  head(n = 150)

print(colData(cds))

colData(cds) %>%
  as.data.frame() %>%
  dplyr::filter(!is.na(gene)) %>%
  dplyr::group_by(gene) %>%
  dplyr::summarise(n = dplyr::n()) %>%
  dplyr::select(n) %>%
  dplyr::pull() %>%
  stats::median()  # 350

print(unique(colData(cds)$gene))

cds <- estimate_size_factors(cds)

# 4.2. Filter for Untreated and NTC
Untreated_cds <- cds[, colData(cds)$treatment %in% c("Untreated")]
# NOTE: targets must be defined somewhere (e.g., from protein_to_gene_map$gene_short_name)
# Here we assume targets already exist in the environment.
targets <- unique(colData(cds)$gene)

filtered_cds <- subset(Untreated_cds, rowData(Untreated_cds)$gene_short_name %in% targets)
NTC_filtered_cds <- filtered_cds[, filtered_cds$gene %in% c("NTC")]
print(NTC_filtered_cds)

# Get NTC gene expression
gene_matrix <- exprs(NTC_filtered_cds)
rownames(gene_matrix) <- rowData(NTC_filtered_cds)$gene_short_name
print(head(rowMeans(gene_matrix)))

# 4.3. Calculate mean knockdown
colData(Untreated_cds)$gene_id <- colData(Untreated_cds)$gene
colData(Untreated_cds)$cell <- colData(Untreated_cds)$cell_ID

Untreated_cds_f <- Untreated_cds[, Untreated_cds$total_sgrna_read_per_cell > 5]
Untreated_cds <- Untreated_cds_f[, Untreated_cds_f$sgRNA_proportion > 0.3]

mean_expres <- calculate_mean_target_knockdown(Untreated_cds, targets, protein_to_gene_map)
qc_path <- file.path(output_folder, "QC_plots/")
if (!dir.exists(qc_path)) dir.create(qc_path, recursive = TRUE)
saveRDS(mean_expres, file.path(qc_path, "mean_target_knockdown.rds"))

print(stats::median(stats::na.omit(mean_expres$mean_knockdown)))            # 0.537
print(stats::median(stats::na.omit(mean_expres$mean_random_knockdown)))     # 0.905
print(stats::median(stats::na.omit(mean_expres$mean_knockdown_10_cutoff)))  # 0.464
print(stats::median(stats::na.omit(mean_expres$mean_random_knockdown_10_cutoff))) # 0.914
print(stats::median(stats::na.omit(mean_expres$mean_knockdown_5_cutoff)))   # 0.499
print(stats::median(stats::na.omit(mean_expres$mean_random_knockdown_5_cutoff)))  # 0.876

print(
  stats::median(stats::na.omit(mean_expres$mean_knockdown)) /
    stats::median(stats::na.omit(mean_expres$mean_random_knockdown))
)  # 0.593

print(
  stats::median(stats::na.omit(mean_expres$mean_knockdown_10_cutoff)) /
    stats::median(stats::na.omit(mean_expres$mean_random_knockdown_10_cutoff))
)  # 0.507

mean_expres_long <- mean_expres %>%
  tidyr::pivot_longer(
    cols = dplyr::starts_with("mean"),
    names_to = "mean_knockdown",
    values_to = "Value"
  )
mean_expres_long$realorrandom <- "real"
mean_expres_long <- mean_expres_long %>%
  dplyr::mutate(realorrandom = ifelse(grepl("mean_random", mean_knockdown), "random", realorrandom))

mean_expres_long <- mean_expres_long %>%
  dplyr::mutate(mean_knockdown = ifelse(grepl("2_cutoff", mean_knockdown), "2_cutoff", mean_knockdown))

mean_expres_long <- mean_expres_long %>%
  dplyr::mutate(mean_knockdown = ifelse(grepl("5_cutoff", mean_knockdown), "5_cutoff", mean_knockdown))

mean_expres_long <- mean_expres_long %>%
  dplyr::mutate(mean_knockdown = ifelse(grepl("0.2gpro", mean_knockdown), "0.2gpro_10_cutff", mean_knockdown))

mean_expres_long <- mean_expres_long %>%
  dplyr::mutate(mean_knockdown = ifelse(grepl("10_cutoff", mean_knockdown), "10_cutoff", mean_knockdown))

mean_expres_long <- mean_expres_long %>%
  dplyr::mutate(mean_knockdown = ifelse(grepl("random", mean_knockdown), "0_cutoff", mean_knockdown))

mean_expres_long <- mean_expres_long %>%
  dplyr::mutate(mean_knockdown = ifelse(grepl("mean", mean_knockdown), "0_cutoff", mean_knockdown))

mean_expres_long$realorrandom <- factor(mean_expres_long$realorrandom, levels = c("real", "random"))
mean_expres_long$mean_knockdown <- factor(mean_expres_long$mean_knockdown, levels = c("0.2gpro_10_cutff", "10_cutoff", "5_cutoff", "2_cutoff", "0_cutoff"))

# Plot median knockdown across perturbations
ggplot2::ggplot(
  mean_expres_long %>% dplyr::filter(!is.na(Value)),
  ggplot2::aes(x = realorrandom, y = Value)
) +
  ggplot2::geom_bar(stat = "summary", fun = "median", position = "dodge", color = "black", size = 0.1) +
  ggplot2::facet_wrap(~mean_knockdown, scales = "fixed", ncol = 3) +
  monocle3:::monocle_theme_opts() +
  ggplot2::theme(text = ggplot2::element_text(size = 6)) +
  ggplot2::ylab("Median knockdown level\nacross perturbations") +
  ggplot2::xlab("sgRNA assignment")

ggplot2::ggsave(file.path(qc_path, "Median_knockdown_levels_by_cutoff_Figure.png"), width = 2.5, height = 2, dpi = 600)

# Focused on 5_cutoff
mean_exprs_cut_off <- mean_expres_long[mean_expres_long$mean_knockdown == "5_cutoff", ]

ggplot2::ggplot(
  mean_exprs_cut_off %>% dplyr::filter(!is.na(Value)),
  ggplot2::aes(x = realorrandom, y = Value)
) +
  ggplot2::geom_bar(stat = "summary", fun = "median", position = "dodge", color = "black", size = 0.1) +
  monocle3:::monocle_theme_opts() +
  ggplot2::theme(text = ggplot2::element_text(size = 12)) +
  ggplot2::ylab("Median knockdown level \nacross perturbations") +
  ggplot2::xlab("sgRNA assignment")

ggplot2::ggsave(file.path(qc_path, "Median_knockdown_levels_by_0.3g_5_cutoff_Figure.png"), width = 2.5, height = 2, dpi = 600)

# 5. Save the CDS
saveRDS(cds, file.path(base_dir, "/11-final-output/CDS_crop_hash_FINAL_CRISPR_GBM_T_cells_500_cutoff_after_qc.rds"))
