# ==============================================================================
# NicheNet Analysis Pipeline: GBM-T Cell Communication
# ==============================================================================

# Load required libraries
rm(list = ls())
gc()

library(SingleCellExperiment)
library(dplyr)
library(ggplot2)
library(nichenetr)
library(multinichenetr)
library(SummarizedExperiment)
library(monocle3)
library(tidyr)
library(scater)
library(patchwork)
library(igraph)
library(ggraph)

# Set options
options(stringsAsFactors = FALSE)
options(future.globals.maxSize = 8000 * 1024^2)
options(timeout = 360)
future::plan("sequential")

# Output directory for all figures and tables
output_dir <- "/home/user/Documents/Kinase_project/figures/GBM_Tcell_communication/"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# ==============================================================================
# Load NicheNet's ligand-receptor network and ligand-target matrix
# ==============================================================================

organism = "human"

if(organism == "human"){
  
  lr_network_all = 
    readRDS(url(
      "https://zenodo.org/record/10229222/files/lr_network_human_allInfo_30112033.rds"
    )) %>% 
    mutate(
      ligand = convert_alias_to_symbols(ligand, organism = organism), 
      receptor = convert_alias_to_symbols(receptor, organism = organism))
  
  lr_network_all = lr_network_all  %>% 
    mutate(ligand = make.names(ligand), receptor = make.names(receptor)) 
  
  lr_network = lr_network_all %>% 
    distinct(ligand, receptor)
  
  ligand_target_matrix = readRDS(url(
    "https://zenodo.org/record/7074291/files/ligand_target_matrix_nsga2r_final.rds"
  ))
  
  colnames(ligand_target_matrix) = colnames(ligand_target_matrix) %>% 
    convert_alias_to_symbols(organism = organism) %>% make.names()
  rownames(ligand_target_matrix) = rownames(ligand_target_matrix) %>% 
    convert_alias_to_symbols(organism = organism) %>% make.names()
  
  lr_network = lr_network %>% filter(ligand %in% colnames(ligand_target_matrix))
  ligand_target_matrix = ligand_target_matrix[, lr_network$ligand %>% unique()]
  
} else if(organism == "mouse"){
  
  lr_network_all = readRDS(url(
    "https://zenodo.org/record/10229222/files/lr_network_mouse_allInfo_30112033.rds"
  )) %>% 
    mutate(
      ligand = convert_alias_to_symbols(ligand, organism = organism), 
      receptor = convert_alias_to_symbols(receptor, organism = organism))
  
  lr_network_all = lr_network_all  %>% 
    mutate(ligand = make.names(ligand), receptor = make.names(receptor)) 
  lr_network = lr_network_all %>% 
    distinct(ligand, receptor)
  
  ligand_target_matrix = readRDS(url(
    "https://zenodo.org/record/7074291/files/ligand_target_matrix_nsga2r_final_mouse.rds"
  ))
  
  colnames(ligand_target_matrix) = colnames(ligand_target_matrix) %>% 
    convert_alias_to_symbols(organism = organism) %>% make.names()
  rownames(ligand_target_matrix) = rownames(ligand_target_matrix) %>% 
    convert_alias_to_symbols(organism = organism) %>% make.names()
  
  lr_network = lr_network %>% filter(ligand %in% colnames(ligand_target_matrix))
  ligand_target_matrix = ligand_target_matrix[, lr_network$ligand %>% unique()]
}

# ==============================================================================
# Read in and prepare SingleCellExperiment objects
# ==============================================================================

cds_GBM <- readRDS('/home/user/Documents/Kinase_project/crispra_final/cds_2000.rds')
cds_Tcells <- readRDS('/home/user/Documents/Kinase_project_backup/GBM_T_Cells/GBM_Tcell_joint_cds.rds')

colData(cds_GBM)$celltype <- "GBM"
colData(cds_GBM)$batch <- "GBM"
colData(cds_GBM)$sample_id <- "GBM"
sce_GBM <- as(cds_GBM, "SingleCellExperiment")
reducedDims(cds_GBM) <- list()
reducedDims(cds_Tcells) <- list()

colData(cds_Tcells)$dose <- colData(cds_Tcells)$Tcell_to_GBM_ratio
colData(cds_Tcells)$celltype <- "Tcells"
colData(cds_Tcells)$batch <- "Tcells"
colData(cds_Tcells)$sample_id <- "Tcells"

sce_Tcells <- as(cds_Tcells, "SingleCellExperiment")

sce_GBM = alias_to_symbol_SCE(sce_GBM, "human") %>% makenames_SCE()
sce_Tcells = alias_to_symbol_SCE(sce_Tcells, "human") %>% makenames_SCE()

# Fix rowData - remove conflicting 'id' column
rowData(sce_GBM)$id <- NULL
rowData(sce_Tcells)$id <- NULL

# Fix colData - add missing columns
all_cols <- union(colnames(colData(sce_GBM)), colnames(colData(sce_Tcells)))

for(col in setdiff(colnames(colData(sce_Tcells)), colnames(colData(sce_GBM)))) {
  colData(sce_GBM)[[col]] <- NA
}

for(col in setdiff(colnames(colData(sce_GBM)), colnames(colData(sce_Tcells)))) {
  colData(sce_Tcells)[[col]] <- NA
}

rownames(sce_GBM) <- rowData(sce_GBM)$gene_short_name
rownames(sce_Tcells) <- rowData(sce_Tcells)$gene_short_name

rowData(sce_GBM)$gene_short_name <- NULL
rowData(sce_Tcells)$gene_short_name <- NULL

colData(sce_GBM) <- colData(sce_GBM)[, all_cols]
colData(sce_Tcells) <- colData(sce_Tcells)[, all_cols]

# Find common genes
common_genes <- intersect(rownames(sce_GBM), rownames(sce_Tcells))

print(paste("GBM genes:", nrow(sce_GBM)))
print(paste("Tcells genes:", nrow(sce_Tcells)))
print(paste("Common genes:", length(common_genes)))

# Subset to common genes
sce_GBM_common <- sce_GBM[common_genes, ]
sce_Tcells_common <- sce_Tcells[common_genes, ]

reducedDims(sce_GBM_common) <- NULL
reducedDims(sce_Tcells_common) <- NULL

# Combine objects
sce <- cbind(sce_GBM_common, sce_Tcells_common)

print("Unique doses:")
print(unique(colData(sce)$dose))
print("Unique cell types:")
print(unique(colData(sce)$celltype))

# ==============================================================================
# Balanced Cell Sampling
# ==============================================================================

doses_to_compare <- c(0.00, 1.00)
sce_filtered <- sce[, colData(sce)$dose %in% doses_to_compare]

print("Original cell distribution:")
print(table(colData(sce_filtered)$dose, colData(sce_filtered)$celltype))

cells_to_keep <- c()

for(dose_val in doses_to_compare) {
  tcell_indices <- which(colData(sce_filtered)$dose == dose_val & 
                        colData(sce_filtered)$celltype == "Tcells")
  n_tcells <- length(tcell_indices)
  
  gbm_indices <- which(colData(sce_filtered)$dose == dose_val & 
                      colData(sce_filtered)$celltype == "GBM")
  
  set.seed(123)
  sampled_gbm_indices <- sample(gbm_indices, n_tcells * 2)
  
  dose_cells <- c(tcell_indices, sampled_gbm_indices)
  cells_to_keep <- c(cells_to_keep, dose_cells)
  
  print(paste("Dose", dose_val, "- T cells:", n_tcells, ", GBM:", n_tcells * 2, 
              "(sampled from", length(gbm_indices), ")"))
}

sce_balanced <- sce_filtered[, cells_to_keep]

print("Balanced cell distribution:")
balanced_counts <- table(colData(sce_balanced)$dose, colData(sce_balanced)$celltype)
print(balanced_counts)

# Filter genes
print("Filtering genes...")
data_matrix <- as.matrix(counts(sce_balanced))
gene_detection_rate <- rowMeans(data_matrix > 0)

keep_genes <- gene_detection_rate > 0.03
data_matrix_filtered <- data_matrix[keep_genes, ]
sce_balanced <- sce_balanced[keep_genes, ]

print(paste("Kept", sum(keep_genes), "out of", length(keep_genes), "genes"))
print(paste("Final dataset:", nrow(sce_balanced), "genes x", ncol(sce_balanced), "cells"))

# ==============================================================================
# Analysis 1: GBM → T cells
# ==============================================================================

print("\n=== ANALYSIS 1: GBM → T CELLS ===\n")

print("Cell distribution:")
print(table(colData(sce_balanced)$dose, colData(sce_balanced)$celltype))

metadata <- as.data.frame(colData(sce_balanced))
receiver = "Tcells"
sender = "GBM" 
condition_oi = "1"
condition_reference = "0"

sce_balanced <- logNormCounts(sce_balanced)
expression <- as.matrix(logcounts(sce_balanced))

# Define expressed genes
expressed_genes_sender <- rownames(expression)[
  rowMeans(expression[, metadata$celltype == sender] > 0) >= 0.10
]
expressed_genes_receiver <- rownames(expression)[
  rowMeans(expression[, metadata$celltype == receiver] > 0) >= 0.10
]

print(paste("Expressed genes in", sender, ":", length(expressed_genes_sender)))
print(paste("Expressed genes in", receiver, ":", length(expressed_genes_receiver)))

# DE analysis
receiver_cells_oi <- metadata$celltype == receiver & metadata$dose == condition_oi
receiver_cells_ref <- metadata$celltype == receiver & metadata$dose == condition_reference

DE_genes <- sapply(rownames(expression), function(gene) {
  expr_oi <- expression[gene, receiver_cells_oi]
  expr_ref <- expression[gene, receiver_cells_ref]
  
  if(mean(c(expr_oi, expr_ref)) == 0) return(c(pval = 1, logFC = 0))
  
  test <- wilcox.test(expr_oi, expr_ref)
  logFC <- mean(expr_oi) - mean(expr_ref)
  
  c(pval = test$p.value, logFC = logFC)
})

DE_table <- data.frame(
  gene = colnames(DE_genes),
  p_val = DE_genes["pval",],
  logFC = DE_genes["logFC",],
  stringsAsFactors = FALSE
)
DE_table$p_val_adj <- p.adjust(DE_table$p_val, method = "BH")

geneset_oi <- DE_table %>%
  filter(logFC > 0.25 & p_val_adj <= 0.05) %>%
  pull(gene) %>%
  intersect(rownames(ligand_target_matrix))

if(length(geneset_oi) < 20) {
  geneset_oi <- DE_table %>%
    filter(logFC > 0) %>%
    arrange(desc(logFC)) %>%
    slice_head(n = 100) %>%
    pull(gene) %>%
    intersect(rownames(ligand_target_matrix))
}

print(paste("Gene set of interest size:", length(geneset_oi)))

background_expressed_genes <- expressed_genes_receiver %>%
  intersect(rownames(ligand_target_matrix))

# Define potential ligands
ligands <- lr_network %>% pull(ligand) %>% unique()
receptors <- lr_network %>% pull(receptor) %>% unique()

expressed_ligands <- intersect(ligands, expressed_genes_sender)
expressed_receptors <- intersect(receptors, expressed_genes_receiver)

print(paste("Expressed ligands:", length(expressed_ligands)))
print(paste("Expressed receptors:", length(expressed_receptors)))

potential_ligands <- lr_network %>%
  filter(ligand %in% expressed_ligands & receptor %in% expressed_receptors) %>%
  pull(ligand) %>%
  unique()

print(paste("Potential ligands:", length(potential_ligands)))

# Ligand activity analysis
ligand_activities <- predict_ligand_activities(
  geneset = geneset_oi,
  background_expressed_genes = background_expressed_genes,
  ligand_target_matrix = ligand_target_matrix,
  potential_ligands = potential_ligands
)

ligand_activities <- ligand_activities %>%
  arrange(-aupr_corrected) %>%
  mutate(rank = rank(desc(aupr_corrected)))

print("Top ligands by activity:")
print(head(ligand_activities, 20))

best_upstream_ligands <- ligand_activities %>%
  top_n(20, aupr_corrected) %>%
  arrange(-aupr_corrected) %>%
  pull(test_ligand)

# Infer ligand-target links
active_ligand_target_links_df <- best_upstream_ligands %>%
  lapply(get_weighted_ligand_target_links,
         geneset = geneset_oi,
         ligand_target_matrix = ligand_target_matrix,
         n = 200) %>%
  bind_rows() %>%
  drop_na()

lr_network_top <- lr_network %>%
  filter(ligand %in% best_upstream_ligands & receptor %in% expressed_receptors) %>%
  select(ligand, receptor)

# Calculate expression
ligand_expression <- sapply(best_upstream_ligands, function(lig) {
  if(lig %in% rownames(expression)) {
    mean(expression[lig, metadata$celltype == sender & metadata$dose == condition_oi])
  } else NA
})

receptor_list <- lr_network_top %>% pull(receptor) %>% unique()
receptor_expression <- sapply(receptor_list, function(rec) {
  if(rec %in% rownames(expression)) {
    mean(expression[rec, metadata$celltype == receiver & metadata$dose == condition_oi])
  } else NA
})

ligand_receptor_df <- lr_network_top %>%
  left_join(ligand_activities %>% select(test_ligand, aupr_corrected, rank),
            by = c("ligand" = "test_ligand")) %>%
  arrange(-aupr_corrected)

ligand_receptor_df$ligand_expression <- ligand_expression[ligand_receptor_df$ligand]
ligand_receptor_df$receptor_expression <- receptor_expression[ligand_receptor_df$receptor]

print("=== TOP LIGAND-RECEPTOR PAIRS ===")
print(head(ligand_receptor_df, 20))

# Save results
write.csv(ligand_receptor_df,
          file.path(output_dir, "nichenet_results_corrected_GBM_Tcells.csv"),
          row.names = FALSE)
write.csv(ligand_activities,
          file.path(output_dir, "ligand_activities_corrected_GBM_Tcells.csv"),
          row.names = FALSE)
write.csv(active_ligand_target_links_df,
          file.path(output_dir, "ligand_target_links_corrected_GBM_Tcells.csv"),
          row.names = FALSE)

print("=== ANALYSIS COMPLETE ===")
print(paste("- Identified", length(geneset_oi), "genes of interest"))
print(paste("- Tested", length(potential_ligands), "potential ligands"))
print(paste("- Top", length(best_upstream_ligands), "active ligands selected"))
print(paste("- Found", nrow(ligand_receptor_df), "ligand-receptor pairs"))

# ==============================================================================
# Analysis 2: T cells → GBM
# ==============================================================================

print("\n=== ANALYSIS 2: T CELLS → GBM ===\n")

receiver = "GBM"
sender = "Tcells"

expressed_genes_sender <- rownames(expression)[
  rowMeans(expression[, metadata$celltype == sender] > 0) >= 0.10
]
expressed_genes_receiver <- rownames(expression)[
  rowMeans(expression[, metadata$celltype == receiver] > 0) >= 0.10
]

print(paste("Expressed genes in", sender, ":", length(expressed_genes_sender)))
print(paste("Expressed genes in", receiver, ":", length(expressed_genes_receiver)))

receiver_cells_oi <- metadata$celltype == receiver & metadata$dose == condition_oi
receiver_cells_ref <- metadata$celltype == receiver & metadata$dose == condition_reference

DE_genes <- sapply(rownames(expression), function(gene) {
  expr_oi <- expression[gene, receiver_cells_oi]
  expr_ref <- expression[gene, receiver_cells_ref]
  
  if(mean(c(expr_oi, expr_ref)) == 0) return(c(pval = 1, logFC = 0))
  
  test <- wilcox.test(expr_oi, expr_ref)
  logFC <- mean(expr_oi) - mean(expr_ref)
  
  c(pval = test$p.value, logFC = logFC)
})

DE_table <- data.frame(
  gene = colnames(DE_genes),
  p_val = DE_genes["pval",],
  logFC = DE_genes["logFC",],
  stringsAsFactors = FALSE
)
DE_table$p_val_adj <- p.adjust(DE_table$p_val, method = "BH")

geneset_oi <- DE_table %>%
  filter(logFC > 0.25 & p_val_adj <= 0.05) %>%
  pull(gene) %>%
  intersect(rownames(ligand_target_matrix))

if(length(geneset_oi) < 20) {
  geneset_oi <- DE_table %>%
    filter(logFC > 0) %>%
    arrange(desc(logFC)) %>%
    slice_head(n = 100) %>%
    pull(gene) %>%
    intersect(rownames(ligand_target_matrix))
}

print(paste("Gene set of interest size:", length(geneset_oi)))

background_expressed_genes <- expressed_genes_receiver %>%
  intersect(rownames(ligand_target_matrix))

expressed_ligands <- intersect(ligands, expressed_genes_sender)
expressed_receptors <- intersect(receptors, expressed_genes_receiver)

print(paste("Expressed ligands:", length(expressed_ligands)))
print(paste("Expressed receptors:", length(expressed_receptors)))

potential_ligands <- lr_network %>%
  filter(ligand %in% expressed_ligands & receptor %in% expressed_receptors) %>%
  pull(ligand) %>%
  unique()

print(paste("Potential ligands:", length(potential_ligands)))

ligand_activities <- predict_ligand_activities(
  geneset = geneset_oi,
  background_expressed_genes = background_expressed_genes,
  ligand_target_matrix = ligand_target_matrix,
  potential_ligands = potential_ligands
)

ligand_activities <- ligand_activities %>%
  arrange(-aupr_corrected) %>%
  mutate(rank = rank(desc(aupr_corrected)))

print("Top ligands by activity:")
print(head(ligand_activities, 20))

best_upstream_ligands <- ligand_activities %>%
  top_n(20, aupr_corrected) %>%
  arrange(-aupr_corrected) %>%
  pull(test_ligand)

active_ligand_target_links_df <- best_upstream_ligands %>%
  lapply(get_weighted_ligand_target_links,
         geneset = geneset_oi,
         ligand_target_matrix = ligand_target_matrix,
         n = 200) %>%
  bind_rows() %>%
  drop_na()

lr_network_top <- lr_network %>%
  filter(ligand %in% best_upstream_ligands & receptor %in% expressed_receptors) %>%
  select(ligand, receptor)

ligand_expression <- sapply(best_upstream_ligands, function(lig) {
  if(lig %in% rownames(expression)) {
    mean(expression[lig, metadata$celltype == sender & metadata$dose == condition_oi])
  } else NA
})

receptor_list <- lr_network_top %>% pull(receptor) %>% unique()
receptor_expression <- sapply(receptor_list, function(rec) {
  if(rec %in% rownames(expression)) {
    mean(expression[rec, metadata$celltype == receiver & metadata$dose == condition_oi])
  } else NA
})

ligand_receptor_df <- lr_network_top %>%
  left_join(ligand_activities %>% select(test_ligand, aupr_corrected, rank),
            by = c("ligand" = "test_ligand")) %>%
  arrange(-aupr_corrected)

ligand_receptor_df$ligand_expression <- ligand_expression[ligand_receptor_df$ligand]
ligand_receptor_df$receptor_expression <- receptor_expression[ligand_receptor_df$receptor]

print("=== TOP LIGAND-RECEPTOR PAIRS ===")
print(head(ligand_receptor_df, 40))

write.csv(ligand_receptor_df,
          file.path(output_dir, "nichenet_results_corrected_Tcells_GBM.csv"),
          row.names = FALSE)
write.csv(ligand_activities,
          file.path(output_dir, "ligand_activities_corrected_Tcells_GBM.csv"),
          row.names = FALSE)
write.csv(active_ligand_target_links_df,
          file.path(output_dir, "ligand_target_links_corrected_Tcells_GBM.csv"),
          row.names = FALSE)

print("=== ANALYSIS COMPLETE ===")
print(paste("- Identified", length(geneset_oi), "genes of interest"))
print(paste("- Tested", length(potential_ligands), "potential ligands"))
print(paste("- Top", length(best_upstream_ligands), "active ligands selected"))
print(paste("- Found", nrow(ligand_receptor_df), "ligand-receptor pairs"))

# ==============================================================================
# Network Visualization Function
# ==============================================================================

create_network_plot <- function(ligand_receptor_df, top_n = 15, 
                                layout_type = "circle", 
                                save_plot = TRUE,
                                filename = "nichenet_network.png",
                                sender_name = "Sender",
                                receiver_name = "Receiver") {
  
  print(paste("Creating network plot with top", top_n, "interactions..."))
  
  network_data <- ligand_receptor_df %>%
    arrange(desc(aupr_corrected)) %>%
    slice_head(n = top_n) %>%
    select(ligand, receptor, aupr_corrected) %>%
    rename(from = ligand, to = receptor, weight = aupr_corrected) %>%
    filter(!is.na(weight))
  
  print(paste("Network has", nrow(network_data), "edges"))
  
  g <- graph_from_data_frame(network_data, directed = TRUE)
  
  V(g)$type <- ifelse(V(g)$name %in% network_data$from, "Ligand", "Receptor")
  V(g)$node_color <- ifelse(V(g)$type == "Ligand", "#87CEEB", "#FFA500")
  V(g)$node_size <- ifelse(V(g)$type == "Ligand", 8, 6)
  
  E(g)$edge_width <- scales::rescale(E(g)$weight, to = c(0.5, 3))
  E(g)$edge_alpha <- scales::rescale(E(g)$weight, to = c(0.3, 0.8))
  
  print(paste("Graph has", vcount(g), "nodes and", ecount(g), "edges"))
  
  # Create dynamic labels
  ligand_label <- paste0("Ligand (", sender_name, ")")
  receptor_label <- paste0("Receptor (", receiver_name, ")")
  
  p_network <- ggraph(g, layout = layout_type) +
    geom_edge_arc(aes(width = edge_width, alpha = edge_alpha), 
                  arrow = arrow(length = unit(4, 'mm'), type = "closed"),
                  start_cap = circle(4, 'mm'),
                  end_cap = circle(4, 'mm'),
                  color = "gray40",
                  strength = 0.1) +
    geom_node_point(aes(color = type, size = type), alpha = 0.8) +
    geom_node_text(aes(label = name), 
                   size = 3.5, 
                   repel = TRUE,
                   point.padding = unit(0.3, "lines"),
                   fontface = "bold") +
    scale_color_manual(name = "Cell Type",
                      values = c("Ligand" = "#4682B4", "Receptor" = "#FF8C00"),
                      labels = c("Ligand" = ligand_label, "Receptor" = receptor_label)) +
    scale_size_manual(name = "Cell Type",
                     values = c("Ligand" = 6, "Receptor" = 4),
                     labels = c("Ligand" = ligand_label, "Receptor" = receptor_label)) +
    scale_edge_width_identity() +
    scale_edge_alpha_identity() +
    theme_void() +
    labs(title = "Ligand-Receptor Communication Network",
         subtitle = paste("Top", top_n, sender_name, "→", receiver_name, 
                         "interactions ranked by ligand activity (AUPR)"),
         caption = "Edge thickness = Ligand activity score\nNode size indicates cell type") +
    theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
          plot.subtitle = element_text(size = 12, hjust = 0.5),
          plot.caption = element_text(size = 10, hjust = 0.5),
          legend.position = "bottom",
          legend.title = element_text(size = 12, face = "bold"),
          legend.text = element_text(size = 10),
          plot.margin = margin(20, 20, 20, 20))
  
  print(p_network)
  
  if(save_plot) {
    ggsave(file.path(output_dir, filename), p_network, 
           width = 8, height = 7, dpi = 300, bg = "white")
    print(paste("Network plot saved as:", filename))
  }
  
  return(p_network)
}

# ==============================================================================
# Create Network Plots for T cells → GBM
# ==============================================================================

print("\n=== CREATING NICHENET NETWORK PLOTS (T CELLS → GBM) ===\n")

ligand_receptor_df <- read.csv(file.path(output_dir, "nichenet_results_corrected_Tcells_GBM.csv"))

print(paste("Total ligand-receptor pairs:", nrow(ligand_receptor_df)))
print(paste("Top AUPR score:", round(max(ligand_receptor_df$aupr_corrected, na.rm = TRUE), 3)))

p_circle <- create_network_plot(ligand_receptor_df, 
                               top_n = 15, 
                               layout_type = "circle",
                               filename = "T_GBM_network_circle.png",
                               sender_name = "T cells",
                               receiver_name = "GBM")

p_stress <- create_network_plot(ligand_receptor_df, 
                               top_n = 15, 
                               layout_type = "stress",
                               filename = "T_GBM_network_stress.png",
                               sender_name = "T cells",
                               receiver_name = "GBM")

p_large <- create_network_plot(ligand_receptor_df, 
                              top_n = 25, 
                              layout_type = "fr",
                              filename = "T_GBM_network_large.png",
                              sender_name = "T cells",
                              receiver_name = "GBM")

print("\n=== NETWORK VISUALIZATION COMPLETE ===")
print(paste("All files saved in:", output_dir))

top_interactions <- ligand_receptor_df %>%
  arrange(desc(aupr_corrected)) %>%
  slice_head(n = 10) %>%
  select(ligand, receptor, aupr_corrected)

print("Top 10 ligand-receptor interactions:")
print(top_interactions)

# ==============================================================================
# Create Network Plots for GBM → T cells
# ==============================================================================

print("\n=== CREATING NICHENET NETWORK PLOTS (GBM → T CELLS) ===\n")

ligand_receptor_df <- read.csv(file.path(output_dir, "nichenet_results_corrected_GBM_Tcells.csv"))

print(paste("Total ligand-receptor pairs:", nrow(ligand_receptor_df)))
print(paste("Top AUPR score:", round(max(ligand_receptor_df$aupr_corrected, na.rm = TRUE), 3)))

p_circle <- create_network_plot(ligand_receptor_df, 
                               top_n = 15, 
                               layout_type = "circle",
                               filename = "GBM_Tcells_network_circle.png",
                               sender_name = "GBM",
                               receiver_name = "T cells")

p_stress <- create_network_plot(ligand_receptor_df, 
                               top_n = 15, 
                               layout_type = "stress",
                               filename = "GBM_Tcells_network_stress.png",
                               sender_name = "GBM",
                               receiver_name = "T cells")

p_large <- create_network_plot(ligand_receptor_df, 
                              top_n = 25, 
                              layout_type = "fr",
                              filename = "GBM_Tcells_network_large.png",
                              sender_name = "GBM",
                              receiver_name = "T cells")

print("\n=== ALL ANALYSES COMPLETE ===")
print(paste("All results and plots saved in:", output_dir))

top_interactions <- ligand_receptor_df %>%
  arrange(desc(aupr_corrected)) %>%
  slice_head(n = 10) %>%
  select(ligand, receptor, aupr_corrected)

print("Top 10 ligand-receptor interactions:")
print(top_interactions)