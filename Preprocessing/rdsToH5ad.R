# Convert rds to h5ad 

library(reticulate)
library(Seurat)
library(SeuratObject)
library(SeuratDisk)
library(dplyr)
library(monocle3)

setwd("/burg/iicd/users/qc2358/Kinase_project")

# Change filename
filename <- "11-final-output/CDS_crop_hash_FINAL_CRISPR_GBM_T_cells_500_cutoff" 
cds.pre <- readRDS(paste0(filename, ".rds"))

# Get rownames of rowData to be gene short names for seurat obj
rowDataDF <- rowData(cds.pre)

# Identify the first unique id for each unique gene_short_name
unique_gene_short_names <- unique(rowDataDF$gene_short_name)
unique_ids <- character(length(unique_gene_short_names))
for (i in seq_along(unique_gene_short_names)) {
  gene_short_name <- unique_gene_short_names[i]
  first_unique_id <- rowDataDF$id[rowDataDF$gene_short_name == gene_short_name][1]
  unique_ids[i] <- first_unique_id
}

# Generate exprs matrix
sciNM_gene_short_names <- cds.pre[rowData(cds.pre)$id %in% unique_ids,]
exprs_sciNM_gene_short_names <- exprs(sciNM_gene_short_names)
rownames(rowData(sciNM_gene_short_names)) <- rowData(sciNM_gene_short_names)$gene_short_name
rownames(exprs_sciNM_gene_short_names) <- rowData(sciNM_gene_short_names)$gene_short_name
colnames(exprs_sciNM_gene_short_names) <- colData(cds.pre)$cell_ID

# Make seurat object then make H5 seurat then make an anndata from that
so_for_spectra <- CreateSeuratObject(exprs_sciNM_gene_short_names, meta.data = colData(sciNM_gene_short_names) %>% as.data.frame(), assay = "RNA")
so_for_spectra[["RNA"]] <- CreateAssayObject(counts = exprs_sciNM_gene_short_names)
so_for_spectra@meta.data <- colData(sciNM_gene_short_names) %>% as.data.frame()

SaveH5Seurat(so_for_spectra, filename = paste0(filename, ".h5Seurat"))
Convert(source = paste0(filename, ".h5Seurat"), dest = "h5ad", filename = paste0(filename, ".h5ad"))
