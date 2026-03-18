# use the conda env rds2H5ad
# Convert rds to h5ad
library(reticulate)

# # Install devtools if needed
# if (!requireNamespace("devtools", quietly = TRUE)) {
#     install.packages("devtools")
# }

# # Install SeuratDisk from GitHub
# devtools::install_github("mojaveazure/seurat-disk")

# Then load it
library(SeuratDisk)

# Need to uninstall Matrix package v1.6-1 and reinstall v1.6-4.
library(Seurat)
library(SeuratObject)
library(dplyr)
library(monocle3)

setwd('/home/user/Documents/Kinase_project/20250430_BT333_Tcell_11chems_run2/')
# Change filename
filename <- "preprocessed_cds_with_gene_sig"
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
sciNM_gene_short_names <- cds.pre[rowData(cds.pre)$id %in% unique_ids, ]
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

