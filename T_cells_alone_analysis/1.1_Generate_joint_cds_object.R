library(devtools)
library(ggplot2)
library(dplyr)
library(tidyr)
library(monocle3)

source("/home/user/Documents/github_repos/sci-plex/bin/cell_cycle.R")
source("/home/user/Documents/github_repos/sci-plex/bin/dispersions_functions.R")
cc.genes <- readRDS("/home/user/Documents/github_repos/sci-plex/bin/cc.genes.RDS")

# Set DelayedArray Parameters
DelayedArray:::set_verbose_block_processing(TRUE)
options(DelayedArray.block.size=1000e7)
options(stringsAsFactors = FALSE)
###########

setwd("/home/user/Documents/Kinase_project_backup/GBM_T_Cells/")

untreated_cds <- load_cellranger_data(pipestance_path = "/home/user/Documents/Kinase_project_backup/GBM_T_Cells/Gbm_tcell_crispra_un_8/")
one_point_five_ratio_cds <- load_cellranger_data(pipestance_path ="/home/user/Documents/Kinase_project_backup/GBM_T_Cells/Gbm_tcell_crispra_one_point_five/")
one_one_ratio_cds <- load_cellranger_data(pipestance_path = "/home/user/Documents/Kinase_project_backup/GBM_T_Cells/Gbm_tcell_crispra_one_one/")
one_two_ratio_cds <- load_cellranger_data(pipestance_path = "/home/user/Documents/Kinase_project_backup/GBM_T_Cells/Gbm_tcell_crispra_one_two_6/")

GBM_Tcell_cds <- combine_cds(list(untreated_cds,one_point_five_ratio_cds, one_one_ratio_cds, one_two_ratio_cds))
colData(GBM_Tcell_cds)$cell <- row.names(colData(GBM_Tcell_cds))
colData(GBM_Tcell_cds)$sample <- sapply(colData(GBM_Tcell_cds)$cell,function(x){
  sample_id <- stringr::str_split(x, pattern = "_")[[1]][2]
  if(sample_id == "1")return("untreated_GBM")
  if(sample_id == "2")return("1to05_GBM_Tcell")
  if(sample_id == "3")return("1to1_GBM_Tcell")
  if(sample_id == "4")return("1to2_GBM_Tcell")
  return(NA)})
colData(GBM_Tcell_cds)$Tcell_to_GBM_ratio <- sapply(colData(GBM_Tcell_cds)$cell,function(x){
  sample_id <- stringr::str_split(x, pattern = "_")[[1]][2]
  if(sample_id == "1")return(0)
  if(sample_id == "2")return(0.5)
  if(sample_id == "3")return(1)
  if(sample_id == "4")return(2)
  return(NA)})

saveRDS(GBM_Tcell_cds, "GBM_Tcell_joint_cds.rds")

