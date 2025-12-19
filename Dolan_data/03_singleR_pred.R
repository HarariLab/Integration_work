#Jai Babe Di
#Jai Guru Maa Ji

library(Seurat)
library(Matrix)
library(dplyr)
library(readr)
library(SingleR)
library(qs)
library(scuttle)
library(SingleCellExperiment)


seu_sc <- qread("/home/garga7/Anjali_data/int_data/Dolan_data/iMGLs_depthnorm_liger_seurat_2024_recluster_sce.qs")

obj_sc <- readRDS("~/Anjali_data/int_data/Allfour_SCT_integration_layer_wdedited_meta_sce.rds")

pred <- SingleR(test = seu_sc, ref = obj_sc, labels = obj_sc$Integrated_cell_type)

saveRDS(pred,"/home/garga7/Anjali_data/int_data/Dolan_data/iMGLs_depthnorm_liger_seurat_2024_recluster_singleR_pred.rds")