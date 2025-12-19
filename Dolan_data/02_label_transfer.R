
library(Seurat)
library(qs)
library(dplyr)


reference <- qread("/home/garga7/Anjali_data/int_data/SRT_mapping/input_file/Allfour_SCT_integration_layer_wdedited_cells_newanno_SCT.qs", nthreads = 6 )

# to human only
#reference <- subset(reference, subset = batch %in% c("KD","ROSMAP_MIT"))

query <- qread("/home/garga7/Anjali_data/int_data/Dolan_data/iMGLs_depthnorm_liger_seurat_2024_recluster.qs")

features <- query@assays$SCT@var.features[query@assays$SCT@var.features %in% reference@assays$SCT@var.features]

anchors <- FindTransferAnchors(reference = reference, query = query, dims = 1:30,recompute.residuals =FALSE,scale = FALSE
                              ,normalization.method = "SCT", features = features)

predictions <- TransferData(anchorset = anchors, refdata = reference$Integrated_cell_type, dims = 1:30)
qsave(predictions,"/home/garga7/Anjali_data/int_data/Dolan_data/iMGLs_depthnorm_liger_seurat_2024_seurat_predicted_model.qs")
write.csv(predictions,"/home/garga7/Anjali_data/int_data/Dolan_data/iMGLs_depthnorm_liger_seurat_2024_seurat_predicted_model.csv")

query <- AddMetaData(query, metadata = predictions)

table(query$predicted.id, query$Cell_type)

#qsave(query,"/home/garga7/Anjali_data/int_data/Dolan_data/iMGLs_depthnorm_liger_seurat_2024_seurat_predicted_with_humanonly.qs")


# feature 
# reduction= cca 



