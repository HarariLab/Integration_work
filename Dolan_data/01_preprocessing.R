# Jai Babe Di
# Jai Guru Maa Ji


library(Seurat)
library(qs)
library(dplyr)

obj <- readRDS("/home/garga7/Anjali_data/int_data/Dolan_data/iMGLs_depthnorm_liger_seurat_2024.rds")


obj$Cell_type <- obj$liger_cluster
# Step 1: Define mapping vector
mapping <- c(
  "1" = "Transition",
  "2" = "DAM",
  "3" = "Antigen-presenting",
  "4" = "Antigen-presenting",
  "5" = "Homeostatic",
  "6" = "Proliferation",
  "7" = "Antigen-presenting",
  "8" = "DAM",
  "9" = "Proliferation",
  "10" = "Proliferation",
  "11" = "Interferon-responsive"
)

obj@meta.data$Cell_type <- as.character(obj@meta.data$Cell_type)
obj@meta.data$Cell_type <- mapping[obj@meta.data$Cell_type]

DefaultAssay(obj) <- "RNA"
obj@meta.data$batch <- "Dolan"

obj_mdata <- obj@meta.data %>%
  dplyr::select(
    orig.ident,
    liger_cluster,
    batch,
    Cell_type
  )

obj <- CreateSeuratObject(
  counts = obj@assays$RNA$counts, meta.data = obj@meta.data
)

obj <- PercentageFeatureSet(obj, "^MT-", col.name = "percent.mt")
obj <- PercentageFeatureSet(obj, "^RPS", col.name = "percent.rbs")


sobj_new <- SCTransform(
  obj,
  vars.to.regress = c('nCount_RNA','nFeature_RNA','percent.mt'),
  verbose=TRUE
)

sobj_new <- RunPCA(sobj_new, npcs = 30, verbose = F)
sobj_new <- FindNeighbors(sobj_new, dims = 1:30, reduction = "pca")
sobj_new <- FindClusters(sobj_new, resolution = 1, cluster.name = "clusters")
sobj_new <- RunUMAP(sobj_new, dims = 1:30, reduction = "pca", reduction.name = "umap")
qsave(sobj_new,"/home/garga7/Anjali_data/int_data/Dolan_data/iMGLs_depthnorm_liger_seurat_2024_recluster.qs" , nthreads = 6)
