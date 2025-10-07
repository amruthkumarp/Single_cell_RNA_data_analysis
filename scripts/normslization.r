#Load the packages
library(tidyverse) # dplyr and ggplot2
library(Seurat) # Seurat toolkit
library(hdf5r) # for data import
library(patchwork) # for plotting
library(presto) # for differential expression
library(glmGamPoi) # for sctransform
#load filtered data
adp_filt<-readRDS("adp_merge_filt.rds")


#applying the normalization
#sctransform
#adp_filt <- SCTransform(adp_filt, vars.to.regress = "percent.mt",
  #                     method = "glmGamPoi",return.only.var.genes = TRUE
# ,verbose = FALSE)

adp_filt <- NormalizeData(adp_filt)
adp_filt <- FindVariableFeatures(adp_filt)
adp_filt <- ScaleData(adp_filt, vars.to.regress = "percent.mt")

#Linear dimension reduction
Assays(adp_filt)
adp_filt <- RunPCA(adp_filt, verbose = FALSE,assay="RNA")

#Exploring the PCA results
#using dimheatmap
# dimensions 1 to 9
png("PCA_heatmap.png", width = 2000, height = 1500, res = 150)
DimHeatmap(adp_filt, dims = 1:9, cells = 500, balanced = TRUE, ncol = 3)
dev.off()

# dimensions 20 to 30
png("PCA_heatmap20_30.png", width = 2000, height = 1500, res = 150)
DimHeatmap(adp_filt, dims = 20:30, cells = 500, balanced = TRUE, ncol = 3)
dev.off()

#elbowplot
ElbowPlot(adp_filt, ndims = 40)

#Clustering
#knn distance
adp_filt <- FindNeighbors(adp_filt, dims = 1:30)
adp_filt <- FindClusters(adp_filt, resolution = 0.1)

#umap
adp_filt <- RunUMAP(adp_filt,dims = 1:30)

DimPlot(adp_filt, reduction = "umap", group.by = c("orig.ident", "seurat_clusters","time_point","condition","cond_tp"),
        alpha=0.4, ncol=2)

adp_filt

#integration step is skipped
#saving the seurat object
saveRDS(adp_filt, "adp_merge_filt_Nor_clust0.1.rds")


