#load packages

library(tidyverse) # dplyr and ggplot2; CRAN
library(Seurat) # Seurat toolkit; CRAN
library(hdf5r) # for data import; CRAN
library(patchwork) # for plotting; CRAN 
library(presto) # for differential expression; Github
library(glmGamPoi) # for sctransform; Bioconductor
library(ggplot2) # for visualization; installed with tidyverse
library(dplyr) #just to make sure it was called; installed with tidyverse


library(SingleR) # for cell type annotation; Bioconductor
library(celldex) # for cell type annotation reference; Bioconductor
library(MAST) # for differential expression; Bioconductor

#load seruat object
adp <- readRDS("adp_merge_filt_Nor_clust0.1.rds")
mem.maxVSize()

mem.maxVSize(vsize=16384*4)

#examine the object
glimpse(adp)

#reviving data
table(adp$cond_tp)

#2d data
table(adp$cond_tp,adp$orig.ident)

DimPlot(adp, reduction = "umap", group.by = c("orig.ident", "seurat_clusters","time_point","condition","cond_tp"),
        alpha=0.4, ncol=2)

Idents(adp) = "SCT_snn_res.0.1"
table(Idents(adp))


plot(density(sample(JoinLayers(adp@assays$RNA)$count["Gapdh",],2500)),cex=0,lty=1, main="Density of Gapdh in 2500 random cells")

#checking layers
Layers(adp[["RNA"]])


# Join all sample layers in the RNA assay
adp <- JoinLayers(adp, assay = "RNA")

#FindAllMarkers

de_allClusters = FindAllMarkers(adp, test.use="wilcox", min.pct=0.1, only.pos=TRUE)
