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


set.seed(123)

# parameters
N <- 300   # target cells per sample-per-cluster (adjust 100-500 as needed)

# metadata column names (change if different)
sample_col <- "orig.ident"        # e.g. "D10_WT_rep1" or "D10"
cluster_col <- "seurat_clusters"  # or "cluster"

head(adp@meta.data)
# build a vector of cells to keep (downsampled)
cells_keep <- adp@meta.data %>%
  mutate(cell = rownames(.)) %>%
  group_by(.data[[sample_col]], .data[[cluster_col]]) %>%
  group_map(~ {
    cells <- .x$cell
    if(length(cells) <= N) return(cells)
    sample(cells, N)
  }) %>% unlist()

# subset Seurat object if you performed downsampling
adp_sub <- subset(adp, cells = cells_keep)
# If you skip downsampling, set adp_sub <- adp_filt

saveRDS(adp_sub,"adp_filt_downsampling_clustered.rds")
