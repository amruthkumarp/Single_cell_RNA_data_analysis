#reading single cell data into seurat object
#loading library
library(Seurat)
library(SeuratDisk)

#Different formats of single cell data can be stored
#.RDS format
rds_obj<-readRDS("ependymal_cells.rds")
str(rds_obj)

#.hdf5 format

hd5_obj <- Read10X_h5(
  filename = "20k_PBMC_3p_HT_nextgem_Chromium_X_filtered_feature_bc_matrix.h5",
  use.names = TRUE,
  unique.features = TRUE
)
head(hd5_obj)

#convert seurat obj
seurat_hdf5<-CreateSeuratObject(counts=hd5_obj)
seurat_hdf5

#.mtx file
mtx_obj<-ReadMtx(mtx="/home/ibab/data_analysis/Single_cell_rna_seq/seurat/raw_feature_bc_matrix/matrix/matrix.mtx.gz",
        features="/home/ibab/data_analysis/Single_cell_rna_seq/seurat/raw_feature_bc_matrix/matrix/features.tsv.gz",
        cells="/home/ibab/data_analysis/Single_cell_rna_seq/seurat/raw_feature_bc_matrix/matrix/barcodes.tsv.gz")
seurat_mtx<-CreateSeuratObject(counts=mtx_obj)
seurat_mtx

#.h5ad files
#step 1 :convert Anndata object to an h5 seurat file
Convert("adata_SS2_for_download.h5ad",dest="h5seurat",overwrite=TRUE)

#step2:load h5seurat file into seurat object
seurat_ann<-LoadH5Seurat("adata_SS2_for_download.h5seurat")
