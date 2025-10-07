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

#pseduo aggregate
pseudo_adp=AggregateExpression(adp,assays="RNA", return.seurat=T, group.by=c("orig.ident","time_point","condition","cond_tp"))

head(pseudo_adp@assays$RNA$counts)
pseudo_adp@meta.data

#just to clean up the look a little bit
pseudo_adp = RenameCells(pseudo_adp,new.names=gsub("_.*","",pseudo_adp$orig.ident))
pseudo_adp$orig.ident=gsub("_.*","",pseudo_adp$orig.ident)
head(pseudo_adp@assays$RNA$counts)

#just to look for
pseudo_adp@meta.data

#running Deseq2
Idents(pseudo_adp)="time_point"
bulk_adp_de = FindMarkers(pseudo_adp, ident.1="Day 6", ident.2="Day 0", test.use="DESeq2")

head(bulk_adp_de)

# Cell Annotation
#From the paper where this data was obtained, the following (incomplete) list of gene markers was obtained:
#Mmp3: preadipocytes
#Mki67: proliferating cells
#Fabp4: differentiating beige adipocytes and differentiated beige adipocytes
#Scd1: differentiated beige adipocytes
#Ucp1: differentiated beige adipocytes
#Ppargc1a: differentiated beige adipocytes
#Elovl3: differentiated beige adipocytes
#Cidea: differentiated beige adipocytes

#average expression of these genes 
markers=c("Mmp3","Mki67","Fabp4","Scd1","Ucp1","Ppargc1a","Elovl3","Cidea")
Idents(adp)="SCT_snn_res.0.1"
avgExp = AverageExpression(adp,markers,assay="SCT")$SCT

