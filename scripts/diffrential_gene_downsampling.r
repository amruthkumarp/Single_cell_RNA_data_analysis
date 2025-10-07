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
adp <- readRDS("adp_filt_downsampling_clustered.rds")
mem.maxVSize()

glimpse(adp)

head(adp@meta.data)

#Reviewing the Data
table(adp$cond_tp)

table(adp$cond_tp,adp$orig.ident)

#UMAP
DimPlot(adp, reduction = "umap", group.by = c("orig.ident", "seurat_clusters","time_point","condition","cond_tp"),
        alpha=0.4, ncol=2)

#checking cluster
Idents(adp) = "RNA_snn_res.0.1"
table(Idents(adp))

#Differential Expression
#density plot of Gadph random 2500 cells
plot(density(sample(JoinLayers(adp@assays$RNA)$count["Gapdh",],2500)),cex=0,lty=1, main="Density of Gapdh in 2500 random cells")

#histogram plot
hist(sample(JoinLayers(adp@assays$RNA)$count["Gapdh",],2500),breaks=99,main="Histogram of Gapdh in 2500 random cells",ylab="Frequency",xlab="Gene counts")

#checking layers
Layers(adp[["RNA"]])


# Join all sample layers in the RNA assay
adp <- JoinLayers(adp, assay = "RNA")

#FindAllMarkers

de_allClusters = FindAllMarkers(adp, test.use="wilcox", min.pct=0.1, only.pos=TRUE)

#examining the resulta
head(de_allClusters)

nrow(de_allClusters)

#top 5 differentially expressed genes
top5PerCluster = matrix(ncol=7)
colnames(top5PerCluster) = colnames(de_allClusters)
for (i in 0:7){
  top5PerCluster = rbind(top5PerCluster, head(de_allClusters[which(de_allClusters$cluster==i),], 5))
}
top5PerCluster=top5PerCluster[-1,]
top5PerCluster

#do heatmap
DoHeatmap(adp,features = top5PerCluster$gene,slot="scale.data")

#we will be drawing a comparison between the two timepoints.
Idents(adp) = "time_point"
table(Idents(adp))

#And then we set up FindMarkers to compare our two conditions:
day6_day0_de=FindMarkers(adp,ident.1="Day 6",ident.2="Day 0",test.use="wilcox")
head(day6_day0_de, 10)

#visualization of the differential expression can be performed using the FeaturePlot function:
fig1 = DimPlot(adp,group.by="time_point")
fig2 = FeaturePlot(adp,features="Acta2",order=T)
fig3 = FeaturePlot(adp,features="Cd36",order=T)

fig1/(fig2|fig3)

# Create an output directory if not already existing
dir.create("./DEG_results", showWarnings = FALSE)
# Save cluster-wise top DEGs
write.csv(de_allClusters, 
          file = "./DEG_results/DEG_all_clusters.csv", 
          row.names = TRUE)

# Save top 5 per cluster
write.csv(top5PerCluster, 
          file = "./DEG_results/top5_genes_per_cluster.csv", 
          row.names = TRUE)

# Save time-point comparison (Day6 vs Day0)
write.csv(day6_day0_de, 
          file = "./DEG_results/DEG_Day6_vs_Day0.csv", 
          row.names = TRUE)

#load seruat object
adp_1 <- readRDS("adp_merge_filt_Nor_clust0.1.rds")

#pseduo aggregate
pseudo_adp=AggregateExpression(adp_1,assays="RNA", return.seurat=T, group.by=c("orig.ident","time_point","condition","cond_tp"))

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

scDE.genes = rownames(day6_day0_de)[which(day6_day0_de$p_val_adj<0.05)]
bulkDE.genes = rownames(bulk_adp_de)[which(bulk_adp_de$p_val_adj<0.05)]
length(scDE.genes)
length(bulkDE.genes)
length(intersect(scDE.genes,bulkDE.genes))
#checking for two targeted genes

bulk_adp_de[c("Acta2","Cd36"),]

#visualizing differential expressed genes
Idents(adp)="RNA_snn_res.0.1"
DotPlot(adp,features=unique(top5PerCluster$gene),dot.scale = 3)+coord_flip()

#volin plots
Idents(adp) = "time_point"
VlnPlot(adp,features=c("Acta2","Cd36"),alpha = 0.1)

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
Idents(adp)="RNA_snn_res.0.1"
avgExp = AverageExpression(adp,markers,assay="RNA")$RNA

avgExp

#dimensional plot
DimPlot(adp,label=T)

#feature plots
FeaturePlot(adp,features=markers,ncol=3,order=T)

#annotation
adipocyte=vector(length=ncol(adp))
adipocyte[which(adp$RNA_snn_res.0.1 %in% c(0,5))]="Preadipocytes"
adipocyte[which(adp$RNA_snn_res.0.1 %in% c(2,6))]="Proliferating cells"
adipocyte[which(adp$RNA_snn_res.0.1 %in% c(1,3))]="Differentiating beige adipocytes"
adipocyte[which(adp$RNA_snn_res.0.1 %in% c(4))]="Differentiated beige adipocytes"
adipocyte[which(adp$RNA_snn_res.0.1 %in% c(7))]="Unclassified"
adp$adipocyte = adipocyte

f1 = DimPlot(adp,group.by="RNA_snn_res.0.1",label=T) + NoLegend()
f2 = DimPlot(adp,group.by="time_point") + NoLegend()
f3 = DimPlot(adp,group.by="adipocyte",label=T)+NoLegend()

(f1|f2)/f3

#Annotation using SingleR
