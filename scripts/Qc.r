#Load the packages
library(tidyverse) # dplyr and ggplot2
library(Seurat) # Seurat toolkit
library(hdf5r) # for data import
library(patchwork) # for plotting
library(presto) # for differential expression
library(glmGamPoi) # for sctransform

#Load the Seurat Object

adp <- readRDS("/home/ibab/Downloads/GettingStarted_scRNASeq/GettingStarted_scRNASeq/outputs/merged_Seurat_adp.rds")

#Examine the object

glimpse(adp)

#Add percent mitochondrial
adp[["percent.mt"]] <- PercentageFeatureSet(adp, pattern = "^mt-")

#extract metadata
metadata <- adp@meta.data

#Visualize the data
#nCount_RNA
# set colors
cnames<-setNames(rep(c("cyan3","darkgoldenrod1"),each=4),levels(factor(metadata$orig.ident))) 

# plot
VlnPlot(adp, features = "nCount_RNA", layer="counts", group.by="orig.ident",raster=FALSE,alpha=0.2) +
  scale_fill_manual(values=cnames) 
#using ggplot
metadata %>% 
  ggplot(aes(color=orig.ident, x=nCount_RNA, fill= orig.ident)) + 
  geom_density(alpha = 0.2) + 
  theme_classic() +
  scale_x_log10() + 
  geom_vline(xintercept = 650,color="red",linetype="dotted")

#Number of Genes
VlnPlot(adp, features = "nFeature_RNA", group.by="orig.ident") +
  scale_fill_manual(values=cnames) 

ggplot(metadata, aes(x=nFeature_RNA,fill=orig.ident)) +
  geom_density(alpha = 0.2) + 
  theme_classic() +
  scale_x_log10() + 
  geom_vline(xintercept = 350,color="red",linetype="dotted")

#Percent mitochondrial
VlnPlot(adp, features = "percent.mt", group.by="orig.ident") +
  scale_fill_manual(values=cnames) +
  geom_hline(yintercept=10,color="red")
ggplot(metadata, aes(x=percent.mt,fill=orig.ident)) +
  geom_density(alpha = 0.2) + 
  scale_x_log10()+
  theme_classic()  

#Examine these metrics together
FeatureScatter(adp, feature1 = "percent.mt", feature2 ="nFeature_RNA" , group.by="orig.ident",split.by="time_point") 

FeatureScatter(adp, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",group.by="orig.ident",split.by="time_point",log=TRUE)

#ggplot
ggplot(metadata) +
  geom_point(aes(x=nCount_RNA,y=nFeature_RNA,fill=percent.mt > 10),shape=21,alpha=0.4) + 
  theme_classic() +
  scale_x_log10()+
  scale_y_log10()+
  facet_grid(.~cond_tp) +
  geom_vline(xintercept = 650,color="red",linetype="dotted")+
  geom_hline(yintercept=350,color="red", linetype="dotted")

#thershold for mt.genes is above 25
ggplot(metadata) +
  geom_point(aes(x=nCount_RNA,y=nFeature_RNA,fill=percent.mt > 25),shape=21,alpha=0.4) + 
  theme_classic() +
  scale_x_log10()+
  scale_y_log10()+http://127.0.0.1:34653/graphics/efb567a4-5620-411d-bdf9-bbcb537337fc.png
  facet_grid(.~cond_tp) +
  geom_vline(xintercept = 650,color="red",linetype="dotted")+
  geom_hline(yintercept=350,color="red", linetype="dotted")

#Decide on thresholds and apply filtering

  # Set one set of parameters for Day 0 samples; 
  # keep the rownames (Cell barcodes)
  t0<- metadata |> filter(time_point=="Day 0", nFeature_RNA > 350, 
                          nCount_RNA >650, percent.mt <10 ) |> 
    rownames_to_column("Cell") |> pull(Cell) 
  
  # Set an alternative set of thresholds for Day 6 samples; 
  # keep the rownames (Cell barcodes)
  t6<- metadata |> filter(time_point=="Day 6", nFeature_RNA > 350, 
                          nCount_RNA >650, percent.mt <25 ) |> 
    rownames_to_column("Cell") |> pull(Cell)
  
  keep<-c(t0,t6)

# Filter
  
  # use different parameters; established above
  adp_filt<-subset(adp, cells=keep)

#running again the same plots for qc check up
  
metadata_filt<-adp_filt@meta.data

# set colors
cnames<-setNames(rep(c("cyan3","darkgoldenrod1"),each=4),levels(factor(metadata_filt$orig.ident))) 

# plot
VlnPlot(adp_filt, features = "nCount_RNA", layer="counts", group.by="orig.ident",raster=FALSE,alpha=0.2) +
  scale_fill_manual(values=cnames) 
#using ggplot
metadata_filt %>% 
  ggplot(aes(color=orig.ident, x=nCount_RNA, fill= orig.ident)) + 
  geom_density(alpha = 0.2) + 
  theme_classic() +
  scale_x_log10() + 
  geom_vline(xintercept = 650,color="red",linetype="dotted")

#Number of Genes
VlnPlot(adp_filt, features = "nFeature_RNA", group.by="orig.ident") +
  scale_fill_manual(values=cnames) 

ggplot(metadata_filt, aes(x=nFeature_RNA,fill=orig.ident)) +
  geom_density(alpha = 0.2) + 
  theme_classic() +
  scale_x_log10() + 
  geom_vline(xintercept = 350,color="red",linetype="dotted")

#Percent mitochondrial
VlnPlot(adp_filt, features = "percent.mt", group.by="orig.ident") +
  scale_fill_manual(values=cnames) +
  geom_hline(yintercept=10,color="red")
#ggplot
ggplot(metadata, aes(x=percent.mt,fill=orig.ident)) +
  geom_density(alpha = 0.2) + 
  scale_x_log10()+
  theme_classic()  

#Examine these metrics together
FeatureScatter(adp_filt, feature1 = "percent.mt", feature2 ="nFeature_RNA" , group.by="orig.ident",split.by="time_point") 

FeatureScatter(adp_filt, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",group.by="orig.ident",split.by="time_point",log=TRUE)

#ggplot
ggplot(metadata_filt) +
  geom_point(aes(x=nCount_RNA,y=nFeature_RNA,fill=percent.mt > 10),shape=21,alpha=0.4) + 
  theme_classic() +
  scale_x_log10()+
  scale_y_log10()+
  facet_grid(.~cond_tp) +
  geom_vline(xintercept = 650,color="red",linetype="dotted")+
  geom_hline(yintercept=350,color="red", linetype="dotted")



#save the object
saveRDS(adp_filt,"adp_merge_filt.rds")


  

