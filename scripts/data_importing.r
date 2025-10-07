#loading library
library(Seurat)
library(tidyverse)

#data import
#Loading hdf5 files.
W10 <- Read10X_h5("20k_PBMC_3p_HT_nextgem_Chromium_X_filtered_feature_bc_matrix.h5")
W10[1:5,1:5]

#Creating a Seurat Object
# Single sample
W10 <- CreateSeuratObject(counts=W10, project="W10", min.cells = 3,
                          min.features = 200)
# View W10
W10

#Inspect and work with a Seurat object
glimpse(W10)

#To return the names of the layers:
Layers(W10)

#To access a specific count layer use:
  
  
  W10[["RNA"]]$counts |> head()
# or
W10@assays$RNA$counts  |> head()
# or
LayerData(W10, assay="RNA", layer='counts') |> head()

#You can always check the version of Seurat that was used to build the object with

W10@version

#And check commands and parameters used to generate results with:

W10@commands

#We haven't worked on the object, so there is nothing here. Compare with pbmc_small a built in data set in the SeuratObject package.
head(pbmc_small@commands, 1)

#Examine metadata
head(W10@meta.data) #using head to return only the first 6 rows  
#or
head(W10[[]])  

#Accessing a column or columns:
#Access a single column. 
head(W10$orig.ident)
head(W10[["orig.ident"]])

##Access multiple columns, rows. 
head(W10[[c("orig.ident", "nCount_RNA")]])[1:3,]

#Add to metadata.

#Add condition to metadata (either wildtype of double knockout).
W10$condition <- ifelse(str_detect(W10@meta.data$orig.ident, "^W"),
                        "WT","DKO")
#Add time information to the metadata.
W10$time_point <- ifelse(str_detect(W10@meta.data$orig.ident, "0"),
                         "Day 0","Day 6")

#Add Condition + Time to the metadata.
W10$cond_tp <- paste(W10$condition, W10$time_point)

#Other Useful information from the Seurat object
head(Cells(W10,layer="counts"))

head(colnames(W10))

#To return feature/gene names:
head(Features(W10))

head(rownames(W10))

#For multimodal data, list the assays using
Assays(W10)

#To return the number of cells or features use:
#return number of cells across all layers
num_cells <- ncol(W10)
num_cells
#return number of features across all layers
num_features <- nrow(W10)
num_features
#to return row by column information
dim(W10)

# Visualize the number of cell counts per sample
W10@meta.data %>% 
  ggplot(aes(x=orig.ident, fill=orig.ident)) + 
  geom_bar(color="black") +
  stat_count(geom = "text", colour = "black", size = 3.5, 
             aes(label = ..count..),
             position=position_stack(vjust=0.5))+
  theme_classic() +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("Number of Cells per Sample")
