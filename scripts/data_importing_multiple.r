#load libraries
library(Seurat)
library(tidyverse)

#Importing Data

#To read multiple files
#List files 
files <- list.files(path="/home/ibab/data_analysis/Single_cell_rna_seq/seurat/data/",recursive=T,pattern="*.h5")

files
#Create a list of count matrices
h5_read <- lapply(paste0("/home/ibab/data_analysis/Single_cell_rna_seq/seurat/data/",
                         files), Read10X_h5)

#Automated way to assign names; modify for your purposes   
names(h5_read)<- sapply(files, 
                        function(x){str_split_1(x,"_")[1]},
                        USE.NAMES = FALSE)

#Assign names manually
names(h5_read)<-c("D1","D2","W1","W2")

# All samples
adp <- mapply(CreateSeuratObject,counts=h5_read,  
              project=names(h5_read),
              MoreArgs = list(min.cells = 3, min.features = 200))
# View adp 
adp

#remove the original sparse matrices
rm(h5_read)

#merging across seurat object
adp <- merge(adp[[1]], y = adp[2:length(adp)], 
             add.cell.ids = names(adp),project="Adipose")

adp

#Examine metadata
head(adp@meta.data) #using head to return only the first 6 rows  

#Accessing a column or columns:
#Access a single column. 
head(adp$orig.ident)

head(adp[["orig.ident"]])

#Access multiple columns, rows. 
head(adp[[c("orig.ident", "nCount_RNA")]])[1:3,]

#Add to metadata
adp$condition <- ifelse(str_detect(adp@meta.data$orig.ident, "^W"),
                        "WT","DKO")
#Add time information to the metadata.
adp$time_point <- ifelse(str_detect(adp@meta.data$orig.ident, "0"),
                         "Day 0","Day 6")

#Other Useful information from the Seurat object
head(Cells(adp,layer="counts.W1"))

head(colnames(adp))

#To return feature/gene names:
head(Features(adp))

head(rownames(adp))

#For multimodal data, list the assays using
Assays(adp)

#return number of cells across all layers
num_cells <- ncol(adp)
#return number of features across all layers
num_features <- nrow(adp)
#to return row by column information
dim(adp)

# Visualize the number of cell counts per sample
adp@meta.data %>% 
  ggplot(aes(x=orig.ident, fill=orig.ident)) + 
  geom_bar(color="black") +
  stat_count(geom = "text", colour = "black", size = 3.5, 
             aes(label = ..count..),
             position=position_stack(vjust=0.5))+
  theme_classic() +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("Number of Cells per Sample")

#save the object
saveRDS(adp,"merged_Seurat_adp.rds")
