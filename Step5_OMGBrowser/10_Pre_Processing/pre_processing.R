#This file is for down sampling the raw RDS files of each species and 
#to plot the UMAP of down sampled data

#Loading libraries
library(shiny)
library(Seurat)
library(SeuratData)
library(SeuratDisk)
library(ggplot2)
library(dplyr)

#For ATH Data

# Load RAW RDS data
ath_data = readRDS("/Users/saipavanbathala/Documents/Li_lab/shinyOMG/Raw_RDS_files/ATH_wLabel.rds")

#Normalizing data before processing
ath_data <- NormalizeData(ath_data)

#Finding top 3000 genes that covers most variance
ath_data_2 = FindVariableFeatures(ath_data, selection.method = "vst", nfeatures = 3000)

# Subset Seurat object using top 1000 variable features
top_features = head(VariableFeatures(ath_data_2), n = 1000)

#Adding OMG genes
gene_file <- readxl::read_excel("/Users/saipavanbathala/Documents/Li_lab/shinyOMG/Raw_OMG_Mapping_files/OMG.xlsx", sheet = "ATH_Seurat")

# Extract the genes from the desired column (assuming it is named "v1")
additional_genes <- gene_file$V1

# Combine the additional genes with the existing top_features
combined_features <- c(top_features, additional_genes)

# Subset Seurat object using the combined features
ath_data_3 <- ath_data_2[combined_features, ]

# Subset 300 cells
set.seed(123)
ath_downsampled_data <- ath_data_3[, sample(ncol(ath_data_3), size = 300)]






# Save downsampled Seurat object as RDS
saveRDS(ath_downsampled_data, file = "/Users/saipavanbathala/Documents/Li_lab/shinyOMG/10_Pre_Processing/downsampled/ath_data.rds")

#plotting for verification of overall structure
ath_data_plot <- readRDS("/Users/saipavanbathala/Documents/Li_lab/shinyOMG/10_Pre_Processing/downsampled/ath_data.rds")
nrow(ath_data_plot)
ncol(ath_data_plot)
DimPlot(ath_data_plot, reduction = "umap", label = TRUE, pt.size = 4)



#For rice Data

# Load RDS data
rice_data = readRDS("/Users/saipavanbathala/Documents/Li_lab/shinyOMG/Raw_RDS_files/Rice_wLabel.rds")

#Normalizing data before processing
rice_data <- NormalizeData(rice_data)

#Finding top 3000 genes that covers most varinace
rice_data_2 = FindVariableFeatures(rice_data, selection.method = "vst", nfeatures = 3000)

# Subset Seurat object using top variable features
top_features = head(VariableFeatures(rice_data_2), n = 1000)

#Adding OMG genes
gene_file <- readxl::read_excel("/Users/saipavanbathala/Documents/Li_lab/shinyOMG/Raw_OMG_Mapping_files/OMG.xlsx", sheet = "Rice_Seurat")

# Extract the genes from the desired column (assuming it is named "v1")
additional_genes <- gene_file$V1

# Combine the additional genes with the existing top_features
combined_features <- c(top_features, additional_genes)

rice_data_3 = rice_data_2[combined_features, ]

# Subset 300 cells
set.seed(123)
rice_downsampled_data <- rice_data_3[, sample(ncol(rice_data_3), size = 300)]

# Save downsampled Seurat object as RDS
saveRDS(rice_downsampled_data, file = "/Users/saipavanbathala/Documents/Li_lab/shinyOMG/10_Pre_Processing/downsampled/rice_data.rds")

#plotting
rice_data_plot <- readRDS("/Users/saipavanbathala/Documents/Li_lab/shinyOMG/10_Pre_Processing/downsampled/rice_data.rds")
nrow(rice_data_plot)
ncol(rice_data_plot)
DimPlot(rice_data_plot, reduction = "umap", label = TRUE, pt.size = 4)



#For maize Data

# Load RDS data
maize_data = readRDS("/Users/saipavanbathala/Documents/Li_lab/shinyOMG/Raw_RDS_files/Maize_wLabel.rds")

#Normalizing data before processing
maize_data <- NormalizeData(maize_data)

#Finding top 3000 genes that covers most varinace
maize_data_2 = FindVariableFeatures(maize_data, selection.method = "vst", nfeatures = 3000)

# Subset Seurat object using top variable features
top_features = head(VariableFeatures(maize_data_2), n = 1000)

#Adding OMG genes
gene_file <- readxl::read_excel("/Users/saipavanbathala/Documents/Li_lab/shinyOMG/Raw_OMG_Mapping_files/OMG.xlsx", sheet = "Maize_Seurat")

# Extract the genes from the desired column (assuming it is named "v1")
additional_genes <- gene_file$V1

# Combine the additional genes with the existing top_features
combined_features <- c(top_features, additional_genes)
maize_data_3 = maize_data_2[combined_features, ]

# Subset 300 cells
set.seed(123)
maize_downsampled_data <- maize_data_3[, sample(ncol(maize_data_3), size = 300)]

# Save downsampled Seurat object as RDS
saveRDS(maize_downsampled_data, file = "/Users/saipavanbathala/Documents/Li_lab/shinyOMG/10_Pre_Processing/downsampled/maize_data.rds")

#plotting
maize_data_plot <- readRDS("/Users/saipavanbathala/Documents/Li_lab/shinyOMG/10_Pre_Processing/downsampled/maize_data.rds")
nrow(maize_data_plot)
ncol(maize_data_plot)
DimPlot(maize_data_plot, reduction = "umap", label = TRUE, pt.size = 4)