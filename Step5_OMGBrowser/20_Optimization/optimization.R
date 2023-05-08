#This step is to extract the cluster,UMAP into a single RDS file to save
#time and execution time. We also extract the expression matrix information
#into separate RDS files. For running for multiple species just change the 
#paths in line 16(downsampled_rds).

#libraries
library(shiny)
library(Seurat)
library(SeuratData)
library(SeuratDisk)
library(ggplot2)
library(dplyr)

# Repeat this for all species
# Loading downsampled RDS data
downsampled_rds = readRDS("/Users/saipavanbathala/Documents/Li_lab/shinyOMG/10_Pre_Processing/downsampled/ath_data.rds")

# Extract cell IDs and cluster IDs from Seurat object
cell_ids <- rownames(downsampled_rds@meta.data)
cluster_ids <- downsampled_rds@active.ident

# Extract UMAP coordinates from Seurat object
umap_coords <- Embeddings(object = downsampled_rds, reduction = "umap")

# Combine cell IDs and UMAP coordinates into a data frame
umap_data <- data.frame(Cell = cell_ids, UMAP1 = umap_coords[,1], UMAP2 = umap_coords[,2])

#saving cluster data
cluster_info = merge(umap_data, data.frame(Cell = cell_ids, Cluster = cluster_ids), by = "Cell", all.x = TRUE)
saveRDS(cluster_info, "/Users/saipavanbathala/Documents/Li_lab/shinyOMG/20_Optimization/rds_files/ath_cluster_info.rds")

# Extract expression matrix from Seurat object
expr_mat <- downsampled_rds@assays$RNA@data

# Write expression matrix to text file
saveRDS(expr_mat, "/Users/saipavanbathala/Documents/Li_lab/shinyOMG/20_Optimization/rds_files/ath_expr_data.rds")

# Read count matrix, cluster assignments, and UMAP coordinates
expr_mat <- readRDS("/Users/saipavanbathala/Documents/Li_lab/shinyOMG/20_Optimization/rds_files/ath_expr_data.rds")
cluster_info <- readRDS("/Users/saipavanbathala/Documents/Li_lab/shinyOMG/20_Optimization/rds_files/ath_cluster_info.rds")

# Create UMAP CLUSTER's plot
ggplot(cluster_info, aes(x = UMAP1, y = UMAP2, color = Cluster)) +
  geom_point(size = 4) +
  theme_bw()