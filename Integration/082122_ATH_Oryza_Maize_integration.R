
library(Seurat); library(Matrix); library(tidyverse)

#***************************************#
#
# Load ATH gene expression matrix
fileNames = c( "Sample_WT-WERGFP", "Sample_WT-WERGFP_2", "Sample_WT-WERGFP_3")
data_sample = c()
for(i in fileNames){
  matrix = readMM(file = paste0("/Users/tranchau/Documents/SC_crossSpecies/GSE123013_RAW_matrices/", i, "/filtered_gene_bc_matrices/TAIR10/matrix.mtx"))
  rownames(matrix) = read.delim(file = paste0("/Users/tranchau/Documents/SC_crossSpecies/GSE123013_RAW_matrices/", i, "/filtered_gene_bc_matrices/TAIR10/genes.tsv"), header = FALSE, stringsAsFactors = FALSE)$V1
  colnames(matrix) = read.delim(file = paste0("/Users/tranchau/Documents/SC_crossSpecies/GSE123013_RAW_matrices/", i, "/filtered_gene_bc_matrices/TAIR10/barcodes.tsv"), header = FALSE, stringsAsFactors = FALSE)$V1
  data_sample <- c(data_sample, matrix)
}
ATH_data = cbind(data_sample[[1]], data_sample[[2]], data_sample[[3]]) # Same rownames, different columns,
#
# Load ATH genes and their OG genes
ATHgene_OG = readRDS("/Users/tranchau/Documents/SC_crossSpecies/Species_OG_ATH_082022.RData")[1:2000,]

ATH_w_OG = merge(ATHgene_OG, ATH_data, by.x = "HVG_ara", by.y = 0)[, -c(1,2)] %>% column_to_rownames("Orthogroup")
ATH_seurat = CreateSeuratObject(ATH_w_OG, project = "ATH")

#***************************************#
#
# Load Maize gene expression
matrix = readMM(file = "/Users/tranchau/Documents/SC_crossSpecies/E-ENAD-52-quantification-raw-files/E-ENAD-52.aggregated_filtered_counts.mtx")
rownames(matrix) = read.delim(file = "/Users/tranchau/Documents/SC_crossSpecies/E-ENAD-52-quantification-raw-files/E-ENAD-52.aggregated_filtered_counts.mtx_rows", header = FALSE)$V1
colnames(matrix) = read.delim(file = "/Users/tranchau/Documents/SC_crossSpecies/E-ENAD-52-quantification-raw-files/E-ENAD-52.aggregated_filtered_counts.mtx_cols", header = FALSE)$V1
#
# Load Maize genes with their OG genes
Oryzagene_OG = readRDS("/Users/tranchau/Documents/SC_crossSpecies/Species_OG_Rice_082022.RData")[1:2000,]
#
Oryza_w_OG = merge(Oryzagene_OG, matrix, by.x = "HVG_oryza", by.y = 0)[, -c(1,2)] %>% column_to_rownames("Orthogroup")
Oryza_seurat = CreateSeuratObject(Oryza_w_OG, project = "Oryza")

#***************************************#
#
library(readr)
GSE173087_Maize <- read.csv("/Users/tranchau/Documents/SC_crossSpecies/GSE173087_Maize_cells_expression_matrix.csv", row.names = 1, header= TRUE)
#
# Load Maize genes and their OG genes
Maizegene_OG = readRDS("/Users/tranchau/Documents/SC_crossSpecies/Species_OG_Maize_082022.RData")[1:2000,]

Maize_w_OG = merge(Maizegene_OG, GSE173087_Maize, by.x = "HVG_maize", by.y = 0)[, -c(1,2)] %>% column_to_rownames("Orthogroup")
Maize_seurat = CreateSeuratObject(Maize_w_OG, project = "Maize")
Maize_seurat@meta.data$orig.ident = "Zeamays"


#***************************************#
#
#
crossSpecies <- SplitObject(merge(ATH_seurat, y= c(Oryza_seurat, Maize_seurat)) , split.by = "orig.ident")
for (i in names(crossSpecies)) {
  crossSpecies[[i]] <- SCTransform(crossSpecies[[i]], verbose = FALSE)
}
features <- SelectIntegrationFeatures(object.list = crossSpecies, nfeatures = 10000)
print(length(features))
crossSpecies <- PrepSCTIntegration(object.list = crossSpecies, anchor.features = features)
crossSpecies.anchors <- FindIntegrationAnchors(object.list = crossSpecies, normalization.method = "SCT", anchor.features = features)
crossSpecies.combine <- IntegrateData(anchorset = crossSpecies.anchors, normalization.method = "SCT")
crossSpecies.combine <- RunPCA(crossSpecies.combine, verbose = FALSE)
crossSpecies.combine <- RunUMAP(crossSpecies.combine, reduction = "pca", dims = 1:30)
crossSpecies.combine <- FindNeighbors(crossSpecies.combine, reduction = "pca", dims = 1:30)
crossSpecies.combine <- FindClusters(crossSpecies.combine, resolution = 0.9) #org =0.3
#saveRDS(crossSpecies.combine, "/Users/tranchau/Documents/SC_crossSpecies/042423_AOZ_inteObj.Rds")

#pdf("/Users/tranchau/Documents/SC_crossSpecies/042423_ATH_Oryza_Maize.pdf", width = 7, height = 5)
DimPlot(crossSpecies.combine, reduction = "umap", group.by = "orig.ident", cols = c("blue", "green", "orange"))
#dev.off()

#pdf("/Users/tranchau/Documents/SC_crossSpecies/042423_ATH_Oryza_Maize_sepColor.pdf", width = 7, height = 5)
DimPlot(crossSpecies.combine, reduction = "umap", split.by = "orig.ident", group.by = "orig.ident", cols = c("blue", "gray", "orange"))
#dev.off()

#pdf("/Users/tranchau/Documents/SC_crossSpecies/042423_ATH_Oryza_Maize_lable.pdf", width = 7, height = 5)
DimPlot(crossSpecies.combine, reduction = "umap", label = TRUE)
#dev.off()


#pdf("/Users/tranchau/Documents/SC_crossSpecies/042423_ATH_Oryza_Maize_sep.pdf", width = 7, height = 5)
DimPlot(crossSpecies.combine, reduction = "umap", split.by = "orig.ident")
#dev.off()

