
library(Seurat); library(Matrix); library(tidyverse); library(harmony)

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

#ATH_seurat <- NormalizeData(object = ATH_seurat, verbose = FALSE)
#ATH_seurat <- FindVariableFeatures(object = ATH_seurat, selection.method = "vst", verbose = FALSE)
#ATH_seurat <- ScaleData(ATH_seurat,  features = rownames(ATH_seurat))
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

#Oryza_seurat <- subset(Oryza_seurat, subset = nFeature_RNA > 500 & nFeature_RNA < 6000 & nCount_RNA > 500 & nCount_RNA < 40000 )

#Oryza_seurat <- NormalizeData(object = Oryza_seurat, normalization.method = "LogNormalize")
#Oryza_seurat <- FindVariableFeatures(object = Oryza_seurat, selection.method = "vst")
#Oryza_seurat <- ScaleData(Oryza_seurat,  features = rownames(Oryza_seurat))
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

#Maize_seurat <- NormalizeData(object = Maize_seurat, verbose = FALSE)
#Maize_seurat <- FindVariableFeatures(object = Maize_seurat, selection.method = "vst")
#Maize_seurat <- ScaleData(Maize_seurat, features = rownames(Maize_seurat))
#***************************************#
#
#

crossSpecies <- SplitObject(merge(ATH_seurat, y= c(Oryza_seurat, Maize_seurat)) , split.by = "orig.ident")
for (i in names(crossSpecies)) {
  crossSpecies[[i]] <- SCTransform(crossSpecies[[i]], verbose = FALSE)
}
features <- SelectIntegrationFeatures(object.list = crossSpecies, nfeatures = 10000)

crossSpecies.combine <- merge(crossSpecies[[1]], y= c(crossSpecies[[2]], crossSpecies[[3]]))
VariableFeatures(crossSpecies.combine) <- features

crossSpecies.combine <- RunPCA(crossSpecies.combine, assay = "SCT")


crossSpecies.combine <- RunHarmony(crossSpecies.combine, reduction = "pca", group.by.vars = "orig.ident", assay.use = "SCT", plot_convergence = TRUE) 


crossSpecies.combine <- RunUMAP(crossSpecies.combine, assay = "SCT", reduction = "harmony", dims = 1:30)

# Find neighbors in the integrated data
crossSpecies.combine <- FindNeighbors(crossSpecies.combine, assay = "SCT", reduction = "harmony", dims = 1:30)

# Find clusters based on the integrated data
crossSpecies.combine <- FindClusters(crossSpecies.combine, resolution = 0.2)

saveRDS(crossSpecies.combine, "/Users/tranchau/Documents/SC_crossSpecies/082122_AOZ_inteObj_Harmony_032323.Rds")

pdf("/Users/tranchau/Documents/SC_crossSpecies/082122_ATH_Oryza_Maize_Harmony_032323.pdf", width = 7, height = 5)
DimPlot(crossSpecies.combine, reduction = "umap", group.by = "orig.ident", cols = c("blue", "green", "orange"))
dev.off()

pdf("/Users/tranchau/Documents/SC_crossSpecies/082122_ATH_Oryza_Maize_sepColor_Harmony_032323.pdf", width = 7, height = 5)
DimPlot(crossSpecies.combine, reduction = "umap", split.by = "orig.ident", group.by = "orig.ident", cols = c("blue", "gray", "orange"))
dev.off()

pdf("/Users/tranchau/Documents/SC_crossSpecies/082122_ATH_Oryza_Maize_lable_Harmony_032323.pdf", width = 7, height = 5)
DimPlot(crossSpecies.combine,  reduction = "umap", label = TRUE)
dev.off()


pdf("/Users/tranchau/Documents/SC_crossSpecies/082122_ATH_Oryza_Maize_sep_Harmony_032323.pdf", width = 7, height = 5)
DimPlot(crossSpecies.combine, reduction = "umap", split.by = "orig.ident")
dev.off()

