# Identify marker genes using three methods

This repository contains R script files for identifying marker genes in Arabidopsis using Seurat package

### Requirements: 
* R version 4.1.0 or higher </br>
* Matrix package </br>    
* Seurat package version 4.1.1 </br>
* SPmarker package and its dependencies including Python, pandas, sklearn, shap, keras </br>
`git clone https://github.com/LiLabAtVT/SPMarker.git`

### Preprocessing:
Preprocess data to generate clusters and identify cell types for each cluster: 

```R
# Load the input dataset
matrix = readMM(file = "matrix.mtx")
rownames(matrix) = read.delim(file = "gene.tsv", header = FALSE)
colnames(matrix) = read.delim(file = "barcode.tsv", header = FALSE)

Species_seurat <- CreateSeuratObject(matrix, project = "Plant")
VlnPlot(Species_seurat, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2) # quality control based on these plots
Species_seurat <- subset(Species_seurat, subset = nFeature_RNA > 500 & nFeature_RNA < 6000 & nCount_RNA > 500 & nCount_RNA < 40000)

# Normal and scale the data
Species_seurat <- NormalizeData(object = Species_seurat, normalization.method = "LogNormalize")
Species_seurat <- FindVariableFeatures(object = Species_seurat, selection.method = "vst")
Species_seurat <- ScaleData(Species_seurat,  features = rownames(Species_seurat))

# Reduce the dimention of the data
Species_seurat <- RunPCA(Species_seurat, features = VariableFeatures(object = Species_seurat) ,verbose = FALSE)
Species_seurat <- FindNeighbors(Species_seurat, reduction = "pca", dims = parameter) 
Species_seurat <- FindClusters(Species_seurat, resolution = parameter) 
Species_seurat <- RunUMAP(Species_seurat, reduction = "pca", dims = parameter) 

# Plot UMAP
DimPlot(Species_seurat, reduction = "umap", label = TRUE)

# Visualize marker genes from reference papers to label cell types in the input dataset
DotPlot(object = Species_seurat, features = c("genes from reference papers"), cols = "RdYlBu",  col.min= -2, col.max = 2, dot.scale = 4) + 
  theme(axis.text.x = element_text(size=9, angle = 90, hjust=1), 
        axis.text.y = element_text(size=10, angle = 0, hjust=1), 
        axis.title.y  = element_text(size=15, angle = 90, vjust=-4),
        legend.title = element_text(size=10),
        legend.text = element_text(size=8)) + 
  xlab('Gene') +  
  ylab('Cell type')

FeaturePlot(Species_seurat, features = "genes from reference papers")
```

### Find marker genes by Seurat package:
```R
Species_seurat_marker <- FindAllMarkers(Species_seurat, only.pos = TRUE, min.pct = params, logfc.threshold = params) %>% 
                               group_by(cluster) %>% 
                                arrange(cluster, desc(avg_log2FC))
```
For detailed information, please refer to these scripts above:
* 10_Umap_120622_Ara.Rmd 

### Result:
The output data is stored in the file named "MG_120923_Ath_05.RData".

MG_120923_Maize_05.RData and MG_120923_Rice_05.RData contain marker genes of maize and rice.

### References:
* Ryu, K. H., Huang, L., Kang, H. M., & Schiefelbein, J. (2019). Single-cell RNA sequencing resolves molecular relationships among individual plant cells. Plant physiology, 179(4), 1444-1456. </br>
* Stuart, T. et al. Comprehensive Integration of Single-Cell Data. Cell 177, 1888-1902.e21 (2019). 
