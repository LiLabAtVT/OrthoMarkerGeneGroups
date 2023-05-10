# Identify marker genes using three methods

This repository contains scripts file for identifying marker genes in three plant species using three different methods: Seurat, SHAP+RF, and SVM.

### Install: 
* R version 4.1.0 or higher </br>
* Matrix package </br>    
* Seurat package version 4.0 </br>
* SPmarker package and its dependencies including Python, pandas, sklearn, shap, keras </br>
`git clone https://github.com/LiLabAtVT/SPMarker.git`

### Preprocessing:
Label cell types for each cluster by employing marker genes from reference papers
```R
# Load the input dataset
matrix = readMM(file = "matrix.mtx")
rownames(matrix) = read.delim(file = "gene.tsv", header = FALSE)
colnames(matrix) = read.delim(file = "barcode.tsv", header = FALSE)

Species_seurat <- CreateSeuratObject(matrix, project = "Plant")
VlnPlot(Oryza_seurat, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2) # quality control based on these plots
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

### Find marker genes by Seurat function:
```R
Species_seurat_marker <- FindAllMarkers(Species_seurat, only.pos = TRUE, min.pct = params, logfc.threshold = params) %>% 
                               group_by(cluster) %>% 
                                arrange(cluster, desc(avg_log2FC))
```
For detailed information, please refer to these scripts above:
* 10_Umap_120622_Ara.Rmd 
* 10_Umap_120622_Maize.Rmd
* 10_Umap_120622_Rice.Rmd

### Find marker genes by SHAP+RF and SVM:
Prepare a gene expression matrix file (.csv) and cell meta file (.csv) </br>
`python SPmarker/SPmarker.py \ 
            -d work_directory/ -o work_directory/ \  
            -mtx gene_expression.csv \ 
            -meta cellType.csv`

### Result:
The collected output data is stored in the file named "112122_top200_Seurat_SHAP_SVM.xlsx".

### References:
Ryu, K. H., Huang, L., Kang, H. M., & Schiefelbein, J. (2019). Single-cell RNA sequencing resolves molecular relationships among individual plant cells. Plant physiology, 179(4), 1444-1456.
Ortiz-Ramírez, C., Dias Araujo, P. C., Zhang, S., Demesa-Arevalo, E., Yan, Z., Xu, X., ... & Birnbaum, K. D. (2021). Ground tissue circuitry regulates organ complexity in cereal roots. BioRxiv, 2021-04.
Zhang, T. Q., Chen, Y., Liu, Y., Lin, W. H., & Wang, J. W. (2021). Single-cell transcriptome atlas and chromatin accessibility landscape reveal differentiation trajectories in the rice root. Nature communications, 12(1), 2053.
Yan, H., Lee, J., Song, Q., Li, Q., Schiefelbein, J., Zhao, B., & Li, S. (2022). Identification of new marker genes from plant single‐cell RNA‐seq data using interpretable machine learning methods. New Phytologist, 234(4), 1507-1520.
