# Integration

### Introduction:
This step aims to integrate gene expression data from three plant species to identify homologous cell type clusters among them.

### Requirements:
* R version 4.1.0 or higher </br>
* Matrix package </br>
* Seurat package version 4.1.1 </br>
* Harmony package

### Preprocessing:
* Species-specific genes were converted into orthologous genes, and 5000 high variable genes were identified and ranked from 1 to 5000 using the FindVariable function and nfeature=5000. </br> 
* Since one orthologous gene can consist of one or more species of genes, ranking the list of high-variable genes supported the selection of the top genes in each orthologous group for integration. Top 2000 unique high-variable orthologous genes corresponding to the top high-variable genes for each species were selected to ensure that the most informative genes were used in the integration process. 
* SCTransform was employed to normalize three datasets independently:
```R
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

DimPlot(crossSpecies.combine, reduction = "umap", group.by = "orig.ident", cols = c("blue", "green", "orange"))
```

* For detailed information, please refer to the above scripts.
