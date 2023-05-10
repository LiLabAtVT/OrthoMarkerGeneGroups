# Cross-species single cell cluster annotation with orthologous marker genes

### Introduction:
In this study, we compared single-cell RNA sequencing data from Arabidopsis, rice, and maize to identify conserved orthologous marker genes (OMGs) in various cell types. Our goal was to explore the significance of conserved OMGs for identifying cell types in non-model plants. 

### Requirements:
* R version 4.1.0 or higher </br>
* Matrix package </br>    
* Seurat package version 4.1.1 </br>
* SPmarker package and its dependencies including Python, pandas, sklearn, shap, keras </br>
* OrthoFinder software (v2.5.4) </br>
* Shiny package </br>
* Harmony package

### Processing:
#### Step 1, Find cell type marker genes:
Identifying marker genes in three plant species using three different methods: 
* Seurat </br>
* SHAP+RF </br>
* SVM </br>
#### Step 2, Find orthologous groups by OrthoFinder:
Convert species marker genes into orthologous genes by applying OrthoFinder

#### Step 3, Find conserved orthologous marker genes (OMGs) between species:
Present a comprehensive pairwise comparison of orthologous marker genes for cell type clusters between two plant species, with the aim of generating lists of conserved orthologous marker genes that could be used to identify cell types in other plant species.

#### Step 4, Extract the conserved OMGs of homologous cell types among three species:
Produce Venn diagrams that compare the shared OMGs among three plant species, with the purpose of generating lists of conserved OMGs that can aid in identifying cell types in other plant species.

#### Step 5, Web brower:
Visualize gene expression data across three different plant species </br>


### Reference:
