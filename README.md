# Cross-species single cell cluster annotation with orthologous marker genes

### Introduction:


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

#### Step 4, Extract the conserved OMGs of homologous cell types among three species:

#### Step 5, Web brower:
