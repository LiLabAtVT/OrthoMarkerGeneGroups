# Cross-species single cell cluster annotation with orthologous marker gene groups

### Introduction:
Despite the widespread application of single-cell RNA sequencing (scRNA-seq) inplant biology, the scarcity of known cell-type marker genes and the divergence ofmarker expression patterns limit the accuracy of cell-type identifi cation in manyspecies. To address this challenge, we have devised a novel computational strategycalled Orthologous Marker Gene groups (OMGs) which can identify cell types inboth model and non-model plant species. Our method does not depend on thecomplexity of cross-species data integration, thus is highly effi cient, while stillaccurately determining interspecies cellular similarities of diverse species. Wevalidated our approach by analyzing published single-cell data from Arabidopsis,rice, and maize, and confi rmed its accuracy in identifying cell types in tomato rootand shoot apex tissues. The robustness of our method was further demonstrated bya successful mapping of 268 cell clusters from 1 million cells across 15 diverseplant species and various tissue types. Our fi ndings suggest that the OMGsmethod, informed by reference single-cell maps, can accurately annotate cell typesfor most monocot and dicot species. To facilitate the use of this method by thebroad research community, we have launched a user-friendly web-based tool calledthe OMG browser, which enables effortless identifi cation of cell types in any plantdataset.

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
Identifying marker genes in three plant species using Seurat package
Tutorial: https://youtu.be/oliRmER1rXw

#### Step 2, Find orthologous groups by OrthoFinder:
Convert species marker genes into orthologous genes by applying OrthoFinder

#### Step 3, Find conserved orthologous marker genes (OMGs) between species:
Present a comprehensive pairwise comparison of orthologous marker genes for cell type clusters between two plant species, with the aim of generating lists of conserved orthologous marker genes that could be used to identify cell types in other plant species.

#### Step 4, Predict cell types in a query species:
Predict cell types in Tomato root and shoot.

#### Step 5, Web brower:
Perform cell types comparison across 15 plant species </br>
Browser: http://orthomarkergenes.org/ </br>
Browser tutorial: https://youtu.be/Jb4uMq394Sg

### Contact:
If you have any questions regarding the OMG methods, please contact me at tnchau@vt.edu
