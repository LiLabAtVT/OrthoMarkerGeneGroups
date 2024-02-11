# Orthologous marker groups reveal broad cell identity conservation across single-cell transcriptomes in plants. 

### Introduction:
Despite the widespread application of single-cell RNA sequencing (scRNA-seq) inplant biology, the scarcity of known cell-type marker genes and the divergence of marker expression patterns limit the accuracy of cell-type identification in many species. To address this challenge, we have devised a novel computational strategy called Orthologous Marker Gene groups (OMGs) which can identify cell types in both model and non-model plant species. Our method does not depend on the complexity of cross-species data integration, thus is highly efficient, while still accurately determining inter species cellular similarities of diverse species. We validated our approach by analyzing published single-cell data from Arabidopsis,rice, and maize, and confirmed its accuracy in identifying cell types in tomato root and shoot apex tissues. The robustness of our method was further demonstrated by a successful mapping of 268 cell clusters from 1 million cells across 15 diverse plant species and various tissue types. Our findings suggest that the OMGs method, informed by reference single-cell maps, can accurately annotate cell types for most monocot and dicot species. To facilitate the use of this method by the broad research community, we have launched a user-friendly web-based tool called the OMG browser, which enables effortless identification of cell types in any plant dataset.

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
