# Orthologous marker groups reveal broad cell identity conservation across single-cell transcriptomes in plants. 

### Introduction:
Despite the widespread application of single-cell RNA sequencing (scRNA-seq) inplant biology, the scarcity of known cell-type marker genes and the divergence of marker expression patterns limit the accuracy of cell-type identification in many species. To address this challenge, we have devised a novel computational strategy called Orthologous Marker Gene groups (OMGs) which can identify cell types in both model and non-model plant species. Our method does not depend on the complexity of cross-species data integration, thus is highly efficient, while still accurately determining inter species cellular similarities of diverse species. We validated our approach by analyzing published single-cell data from Arabidopsis,rice, and maize, and confirmed its accuracy in identifying cell types in tomato root and shoot apex tissues. The robustness of our method was further demonstrated by a successful mapping of 268 cell clusters from 1 million cells across 15 diverse plant species and various tissue types. Our findings suggest that the OMGs method, informed by reference single-cell maps, can accurately annotate cell types for most monocot and dicot species. 

![Test](./Figures/pipeline.jpg)

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
Identify marker genes in plant species using Seurat package
```R
# load packages
library(Seurat)

# load Seurat object
Seurat_obj <- readRDS("Seurat_object.rds")

# Find marker genes
Markers <- FindAllMarkers(Seurat_obj, only.pos = TRUE, min.pct = 0.5, 
              logfc.threshold = 0.5) %>% 
              group_by(cluster) %>% 
              arrange(cluster, desc(avg_log2FC))

# extract data 
saveRDS(Markers, "Markers.RData")
write.csv(Markers, "Markers.csv")
```

Detailed information can be found in this video tutorial: https://youtu.be/oliRmER1rXw
#### Step 2, Find orthologous groups by OrthoFinder:
Run the Orthofinder package and obtain a tsv file OrthoFinder_protein_folder/OrthoFinder/Results_date/Orthogroups/Orthogroups.tsv in linux

- `orthofinder -f OrthoFinder_protein_folder`


Convert species marker genes into orthologous genes
```R
# load the orthologous groups which is the output from the orthofinder package
orthofinder = read.csv("Step2_OrthoFinder/Orthogroups.tsv", header = TRUE, sep = "\t") # Unzip the orthogroups.tsv.zip

source("Step2_OrthoFinder/clean_orthofinder_output.R")
og_ath = clean_orthogroups_gene_names(orthofinder, "Arabidopsis") # EX: Arabidopsis
og_oryza = clean_orthogroups_gene_names(orthofinder, "Oryza") # EX: Rice (Oryza Sativa)

# load marker genes found in step 1, two examples of Arabidospis and Rice 
MG_ara = readRDS("Step1_FindMarkerGenes/MG_120923_ATH_05.RData")
MG_oryza = readRDS("Step1_FindMarkerGenes/MG_120923_Rice_05.RData")

# Merge marker gene and OG genes
source("Step2_OrthoFinder/species_omg.R")
Ath_MG_OG = Species_OMGs(MG_ara, og_ath, "Arabidopsis")  # EX: Arabidopsis 
Oryza_MG_OG = Species_OMGs(MG_oryza, og_oryza, "Oryza")  # EX: Rice 
```

#### Step 3, Find conserved orthologous marker genes (OMGs) between species:
Present a comprehensive pairwise comparison of orthologous marker genes for cell type clusters between two plant species, with the aim of generating lists of conserved orthologous marker genes that could be used to identify cell types in other plant species.  </br>

An example of Arabidopsis and Rice
```R
# Load packages
library(ggplot2); library(reshape); library(hablar)

# Function create table of common OMGs between 2 species
source("Step3_PairwiseComparison/count_comOMGs.R")
df_commonOMGs = test_significant(Ath_MG_OG, Oryza_MG_OG, 0.01)

# Create the plot with ggplot
source("Step3_PairwiseComparison/plot_generate.R")
plot_ATH_Rice = generate_plot_comparison(df_commonOMGs, "Arabidopsis", "Rice")
plot_ATH_Rice
```
![Test](./Figures/030823_Seurat_ATH_Rice_count2.jpg)

Detailed information can be found in the R markdown in step 3 folder

#### Step 4, Predict cell types in a query species:
Predict cell types in Tomato root and shoot.

```R
# Function create table of common OG gene names between 2 species
source("Step3_PairwiseComparison/count_comOMGs.R")

# Load packages
library(ggplot2)
library(reshape)
library(hablar)
df1 = test_significant(Ath_MG_OG, Tomato_MG_OG, 0.01)

# Create the plot with ggplot
source("Step3_PairwiseComparison/plot_generate.R")
plot_ATH_Tom = generate_plot_comparison(df1, "Arabidopsis", "Tomato")

plot_ATH_Tom
```
![Test](./Figures/ATH_Tom_prediction_root.png)
Based on the red boxes highlighted by statisical test, we can predict the cell type for tomato root using arabidopsis marker genes as a reference dataset. </br>
Detailed information for root and shoot prediction can be found in step 4 folder.

#### Mapping cell types across 15 species
This atlas comprises 268 cell clusters, 990,894 cells, and 53,600 marker genes which were obtained from multiple sources. We clustered cell types from different species based on the shared OMGs across species. One outcome of such mapping would be grouping similar cell types into clusters, neglecting any overriding effects of species dominance. 
![Test](./Figures/Map_15species.jpg)
The color bar is a group of a cell type measured by the euclidean distance. The color scale represents the negative logarithm of the FDR-adjusted p-value. Odd Ratios (OR) represent the likelihood of a particular cell type appearing in a specified group relative to its presence in all other groups. Proportion indicates the frequency of the predominant cell type within each group. Groups are named according to the most prevalent cell type, as indicated by their respective OR and proportion values.
Detailed information can be found in this file (Step4_cellType_prediction/Heatmap_cellType_comparison_across15species.Rmd)

#### Step 5, Web brower:
To facilitate the use of this method by the broad research community, we have launched a user-friendly web-based tool called the OMG browser, which enables effortless identification of cell types in any plant dataset.

Browser: http://orthomarkergenes.org/ </br>
Browser tutorial: https://youtu.be/Jb4uMq394Sg

### Contact:
If you have any questions regarding the OMG methods, please contact me at tnchau@vt.edu
