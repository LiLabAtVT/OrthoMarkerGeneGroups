# Pairwise comparison between species cell type clusters

### Introduction:
These scripts present a comprehensive comparison of orthologous marker genes for cell type clusters across three plant species, with the aim of generating lists of conserved orthologous marker genes that could be used to identify cell types in other plant species.

### Requirements:
Use the output generated from step 1 and step 2 as input for this step. </br>
Select the top 200 marker genes identified in step 1:
```R
Species = merge(Marker_genes, Orthologous_group, by.x = "gene", by.y = "Species") %>% 
          arrange( cluster, desc(avg_log2FC)) %>% 
          select("gene", "Orthogroup", "cluster", "avg_log2FC") %>% 
          group_by(cluster) %>% 
          top_n(200)
```
