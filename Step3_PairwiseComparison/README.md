# Pairwise comparison between species cell type clusters

### Introduction:
These scripts present a comprehensive comparison of orthologous marker genes for cell type clusters across three plant species, with the aim of generating lists of conserved orthologous marker genes that could be used to identify cell types in other plant species.

### Requirements:
* Use the output generated from step 1 and step 2 as input for this step. </br>
* Select the top 200 marker genes identified in step 1:
```R
Species = merge(Marker_genes, Orthologous_group, by.x = "gene", by.y = "Species") %>% 
          arrange( cluster, desc(avg_log2FC)) %>% 
          select("gene", "Orthogroup", "cluster", "avg_log2FC") %>% 
          group_by(cluster) %>% 
          top_n(200)
```
* Count the number of common OMGs between two species cell type clusters:
```R
## Function create table of common OG gene names between 2 species
count_com_OG = function(Species1, Species2){
  clusters_S1 = unique(Species1$cluster); clusters_S2 = unique(Species2$cluster)
  two_plants <- matrix(nrow=length(clusters_S1), ncol=length(clusters_S2))
  for(i in 1:length(clusters_S1)){
    for(j in 1:length(clusters_S2)){
      list_overlap = intersect(Species1[Species1$cluster == clusters_S1[i],]$Orthogroup,  
                               Species2[Species2$cluster == clusters_S2[j],]$Orthogroup
                               )
      num_overlap = length(list_overlap) 
      two_plants[i,j] = num_overlap
    }
  }
  rownames(two_plants) = unique(Species1$cluster)
  colnames(two_plants) = unique(Species2$cluster) 
  return(two_plants[order(rownames(two_plants)), order(colnames(two_plants))])
}
```
