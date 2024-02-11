# Pairwise comparison between species cell type clusters

### Introduction:
These scripts present a comprehensive comparison of orthologous marker genes for cell type clusters across three plant species, with the aim of generating lists of conserved orthologous marker genes that could be used to identify cell types in other plant species.

### Processing:
* Use the output generated from step 1 (MG_120923_Ath_05.RData) and step 2 (Orthogroups_091023_cleaned.tsv - the unzipped version of Orthogroups_091023_cleaned.tsv.zip) as input for this step. </br>
* Select the top 200 marker genes identified in step 1:
```R
Species_1 = merge(Marker_genes, Orthologous_group, by.x = "gene", by.y = "Species") %>% 
          arrange( cluster, desc(avg_log2FC)) %>% 
          select("gene", "Orthogroup", "cluster", "avg_log2FC") %>% 
          group_by(cluster) %>% 
          top_n(200)
```
* Count the number of common OMGs between two species cell type clusters:
```R
# Function create table of common OG gene names between 2 species
count_com_OG = function(Species1, Species2){
  clusters_S1 = unique(Species1$cluster); clusters_S2 = unique(Species2$cluster)
  two_plants = matrix(nrow=length(clusters_S1), ncol=length(clusters_S2))
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

df = as.data.frame(count_com_OG(Species_1, Species_2))
```

* Perform the fisher exact test to evaluate the significance of the conserved marker gene list for each cell type cluster:
```R
df$cell_type = rownames(df) # Add column to make 2 variables when using melt function
# Add a row that sums up the values in all other rows
sumrow = df %>% select(-cell_type) %>% colSums() # sum all numeric row in the dataframe
sum_h = c(sumrow, "sum_h") 
df = rbind(df, sum_h) # After merging two data frames, datatype in dataframe will be changed into character
library(hablar)
df = df %>% retype() %>% as.data.frame(df) # This library and function retype will change the data into the correct type

# Add a column to sum up all values in other columns
df$sum_v = df %>% select(-cell_type) %>% rowSums() # sum all numeric columns in the dataframe
rownames(df) = df$cell_type  

p_value_dataframe = df[1: (nrow(df) - 1), 1: (ncol(df) -2)]

for (i in rownames(p_value_dataframe)){
  for (j in colnames(p_value_dataframe)){
    frame = df[c(i, "sum_h"), c(j, "sum_v")]
    frame[2,2] = frame[2,2] - frame[1,2] - frame[2,1] + frame[1,1]
    frame[1,2] = frame[1,2] - frame[1,1]
    frame[2,1] = frame[2,1] - frame[1,1]
    p_value_dataframe[i,j] = fisher.test(frame, alternative = "greater")$p.value
  }
}
```

* Apply FDR correction to the p-values to address the issue of multiple testing: 
```R
adjusted_pvalue_dataframe[] = p.adjust(unlist(p_value_dataframe), method = "BH")
conclusionTable_0.01 = ifelse(adjusted_pvalue_dataframe < 0.01, "Reject", "Fail")
```

* For detailed information, please refer to the above scripts.
