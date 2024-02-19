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
```{r}
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
Run the Orthofinder package and obtain a tsv file OrthoFinder_protein_folder/OrthoFinder/Results_date/Orthogroups/Orthogroups.tsv
```{}
$ orthofinder -f OrthoFinder_protein_folder
```

Convert species marker genes into orthologous genes
```{r}
# load packages
library(tidyr); library(dplyr); library("stringr") 

# load the orthologous groups which is the output from the orthofinder package
orthofinder = read.csv("/Step2_OrthoFinder/Orthogroups_091023_cleaned.tsv", header = TRUE, sep = "\t")

# clean up for each species column
# EX: Arabidopsis
og_ath = orthofinder[,c("Orthogroup", "Arabidopsis")] %>%
         mutate(Arabidopsis = strsplit(as.character(Arabidopsis), ",")) %>% # Split long string in 1 row to multiple rows
         unnest(Arabidopsis) %>%
         mutate(Arabidopsis = str_extract(Arabidopsis, "[^.]+"))  %>% # Extract all character before the first dot
         mutate(across(where(is.character), str_trim)) %>%  # Remove white spaces
         distinct(Arabidopsis, .keep_all = TRUE)  # Remove duplicated rows based on Arabidopsis column

# EX: Rice (Oryza Sativa)
og_oryza = orthofinder[,c("Orthogroup", "Oryza")] %>%
           mutate(Oryza = strsplit(as.character(Oryza), ",")) %>% 
           unnest(Oryza) %>%
           mutate(Oryza = str_extract(Oryza, "[^-]+")) %>% 
           mutate(across(where(is.character), str_trim)) %>%
           distinct(Oryza, .keep_all = TRUE)  %>% 
           mutate(across('Oryza', str_replace, 't', 'g'))

# load marker genes found in step 1, two examples of Arabidospis and Rice (Oryza)
MG_ara = readRDS("/Step1_FindMarkerGenes/MG_120923_ATH_05.RData")
MG_oryza = readRDS("/Step1_FindMarkerGenes/MG_120923_Rice_05.RData")

# Merge marker gene and OG gene table
# EX: Arabidopsis 
Ath_MG_OG = merge(MG_ara, og_ath, by.x = "gene", by.y = "Arabidopsis") %>% 
            arrange( cluster, desc(avg_log2FC)) %>% 
            dplyr:::select("gene", "Orthogroup", "cluster", "avg_log2FC") %>% 
            group_by(cluster) %>% 
            top_n(200)

# Oryza #############################
Oryza_MG_OG = merge(MG_oryza, og_oryza, by.x = "gene", by.y = "Oryza") %>% 
              arrange( cluster, desc(avg_log2FC)) %>% 
              dplyr:::select("gene", "Orthogroup", "cluster", "avg_log2FC") %>% 
              group_by(cluster) %>% 
              top_n(200)

```

#### Step 3, Find conserved orthologous marker genes (OMGs) between species:
Present a comprehensive pairwise comparison of orthologous marker genes for cell type clusters between two plant species, with the aim of generating lists of conserved orthologous marker genes that could be used to identify cell types in other plant species.  </br>

An example of Arabidopsis and Rice
```{r}
# Function create table of common OG gene names between 2 species
count_com_OG = function(Species1, Species2){
  clusters_S1 = unique(Species1$cluster); clusters_S2 = unique(Species2$cluster)
  two_plants <- matrix(nrow=length(clusters_S1), ncol=length(clusters_S2))
  for(i in 1:length(clusters_S1)){
    for(j in 1:length(clusters_S2)){
      list_overlap = intersect(Species1[Species1$cluster == clusters_S1[i],]$Orthogroup, Species2[Species2$cluster == clusters_S2[j],]$Orthogroup)
      num_overlap = length(list_overlap)
    
      two_plants[i,j] = num_overlap
    }
  }
  rownames(two_plants) = unique(Species1$cluster)
  colnames(two_plants) = unique(Species2$cluster) 
  return(two_plants[order(rownames(two_plants)), order(colnames(two_plants))])
}

# load packages
library(ggplot2)
library(reshape)
library(hablar)

# load function and compare Arabidopsis and Rice
df = as.data.frame(count_com_OG(Ath_MG_OG[Ath_MG_OG$cluster %in% c("Xylem", "Cortex",  "Hair1", "Hair2", "Nonhair1", "Nonhair2", "Endodermis", "Phloem","Stele1", "Stele2", "Endocortex"),], Oryza_MG_OG[Oryza_MG_OG$cluster %in%  c("Xylem", "Cortex", "Cortex_like", "Hair1",  "Hair2", "Nonhair",  "Endodermis", "Meristem", "Phloem", "Stele"),]))

df$cell_type = rownames(df) # Add column to make 2 variables when using melt function

# Add a row that sums up the values in all other rows
sumrow = df %>% dplyr:::select(-cell_type) %>% colSums() # sum all numeric row in the dataframe
sum_h = c(sumrow, "sum_h") 
df = rbind(df, sum_h) # After merging two data frames, datatype in dataframe will be changed into character
df = df %>% retype() # %>% as.data.frame(df) # function retype will change the data into the correct type

# Add a column to sum up all values in other columns
df$sum_v = df %>% dplyr:::select(-cell_type) %>% rowSums() # sum all numeric columns in the dataframe
rownames(df) <- df$cell_type  


p_value_dataframe = df[1: (nrow(df) - 1), 1: (ncol(df) -2)]
conclusionTable_0.01 = p_value_dataframe # Reject the null hypothesis if p-value < 0.01
for (i in rownames(p_value_dataframe)){
  for (j in colnames(p_value_dataframe)){
    frame = df[c(i, "sum_h"), c(j, "sum_v")]
    frame[2,2] = frame[2,2] - frame[1,2] - frame[2,1] + frame[1,1]
    frame[1,2] = frame[1,2] - frame[1,1]
    frame[2,1] = frame[2,1] - frame[1,1]
    p_value_dataframe[i,j] = fisher.test(frame, alternative = "greater")$p.value
  }
}
# calculate FDR to correct for multiple tests
adjusted_pvalue_dataframe <- p_value_dataframe
adjusted_pvalue_dataframe[] = p.adjust(unlist(p_value_dataframe), method = "BH") 
conclusionTable_0.01 = ifelse(adjusted_pvalue_dataframe < 0.01, "Reject", "Fail")


table_count = as.matrix(count_com_OG(Ath_MG_OG[Ath_MG_OG$cluster %in% c("Xylem", "Cortex",  "Hair1", "Hair2", "Nonhair1", "Nonhair2", "Endodermis", "Phloem","Stele1", "Stele2", "Endocortex"),], Oryza_MG_OG[Oryza_MG_OG$cluster %in%  c("Xylem", "Cortex", "Cortex_like",  "Hair1",  "Hair2", "Nonhair",  "Endodermis", "Meristem", "Phloem", "Stele"),]))

df1 = melt(table_count)
df1$test = melt(conclusionTable_0.01)$value

# Create the plot with ggplot
plot_ATH_Rice = ggplot(df1, aes(x = factor(X1, levels = c("Xylem", "Cortex",  "Hair1", "Hair2", "Nonhair1", "Nonhair2", "Endodermis", "Phloem","Stele1", "Stele2", "Endocortex")), y = factor(X2, levels = c("Xylem", "Cortex", "Cortex_like",  "Hair1",  "Hair2", "Nonhair",  "Endodermis", "Meristem", "Phloem", "Stele")), fill = value)) +
       geom_tile() + 
       geom_text(aes(label = value), color = "black") +
       scale_fill_gradient(low = "white", high = "darkgreen") +
       geom_rect(data = subset(df1, test == "Reject"), 
            aes(xmin = as.numeric(factor(X1, levels = c("Xylem", "Cortex",  "Hair1", "Hair2", "Nonhair1", "Nonhair2", "Endodermis", "Phloem","Stele1", "Stele2", "Endocortex"))) - 0.5, xmax = as.numeric(factor(X1, levels = c("Xylem", "Cortex",  "Hair1", "Hair2", "Nonhair1", "Nonhair2", "Endodermis", "Phloem","Stele1", "Stele2", "Endocortex"))) + 0.5, 
                ymin = as.numeric(factor(X2, levels = c("Xylem", "Cortex", "Cortex_like", "Hair1",  "Hair2", "Nonhair",  "Endodermis", "Meristem", "Phloem", "Stele"))) - 0.5, ymax = as.numeric(factor(X2, levels = c("Xylem", "Cortex", "Cortex_like", "Hair1",  "Hair2", "Nonhair",  "Endodermis", "Meristem", "Phloem", "Stele"))) + 0.5), 
            fill = NA, color = "red", size = 1) +
       theme(legend.position="right", 
                  axis.line.x = element_line(size = 1),
                  axis.line.y = element_line(size = 1),
                  axis.text.x = element_text(color="black", face = "bold",size=13, angle=90, vjust = 0.4, hjust = 1), 
                  axis.text.y = element_text(color="black", face = "bold",size=13, angle=0, vjust = 0.5, hjust = 1),
                  axis.ticks.length = unit(0.1,"cm"),
                  axis.ticks = element_line(size = 1),
                  legend.title = element_text(size=10, face= "bold"),
                  legend.text = element_text(size=8),
                  plot.title = element_text(size=14, face= "bold", colour= "blue", hjust = 0.5),
                  axis.title.x = element_text(size=18, face="bold", colour = "darkgreen", vjust = 1),    
                  axis.title.y = element_text(size=18, face="bold", colour = "darkgreen")) +
  labs(x = "Arabidopsis", y = "Rice")

plot_ATH_Rice
```
![Test](./Figures/030823_Seurat_ATH_Rice_count2.jpg)

Detailed information can be found in the R markdown in step 3 folder

#### Step 4, Predict cell types in a query species:
Predict cell types in Tomato root and shoot.

```{r}
# Clean orthologous marker genes
orthofinder = read.csv("/Step2_OrthoFinder/Orthogroups_091023_cleaned.tsv", header = TRUE, sep = "\t")
# EX: Arabidopsis root
og_ath = orthofinder[,c("Orthogroup", "Arabidopsis")] %>%
  mutate(Arabidopsis = strsplit(as.character(Arabidopsis), ",")) %>% 
  unnest(Arabidopsis) %>%
  mutate(Arabidopsis = str_extract(Arabidopsis, "[^.]+"))  %>% 
  mutate(across(where(is.character), str_trim)) %>% 
  distinct(Arabidopsis, .keep_all = TRUE)  
  #remove_rownames %>% column_to_rownames(var="Arabidopsis")

# EX: Tomato root
og_tom = orthofinder[,c("Orthogroup", "Solanum_lycopersicum")] %>%
  mutate(Solanum_lycopersicum = strsplit(as.character(Solanum_lycopersicum), ",")) %>% 
  unnest(Solanum_lycopersicum) %>%
  mutate(Solanum_lycopersicum = str_extract(Solanum_lycopersicum, "[^.]+"))  %>% 
  mutate(across(where(is.character), str_trim)) %>% 
  distinct(Solanum_lycopersicum, .keep_all = TRUE) 
  #remove_rownames %>% column_to_rownames(var="Solanum") 

# load marker genes
MG_ara = readRDS("/Step1_FindMarkerGenes/MG_120923_ATH_05.RData")
MG_tom = readRDS("/Step1_FindMarkerGenes/MG_082323_Tom_05_Default15clusters.RData")

# merge marker genes and orthologous file
Ath_MG_OG = merge(MG_ara, og_ath, by.x = "gene", by.y = "Arabidopsis") %>% 
            arrange( cluster, desc(avg_log2FC)) %>% 
            dplyr::select("gene", "Orthogroup", "cluster", "avg_log2FC") %>% 
            group_by(cluster) %>% 
            top_n(200)

Tomato_MG_OG = merge(MG_tom, og_tom, by.x = "gene", by.y = "Solanum_lycopersicum") %>% 
               arrange( cluster, desc(avg_log2FC)) %>% 
               dplyr::select("gene", "Orthogroup", "cluster", "avg_log2FC") %>% 
               group_by(cluster) %>% 
               top_n(200)

# Function create table of common OG gene names between 2 species
count_com_OG = function(Species1, Species2){
  clusters_S1 = unique(Species1$cluster); clusters_S2 = unique(Species2$cluster)
  two_plants <- matrix(nrow=length(clusters_S1), ncol=length(clusters_S2))
  for(i in 1:length(clusters_S1)){
    for(j in 1:length(clusters_S2)){
      list_overlap = intersect(Species1[Species1$cluster == clusters_S1[i],]$Orthogroup, Species2[Species2$cluster == clusters_S2[j],]$Orthogroup)
      num_overlap = length(list_overlap)
    
      two_plants[i,j] = num_overlap
    }
  }
  rownames(two_plants) = unique(Species1$cluster)
  colnames(two_plants) = unique(Species2$cluster) 
  return(two_plants[order(rownames(two_plants)), order(colnames(two_plants))])
}

# perform pairwise comparison
df = as.data.frame(count_com_OG(Ath_MG_OG[Ath_MG_OG$cluster %in% c("Nonhair1", "Hair1", "Hair2", "Cortex", "Endodermis", "Endocortex",  "Xylem", "Phloem", "Stele2","MZ1", "MZ2"),], Tomato_MG_OG))
df$cell_type = rownames(df) # Add column to make 2 variables when using melt function


# Add a row that sums up the values in all other rows
sumrow = df %>% dplyr::select(-cell_type) %>% colSums() # sum all numeric row in the dataframe
sum_h = c(sumrow, "sum_h") 
df = rbind(df, sum_h) # After merging two data frames, datatype in dataframe will be changed into character
library(hablar)
df = df %>% retype() #%>% as.data.frame(df) # This library and function retype will change the data into the correct type

# Add a column to sum up all values in other columns
df$sum_v = df %>% dplyr::select(-cell_type) %>% rowSums() # sum all numeric columns in the dataframe
rownames(df) <- df$cell_type  

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
# Calculate the FDR to correct for multiple tests
adjusted_pvalue_dataframe <- p_value_dataframe
adjusted_pvalue_dataframe[] = p.adjust(unlist(p_value_dataframe), method = "BH")
conclusionTable_0.01 = ifelse(adjusted_pvalue_dataframe < 0.01, "Reject", "Fail")

#-----
table_count = as.matrix(count_com_OG(Ath_MG_OG[Ath_MG_OG$cluster %in% c("Nonhair1", "Hair1", "Hair2", "Cortex", "Endodermis", "Endocortex",  "Xylem", "Phloem", "Stele2","MZ1", "MZ2"),], Tomato_MG_OG))
df1 = melt(table_count)
df1$test = melt(conclusionTable_0.01)$value


# Generate the plot using ggplot
plot_ATH_Tom = ggplot(df1, aes(x = factor(X1, levels = c("Nonhair1", "Hair1", "Hair2", "Cortex", "Endodermis", "Endocortex",  "Xylem", "Phloem", "Stele2","MZ1", "MZ2")), y = factor(X2), fill = value)) +
       geom_tile() + 
       geom_text(aes(label = value), color = "black") +
       scale_fill_gradient(low = "white", high = "darkgreen") +
       geom_rect(data = subset(df1, test == "Reject"), 
            aes(xmin = as.numeric(factor(X1, levels = c("Nonhair1", "Hair1", "Hair2", "Cortex", "Endodermis", "Endocortex",  "Xylem", "Phloem", "Stele2","MZ1", "MZ2"))) - 0.5, xmax = as.numeric(factor(X1, levels = c("Nonhair1", "Hair1", "Hair2", "Cortex", "Endodermis", "Endocortex",  "Xylem", "Phloem", "Stele2","MZ1", "MZ2"))) + 0.5, 
                ymin = as.numeric(X2) + 0.5 , ymax = as.numeric(X2) + 1.5), 
            fill = NA, color = "red", size = 1) +
       theme(legend.position="right", 
                  axis.line.x = element_line(size = 1),
                  axis.line.y = element_line(size = 1),
                  axis.text.x = element_text(color="black", face = "bold",size=12, angle=90, vjust = 0.4, hjust = 1), 
                  axis.text.y = element_text(color="black", face = "bold",size=12, angle=0, vjust = 0.5, hjust = 1),
                  axis.ticks.length = unit(0.1,"cm"),
                  axis.ticks = element_line(size = 1),
                  legend.title = element_text(size=10, face= "bold"),
                  legend.text = element_text(size=8),
                  plot.title = element_text(size=14, face= "bold", colour= "blue", hjust = 0.5),
                  axis.title.x = element_text(size=17, face="bold", colour = "darkgreen", vjust = 1),    
                  axis.title.y = element_text(size=17, face="bold", colour = "darkgreen")) +
  labs(x = "Arabidopsis", y = "Tomato")

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
