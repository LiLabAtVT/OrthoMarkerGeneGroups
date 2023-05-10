# Compare the conserved OMGs among three plant species

### Introduction:
These scripts produce Venn diagrams that compare the shared OMGs among three plant species, with the purpose of generating lists of conserved OMGs that can aid in identifying cell types in other plant species.

### Processing:
* Use the output generated from step 1 and step 2 as input for this step. </br>
* Select the top 200 marker genes identified in step 1:
```R
Oryza_MG_OG = merge(Marker_genes, Orthologous_group, by.x = "gene", by.y = "Species") %>% 
          arrange( cluster, desc(avg_log2FC)) %>% 
          select("gene", "Orthogroup", "cluster", "avg_log2FC") %>% 
          group_by(cluster) %>% 
          top_n(200)
```
* Generate a Venn diagram to determine the number of shared OMGs across the same cell type in three plant species:
```R
# List of items
cortex <- list(A = unique(Oryza_MG_OG[Oryza_MG_OG$cluster == "Cortex",]$Orthogroup), 
               B = unique(Maize_MG_OG[Maize_MG_OG$cluster == "Cortex",]$Orthogroup),
               C = unique(Ath_MG_OG[Ath_MG_OG$cluster == "Cortex",]$Orthogroup)
               )
venn = Venn(cortex)
data = process_data(venn, shape_id == "301f")

cortex_plot = ggplot() +
  ggtitle("Cotex") +  
  geom_sf(aes(fill = count), data = venn_region(data)) +
  geom_sf(color = c("purple2","tomato3","orange2"), size = 1, data = venn_setedge(data), show.legend = FALSE) +
  geom_sf_text(label = c("Oryza", "Zeamays", "Arabidopsis"), 
               size = 6, 
               color = c("purple2","tomato3","orange2"), 
               fontface = "bold", 
               data = venn_setlabel(data)) +
  geom_sf_label(aes(label = count), data = venn_region(data), size = 5) +
  scale_fill_gradient(low="white",high = "green3") +
  theme_void() +
  theme(legend.position="right", 
                  legend.title = element_text(size=12),
                  legend.text = element_text(size=10),
                  plot.title = element_text(size=18, face= "bold", colour= "black", hjust = 0.5)
        )  
```

### Result:
The file named "ComMG_Seurat.xlsx" contains the list of common OMGs among three plant species, as well as their corresponding marker genes for each species.
