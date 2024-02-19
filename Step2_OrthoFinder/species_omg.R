# Merge marker genes and OG genes table:
Species_OMGs = function(Marker_genes, orthogroup, species){
  Species_MG_OG = merge(Marker_genes, orthogroup, by.x = "gene", by.y = species) %>% 
              arrange( cluster, desc(avg_log2FC)) %>% 
              dplyr:::select("gene", "Orthogroup", "cluster", "avg_log2FC") %>% 
              group_by(cluster) %>% 
              top_n(200)
  return(Species_MG_OG)
}