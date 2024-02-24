clean_OG <- function(orthogroups, species_name) {
  # Check if the species_name column exists in the orthogroups DataFrame
  if (!species_name %in% colnames(orthogroups)) {
    stop("Species name is not in the orthogroups table")
  } else {
    # Selecting the Orthogroup and species_name columns
    orthogroups <- orthogroups[, c("Orthogroup", species_name)] %>%
      mutate_at(vars(species_name), list(~strsplit(as.character(.), ","))) %>%
      unnest(cols = species_name) %>%
      mutate(across(where(is.character), str_trim)) %>% # Remove white spaces
      dplyr:::distinct(!!sym(species_name), .keep_all = TRUE)
    
    return(orthogroups)
  }
}

#clean_orthogroups_gene_names(orthofinder, "Arabidopsis")

# Merge marker genes and OG genes table:
merge_MG_OG = function(Marker_genes, orthogroup, species){
  Species_MG_OG = merge(Marker_genes, orthogroup, by.x = "gene", by.y = species) %>% 
              arrange( cluster, desc(avg_log2FC)) %>% 
              dplyr:::select("gene", "Orthogroup", "cluster", "avg_log2FC") %>% 
              group_by(cluster) %>% 
              top_n(200)
  return(Species_MG_OG)
}