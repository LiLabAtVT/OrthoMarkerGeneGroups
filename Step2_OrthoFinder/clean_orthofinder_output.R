clean_orthogroups_gene_names <- function(orthogroups, species_name) {
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