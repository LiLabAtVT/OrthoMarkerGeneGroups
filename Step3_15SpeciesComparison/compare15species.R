
library(reshape2); library(tidyverse)

## Function create table that counts common OG gene names between 2 species
## Save the numbers to the table
count_com_OG_number = function(Species1, Species2){
  clusters_S1 = unique(Species1$clusterName); clusters_S2 = unique(Species2$clusterName)
  two_plants <- matrix(nrow=length(clusters_S1), ncol=length(clusters_S2))
  for(i in 1:length(clusters_S1)){
    for(j in 1:length(clusters_S2)){
      list_overlap = intersect(Species1[Species1$clusterName == clusters_S1[i],]$Orthogroup, Species2[Species2$clusterName == clusters_S2[j],]$Orthogroup)
      num_overlap = length(list_overlap)
      
      two_plants[i,j] = num_overlap
    }
  }
  rownames(two_plants) = unique(Species1$clusterName)
  colnames(two_plants) = unique(Species2$clusterName) 
  return(two_plants[order(rownames(two_plants)), order(colnames(two_plants))])
}

## Function create table of common OG gene names between 2 species
## Save the common OG name to the table
count_com_OG = function(Species1, Species2){
  clusters_S1 = unique(Species1$clusterName); clusters_S2 = unique(Species2$clusterName)
  two_plants <- data.frame(matrix(vector("list", length(clusters_S1) * length(clusters_S2)), 
                                  nrow = length(clusters_S1), 
                                  ncol = length(clusters_S2)))
  for(i in 1:length(clusters_S1)){
    for(j in 1:length(clusters_S2)){
      list_overlap = intersect(Species1[Species1$clusterName == clusters_S1[i],]$Orthogroup, Species2[Species2$clusterName == clusters_S2[j],]$Orthogroup)
      two_plants[[i,j]] = list(list_overlap)
      
    }
  }
  rownames(two_plants) = unique(Species1$clusterName)
  colnames(two_plants) = unique(Species2$clusterName)
  return(two_plants[order(rownames(two_plants)), order(colnames(two_plants))])
}

Heatmap_with_count_comOMG = function(reference_data, query_data){
  # Add a column that combine all common OG in the same row, and add the row that combine all lists in the same column
  df = count_com_OG(reference_data, query_data)
  colnames(df) = gsub(" ", "_",colnames(df))
  
  # Add row that sum for each column
  new_row <- list()
  for (col in names(df)) {
    concatenated_list <- unlist(df[[col]], use.names = FALSE)
    new_row[[col]] <- list(concatenated_list)
  }
  new_row_df <- t(enframe(new_row, name = NULL)) %>% as_tibble() 
  colnames(new_row_df) = colnames(df)
  row.names(new_row_df) <- "new_row"
  df <- rbind(df, new_row_df)
  
  
  # Add column that sum for each row
  new_col <- vector("list", nrow(df))
  for (i in 1:nrow(df)) {
    concatenated_list <- list()
    for (col in names(df)) {
      concatenated_list <- c(concatenated_list, unlist(df[[col]][[i]]))
    }
    new_col[[i]] <- unlist(concatenated_list, recursive = TRUE)
  }
  df$new_col <- new_col
  ########################################
  
  
  
  p_value_dataframe = df[1: (nrow(df) - 1), 1: (ncol(df) -1)] #Create a dataframe p-value for each pair comparison
  
  for (i in rownames(p_value_dataframe)){
    for (j in colnames(p_value_dataframe)){
      frame = df[c(i, "new_row"), c(j, "new_col")]
      # 2 matrices, 1 for common OG in list, 1 count the number of OG
      frame_int = data.frame(Column1 = c(0, 0), Column2 = c(0, 0))
      frame_int[2,2] = length(setdiff(setdiff(setdiff(frame[2,2][[1]], frame[1,2][[1]]), frame[2,1][[1]][[1]]), frame[1,1][[1]][[1]]))
      frame_int[1,2] = length(setdiff(frame[1,2][[1]], frame[1,1][[1]][[1]]))
      frame_int[2,1] = length(setdiff(frame[2,1][[1]][[1]], frame[1,1][[1]][[1]]))
      frame_int[1,1] = length(frame[1,1][[1]][[1]])
      
      p_value_dataframe[i,j] = fisher.test(frame_int, alternative = "greater")$p.value
    }
  }
  adjusted_pvalue_dataframe <- p_value_dataframe #Create another frame for FDR
  adjusted_pvalue_dataframe[] = p.adjust(unlist(p_value_dataframe), method = "BH")
  
  ## Extract table
  frame = melt(count_com_OG_number(reference_data, query_data))
  frame$p_value = melt(adjusted_pvalue_dataframe)$value
  frame$transform  = -log10(frame$p_value)
  return(frame)
}


compare_w_15_species <- function(sample_input_csv){
Ortho = read.delim( "Step3_15SpeciesComparison/ortho_extract.tsv", header = TRUE, sep = "\t")
Marker_Gene = read.delim( "Step3_15SpeciesComparison/MG_all_list.tsv", header = TRUE, sep = "\t")
Marker_Gene$clusterName = paste(Marker_Gene$species, Marker_Gene$tissue, Marker_Gene$clusterName)

query = Marker_Gene[  ((Marker_Gene$species == "arabidopsis_thaliana") & (Marker_Gene$tissue == "Hypocotyl callus")) |
                        ((Marker_Gene$species == "arabidopsis_thaliana") & (Marker_Gene$tissue == "Inflorescence")) |
                        ((Marker_Gene$species == "arabidopsis_thaliana") & (Marker_Gene$tissue == "Leaf")) |
                        ((Marker_Gene$species == "arabidopsis_thaliana") & (Marker_Gene$tissue == "Pollen")) |
                        ((Marker_Gene$species == "arabidopsis_thaliana") & (Marker_Gene$tissue == "Root")) |
                        ((Marker_Gene$species == "arabidopsis_thaliana") & (Marker_Gene$tissue == "Shoot axis apex")) |
                        
                        ((Marker_Gene$species == "oryza_sativa") & (Marker_Gene$tissue == "Inflorescence")) |
                        ((Marker_Gene$species == "oryza_sativa") & (Marker_Gene$tissue == "Leaf")) |
                        ((Marker_Gene$species == "oryza_sativa") & (Marker_Gene$tissue == "Pistil")) |
                        ((Marker_Gene$species == "oryza_sativa") & (Marker_Gene$tissue == "Root")) |
                        ((Marker_Gene$species == "oryza_sativa") & (Marker_Gene$tissue == "Root;Leaf")) |
                        
                        #((Marker_Gene$species == "triticum_aestivum") & (Marker_Gene$tissue == "Root")) |
                        
                        ((Marker_Gene$species == "zea_mays") & (Marker_Gene$tissue == "Ear inflorescence")) |
                        ((Marker_Gene$species == "zea_mays") & (Marker_Gene$tissue == "Leaf")) |
                        ((Marker_Gene$species == "zea_mays") & (Marker_Gene$tissue == "Root")) |
                        ((Marker_Gene$species == "zea_mays") & (Marker_Gene$tissue == "Shoot axis apex")) |
                        
                        ((Marker_Gene$species == "brassica_rapa") & (Marker_Gene$tissue == "Leaf")) |
                        ((Marker_Gene$species == "catharanthus_roseus") & (Marker_Gene$tissue == "Leaf")) |
                        ((Marker_Gene$species == "fragaria_vesca") & (Marker_Gene$tissue == "Leaf")) |
                        ((Marker_Gene$species == "gossypium_bickii") & (Marker_Gene$tissue == "Seed")) |
                        ((Marker_Gene$species == "gossypium_hirsutum") & (Marker_Gene$tissue == "Ovule outer integument")) |
                        ((Marker_Gene$species == "glycine_max") & (Marker_Gene$tissue == "Root; Root nodule")) | 
                        ((Marker_Gene$species == "manihot_esculenta") & (Marker_Gene$tissue == "Tuberous root"))| 
                        ((Marker_Gene$species == "medicago_truncatula") & (Marker_Gene$tissue == "Root"))|
                        ((Marker_Gene$species == "nicotiana_attenuata") & (Marker_Gene$tissue == "Corolla")) |
                        
                        ((Marker_Gene$species == "populus_alba_var_pyramidalis") & (Marker_Gene$tissue == "Stem")) |
                        ((Marker_Gene$species == "populus_alba_x_populus_glandulosa") & (Marker_Gene$tissue == "Stem")) |
                        ((Marker_Gene$species == "solanum_lycopersicum") & (Marker_Gene$tissue == "Root")) |
                        ((Marker_Gene$species == "solanum_lycopersicum") & (Marker_Gene$tissue == "Shoot axis apex")) ,]

## Input
Sample_MG = read.csv(sample_input_csv)
Sample_frame = merge(Sample_MG, Ortho, by.x = "gene", by.y = "MarkerGene", all.x = TRUE) %>% arrange(cluster, desc(avg_log2FC)) %>% group_by(cluster) %>% slice(1:200)
Sample_frame$clusterName = Sample_frame$cluster

###
query_frame = merge(query, Ortho, by.x = "gene", by.y = "MarkerGene", all.x = TRUE) %>% group_by(species) %>% arrange(clusterName, desc(avg_log2FC)) %>% group_by(clusterName) %>% slice(1:200)
frame_output <- Heatmap_with_count_comOMG(query_frame, Sample_frame)
colnames(frame_output) <- c("reference_cell_types", "query_clusters", "number_OMGs", "p_value", "negative_log10_p_value")
return(frame_output)
}
