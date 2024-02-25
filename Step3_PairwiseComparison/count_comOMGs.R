# Function create table of common OMGs between 2 species
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

#####-----#####-----#####-----#####-----#####-----#####

# This function tests which common OMGs is significant
test_significant = function(Species1_OMG, Species2_OMG, FDR_threshold){
  
  df = as.data.frame(count_com_OG(Species1_OMG, Species2_OMG))
  df$cell_type = rownames(df) # Add column to make 2 variables when using melt function

  # Add a row that sums up the values in all other rows
  sumrow = df %>% dplyr:::select(-cell_type) %>% colSums() # sum all numeric row in the dataframe
  sum_h = c(sumrow, "sum_h") 
  df = rbind(df, sum_h) # After merging two data frames, datatype in dataframe will be changed into character
  df = df %>% retype() # %>% as.data.frame(df) # function retype will change the data into the correct type

  # Add a column to sum up all values in other columns
  df$sum_v = df %>% dplyr:::select(-cell_type) %>% rowSums() # sum all numeric columns in the dataframe
  rownames(df) <- df$cell_type  

  # Calculate the p-value from Fisher exact test
  p_value_dataframe = df[1: (nrow(df) - 1), 1: (ncol(df) -2)]
  conclusionTable = p_value_dataframe # Reject the null hypothesis if p-value < 0.01
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
  conclusionTable = ifelse(adjusted_pvalue_dataframe < FDR_threshold, "Reject", "Fail")
  # merge test column to the count common OMGs table
  table_count = as.matrix(count_com_OG(Species1_OMG, Species2_OMG))
  df1 = melt(table_count) 
  #adjusted_pvalue_dataframe = tibble::rownames_to_column(adjusted_pvalue_dataframe, "X1")
  df1$FDR = melt(adjusted_pvalue_dataframe)$value
  df1$test = melt(conclusionTable)$value 
  
  return(df1)
}

# Extract data from the heatmap
extract_table <- function(Species1_OMG, Species2_OMG, FDR_threshold){
  df = test_significant(Species1_OMG, Species2_OMG, FDR_threshold)

##
  Species_x_cluster_numOMG = Species1_OMG %>%
                       group_by(cluster) %>%
                       summarize(Species_x_numOMG = n_distinct(Orthogroup))

  Species_y_cluster_numOMG = Species2_OMG %>%
                        group_by(cluster) %>%
                        summarize(Species_y_numOMG = n_distinct(Orthogroup))
##
  extractTable = merge(merge(df, Species_x_cluster_numOMG, by.x = "X1", by.y = "cluster"), Species_y_cluster_numOMG, by.x = "X2", by.y = "cluster") 
  names(extractTable)[1] = "Species_y_clusters"
  names(extractTable)[2] = "Species_x_clusters"
  names(extractTable)[3] = "comOMG"
  extractTable <- extractTable[, c("Species_x_clusters", "Species_x_numOMG", "Species_y_clusters", "Species_y_numOMG", "comOMG", "FDR", "test")] %>% arrange(FDR)
  return(extractTable)
}

# Extract OMG genes
extract_gene <- function(Species1_OMG, Species2_OMG, Species1_cluster, Species2_cluster){
  Species_1_genes = Species1_OMG[Species1_OMG$cluster == Species1_cluster,]
  Species_2_genes = Species2_OMG[Species2_OMG$cluster == Species2_cluster,]
  common_OMG = merge(Species_1_genes, Species_2_genes, by = "Orthogroup")$Orthogroup
  Species_1_genes_wCommonOMG = Species_1_genes[Species_1_genes$Orthogroup %in% common_OMG, ]
  Species_1_genes_wCommonOMG$SpeciesName = "on_X_axis"
  Species_2_genes_wCommonOMG = Species_2_genes[Species_2_genes$Orthogroup %in% common_OMG, ]
  Species_2_genes_wCommonOMG$SpeciesName = "on_Y_axis"
  df = rbind(Species_1_genes_wCommonOMG, Species_2_genes_wCommonOMG)
  return(df) # Arrange by the avg_log2FC values
}

# Generate the plot for comparison
generate_plot_comparison <- function(comOMGs_test, species_x, species_y) {
  plot_comparison <- ggplot(comOMGs_test, aes(x = factor(X1), y = factor(X2), fill = value)) +
    geom_tile() +
    geom_text(aes(label = value), color = "black") +
    scale_fill_gradient(low = "white", high = "darkgreen") +
    geom_rect(data = subset(comOMGs_test, test == "Reject"),
              aes(xmin = as.numeric(X1) - 0.5, xmax = as.numeric(X1) + 0.5,
                  ymin = as.numeric(X2) - 0.5, ymax = as.numeric(X2) + 0.5),
              fill = NA, color = "red", size = 1) +
        theme(legend.position = "right",
          axis.line.x = element_line(size = 1),
          axis.line.y = element_line(size = 1),
          axis.text.x = element_text(color = "black", face = "bold", size = 13, angle = 90, vjust = 0.4, hjust = 1),
          axis.text.y = element_text(color = "black", face = "bold", size = 13, angle = 0, vjust = 0.5, hjust = 1),
          axis.ticks.length = unit(0.1, "cm"),
          axis.ticks = element_line(size = 1),
          legend.title = element_text(size = 10, face = "bold"),
          legend.text = element_text(size = 8),
          plot.title = element_text(size = 14, face = "bold", colour = "blue", hjust = 0.5),
          axis.title.x = element_text(size = 18, face = "bold", colour = "darkgreen", vjust = 1),
          axis.title.y = element_text(size = 18, face = "bold", colour = "darkgreen")) +
    labs(x = species_x, y = species_y)

  return(plot_comparison)
}
