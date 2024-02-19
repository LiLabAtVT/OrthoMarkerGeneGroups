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
test_significant = function(Species1, Species2, FDR_threshold){
  
  df = as.data.frame(count_com_OG(Species1, Species2))
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
  table_count = as.matrix(count_com_OG(Species1, Species2))
  df1 = melt(table_count) 
  df1$test = melt(conclusionTable)$value 
  
  return(df1)
}

