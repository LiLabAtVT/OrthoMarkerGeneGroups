# Extract marker genes from a tissue from the file MG_example.csv
MarkerGenes = function(species_MG_frame, tissue_MG_frame){
  df = MGs[(MGs$species == species_MG_frame) & (MGs$tissue == tissue_MG_frame),]
  return(df)
}