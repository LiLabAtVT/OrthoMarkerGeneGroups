# This script is just to extract OMG information from the raw OMG mapping
# files and save it in seperate RDS files for faster execution of shiny App.

# Load necessary packages
library(dplyr)

# uncomment these for converting the excel sheets into RDS files

# # Read in the three sheets and sort by V4 column
ath_seurat <- read_excel("/Users/saipavanbathala/Documents/Li_lab/shinyOMG/Raw_OMG_Mapping_files/OMG.xlsx", sheet = "ATH_SHAP") %>% arrange(desc(V4))
maize_seurat <- read_excel("/Users/saipavanbathala/Documents/Li_lab/shinyOMG/Raw_OMG_Mapping_files/OMG.xlsx", sheet = "Maize_SHAP") %>% arrange(desc(V4))
rice_seurat <- read_excel("/Users/saipavanbathala/Documents/Li_lab/shinyOMG/Raw_OMG_Mapping_files/OMG.xlsx", sheet = "Rice_SHAP") %>% arrange(desc(V4))
#
# # Rename columns for each dataframe
ath_seurat <- rename(ath_seurat, ath_gene = V1)
maize_seurat <- rename(maize_seurat, maize_gene = V1)
rice_seurat <- rename(rice_seurat, rice_gene = V1)
#
saveRDS(ath_seurat, file = "/Users/saipavanbathala/Documents/Li_lab/shinyOMG/30_OMG_Processing/omg_rds/ath_shap.rds")
saveRDS(maize_seurat, file = "/Users/saipavanbathala/Documents/Li_lab/shinyOMG/30_OMG_Processing/omg_rds/maize_shap.rds")
saveRDS(rice_seurat, file = "/Users/saipavanbathala/Documents/Li_lab/shinyOMG/30_OMG_Processing/omg_rds/rice_shap.rds")

#reading the rds files
# readRDS("/Users/saipavanbathala/Documents/Li_lab/shinyOMG/omg_rds/ath_seurat.rds")
# readRDS("/Users/saipavanbathala/Documents/Li_lab/shinyOMG/omg_rds/maize_seurat.rds")
# readRDS("/Users/saipavanbathala/Documents/Li_lab/shinyOMG/omg_rds/rice_seurat.rds")