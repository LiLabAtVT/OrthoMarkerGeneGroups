# OMG Browser:

## **Introduction**

This Shiny app is designed to visualize gene expression data across three different plant species (Arabidopsis thaliana, Maize, and Rice). The app is called "OMG Browser" and allows users to select a species of interest and explore gene expression patterns across different clusters of cells. The app also provides the ability to search for a specific gene and display the expression pattern of that gene across all clusters.

## **Required Packages**

The app uses the following R packages: shiny, ggplot2, and dplyr. Please make sure to have these packages installed before running the app.

## **Data**

The app loads data from RDS files containing gene expression data for each species, as well as clustering data and expression matrices for each species. The files must be located in the "omgrds" and "rdsfiles" directories, respectively. You can refer to preprocessing steps to know more about data setup.

## **Running the App**

To run the app, please ensure that all the necessary packages are installed and the data files are in the correct directories. Then, simply run the following code:

```r
library(shiny)
runApp("path/to/app")
```

Replace "path/to/app" with the actual file path where the app files are located.

## **App Interface**

Plots-Seurat: allows users to select a species of interest, select a cluster, and visualize gene expression across all clusters. Users can either select a gene from a dropdown or search for a specific gene using a text input box.

### **Sidebar**

The sidebar consists of the following elements:

- Select input type: users can select whether to input a gene using a dropdown or text input box.
- Select a Species: users can select the species of interest (Arabidopsis thaliana, Maize, or Rice).
- Cluster Dropdown: once a species is selected, the available clusters are displayed in a dropdown menu.
- Gene Dropdown: once a cluster is selected, the available genes are displayed in a dropdown menu.
- Text Input: if "Text Input" is selected, users can enter the gene ID of interest using the text input box.
- Matched OMG: displays the Orthogroup of the selected gene.

### **Main Panel**

The main panel consists of three columns of plots, each containing a UMAP plot and a gene expression plot. The plots display gene expression patterns across different clusters of cells for the selected gene and species.

## UI **Design**

The app interface is designed to be user-friendly and easy to navigate. The app uses a dark color scheme with contrasting text to make it easy to read. The navbar at the top of the app allows users to switch between tabs. The sidebar provides users with the ability to select a species, cluster, and gene of interest, as well as the option to input a gene ID manually. The main panel displays the plots and allows users to explore gene expression patterns. Finally, the footer contains information about the creators of the app and their affiliation.

# Pre-Processing for Data creation:

## Steps:

There are three steps involved in the preprocessing step for this application. These steps output files which are necessary in the final shiny application.

1. **Pre-processing raw rds files for all species:**
    
    The first step R performs down-sampling on raw RDS files for three different species (ATH, rice, maize) and plots the Uniform Manifold Approximation and Projection (UMAP) of the down-sampled data. The following steps are performed for each species:
    
    - Load raw RDS data.
    - Normalize the data before processing.
    - Find top 3000 genes that cover the most variance.
    - Subset Seurat object using top 1000 variable features.
    - Subset 300 cells.
    - Save down-sampled Seurat object as RDS.
    - Plot UMAP for verification of overall structure.
    
2. O**ptimization**
    
    The second step is used to extract cluster and UMAP information into a single RDS file to save time and execution time. Additionally, expression matrix information is also extracted into separate RDS files. This step is repeated for all species. Users need to update the path in line 16(for variable rds_downsampled) to run this script for different species.
    
3. **OMG Mapping file Preprocessing**
    
    The purpose of this third step is to extract information from the raw OMG mapping files and save it in separate RDS files for faster execution of a Shiny App. 
    
    - The script loads the necessary packages, reads in the three sheets from an Excel file.
    - Sorts them by a column, renames columns for each dataframe, and saves them as RDS files.
