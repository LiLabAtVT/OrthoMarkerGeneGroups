# Load required packages
library(shiny)
library(ggplot2)
library(dplyr)

#Reading RDS files of three species
ath_selected <- readRDS("omgrds/ath_seurat.rds")
maize_selected <- readRDS("omgrds/maize_seurat.rds")
rice_selected <- readRDS("omgrds/rice_seurat.rds")

# Read clustering data and expression matrix of ATH
ath_expr_mat <- readRDS("rdsfiles/ath_expr_data.rds")
ath_cluster_info <- readRDS("rdsfiles/ath_cluster_info.rds")

# Read clustering data and expression matrix of Maize
maize_expr_mat <- readRDS("rdsfiles/maize_expr_data.rds")
maize_cluster_info <- readRDS("rdsfiles/maize_cluster_info.rds")

# Read clustering data and expression matrix of Rice
rice_expr_mat <- readRDS("rdsfiles/rice_expr_data.rds")
rice_cluster_info <- readRDS("rdsfiles/rice_cluster_info.rds")


#Creating the UI part of the app
ui <- fluidPage(
  tags$head(
    tags$style(
      HTML("
        .navbar-brand {
          font-weight: bold;
        }
        .nav-link {
          color: white !important;
        }
        .nav-link.active {
          font-weight: bold;
          background-color: white !important;
          color: #4285F4 !important;
        }
        .footer {
          background-color: #343a40;
          color: white;
          position: fixed;
          left: 0;
          bottom: 0;
          width: 100%;
          padding: 10px;
          text-align: center;
        }
        .footer a {
          color: white;
        }
      ")
    )
  ),
  tags$div(
    class = "jumbotron text-center",
    style = "background-color: #333; color: white; padding: 8px;margin: 11px;width: 100%;",
    tags$h1("OMG Browser", style = "font-size: 25px;"),
  ),
  navbarPage(
    title = "",
    theme = "bootstrap",
    tabPanel(
      title = "Plots Across Species",
      id = "plots",
      
      sidebarLayout(
        sidebarPanel(
          width = 2,
          fluidRow(
            column(width = 12, 
                   selectInput("algorithm", "Select a Algorithm:", choices = c("Seurat", "SVM", "SHAP"), selected = "Seurat")
            ),
            column(
              width = 12,
              div(
                #style = "display: none;",
                radioButtons(
                  "gene_input_type", 
                  label = "Select input type", 
                  choices = c("Dropdown", "Text Input"),
                  selected = "Text Input",
                  inline = TRUE
                )
              )
            ),
            

            column(width = 12, 
                   selectInput("species", "Select a Species:", choices = c("ATH", "Maize", "Rice"), selected = "ATH")
            ),
            column(width = 12, 
                   uiOutput("cluster_dropdown")            ),
            column(width = 12, 
                   uiOutput("gene_dropdown")
            ),
            column(width = 12, 
                  div(style = "text-align: center;", "(or)")
            ),
            column(width = 12, br()),
            column(width = 12, 
                   textInput("gene_text", "Enter gene id:", value = "AT4G35350"),
            ),
            column(width = 12, 
                   tags$strong("Matched OMG for this gene id:")),
            column(width = 12, 
                   verbatimTextOutput("matched_omg")
            )
          )
        ),
        mainPanel(
          width = 10,
          fluidRow(
            column(width = 4, plotOutput("umap1"), plotOutput("exp1")),
            column(width = 4, plotOutput("umap2"), plotOutput("exp2")),
            column(width = 4, plotOutput("umap3"), plotOutput("exp3"))
          )
        )
      )
    ),
    tabPanel(
      title = "Additonal Information",
      id = "plots_svm"
    )
  ),
  tags$div(
    class = "footer",
    HTML(paste0("<span style='float: left;'>Made in <a href='https://lilabatvt.github.io/research/' target='_blank'>Li Lab of Applied Machine Learning in Genomics and Phenomics, Virginia Tech</a></span> ",
                "<span style='float: right;'>By: <a href='mailto:tnchau@vt.edu'>Tran Chau</a> and <a href='mailto:bspavan25@vt.edu'>Sai Pavan Bathala</a></span>"))
  )
)


server <- function(input, output, session) {
  
  # Create reactive expression for selected cluster
  output$cluster_dropdown <- renderUI({
    
    # Get cluster info based on selected species
    cluster_info <- switch(input$species,
                           "ATH" = ath_cluster_info,
                           "Maize" = maize_cluster_info,
                           "Rice" = rice_cluster_info)
    
    # Get available clusters for selected species
    available_clusters <- unique(cluster_info$Cluster)

    selectInput("cluster_dropdown", "Select a Cluster:", choices = available_clusters)
  })
  
  #drop down for genes
  output$gene_dropdown <- renderUI({

    # Get cluster info based on dropdown
    dummy <-input$cluster_dropdown

    # Get cluster info based on selected species
    cluster_info <- switch(input$species,
                           "ATH" = ath_cluster_info,
                           "Maize" = maize_cluster_info,
                           "Rice" = rice_cluster_info)

    # Get cluster info based on selected species
    expr_mat_info <- switch(input$species,
                            "ATH" = ath_expr_mat,
                            "Maize" = maize_expr_mat,
                            "Rice" = rice_expr_mat)


    if (is.null(input$cluster_dropdown)) {
    cell_ids <- cluster_info %>%
      filter(Cluster == "Xylem") %>%
      pull(Cell)
    }
    else{
    cell_ids <- cluster_info %>%
      filter(Cluster == input$cluster_dropdown) %>%
      pull(Cell)
    }

    # Select cells with expression greater than zero
    expr_mat_filtered <- expr_mat_info[, cell_ids][apply(expr_mat_info[, cell_ids], 1, function(x) any(x > 0)), ]

    # Select genes with expression greater than zero in the selected cells
    available_genes <- rownames(expr_mat_filtered)[apply(expr_mat_filtered, 2, function(x) any(x > 0))]
    
    # Subset expression matrix for available genes
    expr_mat_filtered <- expr_mat_filtered[available_genes, ]
    
    # Select only the common genes that exist in at least 2 species
    # Select only the columns of interest
    ath_sdf <- select(ath_selected, Orthogroup, ath_gene)
    maize_sdf <- select(maize_selected, Orthogroup, maize_gene)
    rice_sdf <- select(rice_selected, Orthogroup, rice_gene)
    
    # Inner join the data frames based on the Orthogroup column
    combined_df <- inner_join(ath_sdf, maize_sdf, by = "Orthogroup", multiple = "all") %>% 
      inner_join(rice_sdf, by = "Orthogroup", multiple = "all")
    
    # Select only the common genes that exist in at least 2 species
    common_df <- combined_df %>%
      group_by(Orthogroup) %>%
      filter(n() >= 2)
    
    if (length(expr_mat_filtered) > 0) {
      # Filter available_genes with genes that also exist in common_df
      available_genes <- available_genes[available_genes %in% unique(common_df[[switch(input$species, "ATH" = "ath_gene", "Maize" = "maize_gene", "Rice" = "rice_gene")]])]
      
      selectInput("gene_dropdown", "Select a gene:", choices = available_genes)
    }

  })
  
  # observe changes in gene_text and set the selected value of gene_dropdown accordingly
  observeEvent(input$gene_text, {
    updateRadioButtons(session, "gene_input_type", selected = "Text Input") 
    })
  
  # observe changes in gene_text and set the selected value of gene_dropdown accordingly
  observeEvent(input$gene_dropdown, {
    updateRadioButtons(session, "gene_input_type", selected = "Dropdown") 
  })
  
  # Create reactive expression for selected cluster
  observeEvent(input$algorithm, {
    
    if (input$algorithm == "SHAP"){
      print("SHAP")
      ath_selected <- readRDS("omgrds/ath_shap.rds")
      maize_selected <- readRDS("omgrds/maize_shap.rds")
      rice_selected <- readRDS("omgrds/rice_shap.rds")
    }
    else if(input$algorithm == "SVM"){
      print("SVM")
      ath_selected <- readRDS("omgrds/ath_svm.rds")
      maize_selected <- readRDS("omgrds/maize_svm.rds")
      rice_selected <- readRDS("omgrds/rice_svm.rds")
    }else{
      print("Seurat")
      ath_selected <- readRDS("omgrds/ath_seurat.rds")
      maize_selected <- readRDS("omgrds/maize_seurat.rds")
      rice_selected <- readRDS("omgrds/rice_seurat.rds")
    }
    
    #find genes of that orthogroup
    find_corresponding_genes <- reactive ({
      
      if (input$gene_input_type == "Dropdown") {
        if(is.null(input$gene_dropdown)){
          input_gene <- "AT4G35350"
        }
        else{
          input_gene <- input$gene_dropdown
        }
        
      }
      else{
        input_gene <- input$gene_text
      }
      
      # Filter each dataframe to keep only the rows with the input gene
      ath_df <- ath_selected %>% filter(ath_gene == input_gene)
      maize_df <- maize_selected %>% filter(maize_gene == input_gene)
      rice_df <- rice_selected %>% filter(rice_gene == input_gene)
      
      # Get all Orthogroups that contain the input gene
      orthogroups <- unique(c(ath_df$Orthogroup, maize_df$Orthogroup, rice_df$Orthogroup))
      
      if (length(orthogroups) > 0)
      {
        output$matched_omg <- renderText({orthogroups})
      }
      else{
        output$matched_omg <- renderText({"Nothing Matched"}) 
      }
      
      # Initialize empty vectors to store the corresponding genes for each species
      ath_genes <- c()
      maize_genes <- c()
      rice_genes <- c()
      
      # Iterate over each Orthogroup and find the corresponding genes for each species
      for (og in orthogroups) {
        
        # Filter each dataframe to keep only the rows with the current Orthogroup
        ath_og_df <- ath_selected %>% filter(Orthogroup == og)
        maize_og_df <- maize_selected %>% filter(Orthogroup == og)
        rice_og_df <- rice_selected %>% filter(Orthogroup == og)
        
        # If any of the dataframes is empty, skip this Orthogroup
        if (nrow(ath_og_df) == 0 & nrow(maize_og_df) == 0 & nrow(rice_og_df) == 0) {
          next
        }
        
        # Get the first gene for each species from the current Orthogroup
        ath_genes <- c(ath_genes, head(ath_og_df$ath_gene, 1))
        maize_genes <- c(maize_genes, head(maize_og_df$maize_gene, 1))
        rice_genes <- c(rice_genes, head(rice_og_df$rice_gene, 1))
      }
      
      
      # Get the first gene for each species from the current Orthogroup
      if (dim(ath_df)[1] != 0) {
        ath_genes <- c(ath_df$ath_gene)
      }
      if (dim(maize_df)[1] != 0) {
        maize_genes <- c(maize_df$maize_gene)
      }
      
      if (dim(rice_df)[1] != 0) {
        rice_genes <- c(rice_genes, rice_df$rice_gene)
      }
      
      # Return the corresponding genes for each species
      return(list(ath_genes, maize_genes, rice_genes))
    })
    #For ATH
    # Create expression for gene expression data
    ath_gene_expr <- function(){
      ath_gene_of_interest <- find_corresponding_genes()[[1]]
      if (length(ath_gene_of_interest) >= 1){
        ath_gene_of_interest <- ath_gene_of_interest[[1]]
        if ( ath_gene_of_interest %in% rownames(ath_expr_mat)) {
          
          #get the index of the gene
          gene_idx <- which(rownames(ath_expr_mat) == ath_gene_of_interest)
          
          # Extract expression data for the gene
          gene_expr_t  <- ath_expr_mat[gene_idx, ]
          gene_expr_t <- as.data.frame(gene_expr_t)
          gene_expr_t$Cell <- rownames(gene_expr_t)
          colnames(gene_expr_t)[-ncol(gene_expr_t)] <- "expression"
          
          # Merge expression data with cluster info
          merged_data <- merge(ath_cluster_info, gene_expr_t, by = "Cell")
          
        } else {
          data.frame(x = 1)
        }
      }
      else{
        data.frame(x = 1)
      }
      
    }
    
    # Create UMAP plot
    output$umap1 <- renderPlot({
      
      #selected_filters()
      
      ggplot(ath_cluster_info, aes(x = UMAP1, y = UMAP2, color = Cluster)) +
        geom_point(size = 4) +
        ggtitle(paste0("UMAP of Arabidopsis"))+
        theme_bw()
      
    })
    
    output$exp1 <- renderPlot({
      ath_gene_exp_mat <- ath_gene_expr()
      if (length(ath_gene_exp_mat) == 1) {
        ggplot(ath_cluster_info, aes(x = UMAP1, y = UMAP2, fill = "unknown gene"))+
          geom_point(size = 4, shape = 1, color = "grey") +
          ggtitle(paste0("Gene/OMG not present"))+
          theme_bw()
        
        
      }
      else{
        ggplot(ath_gene_exp_mat, aes(x = UMAP1, y = UMAP2, color = expression)) +
          geom_point(size = 4) +
          theme_void() +
          scale_color_gradientn(colors = rev(grey.colors(4))) +
          labs(color = "Expression level") +
          ggtitle(paste0("Feature plot of the input gene"))+
          theme_bw()
        
      }
    })
    
    #For Maize
    # Create expression for gene expression data
    maize_gene_expr <- function(){
      
      maize_gene_of_interest <- find_corresponding_genes()[[2]]
      if (length(maize_gene_of_interest) >= 1){
        maize_gene_of_interest <- maize_gene_of_interest[[1]]
        if ( maize_gene_of_interest %in% rownames(maize_expr_mat)) {
          
          #get the index of the gene
          gene_idx <- which(rownames(maize_expr_mat) == maize_gene_of_interest)
          
          # Extract expression data for the gene
          gene_expr_t  <- maize_expr_mat[gene_idx, ]
          gene_expr_t <- as.data.frame(gene_expr_t)
          gene_expr_t$Cell <- rownames(gene_expr_t)
          colnames(gene_expr_t)[-ncol(gene_expr_t)] <- "expression"
          
          # Merge expression data with cluster info
          merged_data <- merge(maize_cluster_info, gene_expr_t, by = "Cell")
          
        } else {
          data.frame(x = 1)
        }
      }
      else{
        data.frame(x = 1)
      }
      
    }
    
    # Create UMAP plot
    output$umap2 <- renderPlot({
      ggplot(maize_cluster_info, aes(x = UMAP1, y = UMAP2, color = Cluster)) +
        geom_point(size = 4) +
        ggtitle(paste0("UMAP of Maize"))+
        theme_bw()
      
    })
    
    output$exp2 <- renderPlot({
      maize_gene_exp_mat <- maize_gene_expr()
      if (length(maize_gene_exp_mat) == 1) {
        ggplot(maize_cluster_info, aes(x = UMAP1, y = UMAP2, fill = "unknown gene"))+
          geom_point(size = 4, shape = 1, color = "lightgrey") +
          ggtitle(paste0("Gene/OMG not present"))+
          theme_bw()
        
        
      }
      else{
        ggplot(maize_gene_exp_mat, aes(x = UMAP1, y = UMAP2, color = expression)) +
          geom_point(size = 4) +
          theme_void() +
          scale_color_gradientn(colors = rev(grey.colors(4))) +
          labs(color = "Expression level") +
          ggtitle(paste0("Feature plot of the input gene"))+
          theme_bw()
        
      }
    })
    
    
    #For Rice
    # Create expression for gene expression data
    rice_gene_expr <- function(){
      rice_gene_of_interest <- find_corresponding_genes()[[3]]
      if (length(rice_gene_of_interest) >= 1){
        rice_gene_of_interest <- rice_gene_of_interest[[1]]
        if ( rice_gene_of_interest %in% rownames(rice_expr_mat)) {
          
          #get the index of the gene
          gene_idx <- which(rownames(rice_expr_mat) == rice_gene_of_interest)
          
          # Extract expression data for the gene
          gene_expr_t  <- rice_expr_mat[gene_idx, ]
          gene_expr_t <- as.data.frame(gene_expr_t)
          gene_expr_t$Cell <- rownames(gene_expr_t)
          colnames(gene_expr_t)[-ncol(gene_expr_t)] <- "expression"
          
          # Merge expression data with cluster info
          merged_data <- merge(rice_cluster_info, gene_expr_t, by = "Cell")
          
          
        } else {
          data.frame(x = 1)
        }
      }
      else{
        data.frame(x = 1)
      }
      
    }
    
    # Create UMAP plot
    output$umap3 <- renderPlot({
      ggplot(rice_cluster_info, aes(x = UMAP1, y = UMAP2, color = Cluster)) +
        geom_point(size = 4) +
        ggtitle(paste0("UMAP of Rice"))+
        theme_bw()
      
    })
    
    output$exp3 <- renderPlot({
      rice_gene_exp_mat <- rice_gene_expr()
      if (length(rice_gene_exp_mat) == 1) {
        ggplot(rice_cluster_info, aes(x = UMAP1, y = UMAP2, fill = "unknown gene"))+
          geom_point(size = 4, shape = 1, color = "lightgrey") +
          ggtitle(paste0("Gene/OMG not present"))+
          theme_bw()
        
        
      }
      else{
        ggplot(rice_gene_exp_mat, aes(x = UMAP1, y = UMAP2, color = expression)) +
          geom_point(size = 4) +
          theme_void() +
          scale_color_gradientn(colors = rev(grey.colors(4))) +
          labs(color = "Expression level") +
          ggtitle(paste0("Feature plot of the input gene"))+
          theme_bw()
        
      }
    })
    
    
  })

  
}

# Run the app
shinyApp(ui = ui, server = server)
