#!/usr/bin/env Rscript
# =============================================================================
# OMG (Ortholog Marker Gene) Cross-Species Cell Type Annotation Pipeline
# Usage: Rscript Omg_annotation.R <marker_gene.csv> [FDR_threshold] [top_n_genes]
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(reshape2)
  library(pheatmap)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript Omg_annotation.R <marker_gene.csv> [FDR_threshold] [top_n_genes]\n",
       "  marker_gene.csv  - CSV with columns: gene, avg_log2FC, cluster\n",
       "  FDR_threshold     - optional, default 0.01\n",
       "  top_n_genes       - optional, top N marker genes per cluster before orthogroup mapping, default 200")
}

input_csv    <- args[1]
FDR_threshold <- ifelse(length(args) >= 2, as.numeric(args[2]), 0.01)
top_n_genes  <- ifelse(length(args) >= 3, as.integer(args[3]), 200L)

if (!file.exists(input_csv)) stop("Input file not found: ", input_csv)

cat("=== OMG Comparison Pipeline ===\n")
cat("Input file    :", input_csv, "\n")
cat("FDR threshold :", FDR_threshold, "\n")
cat("Top N genes   :", top_n_genes, "\n\n")

# --- Create output directory with date_time ---
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
out_dir <- paste0("output_", timestamp)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# 1. Load reference data
# =============================================================================
setwd("../Omg_annotation")
cat("Loading reference data...\n")
MG_all <- read.delim("reference_markers.tsv", header = TRUE, sep = "\t")
OG_full <- read.delim("orthogroups.tsv", header = TRUE, sep = "\t")

# Build gene-to-orthogroup lookup (long format) from the full orthogroups file
cat("Building gene-to-orthogroup mapping from Orthogroups file...\n")
og_species_cols <- setdiff(colnames(OG_full), "Orthogroup")
Ortho <- do.call(rbind, lapply(og_species_cols, function(sp_col) {
  rows <- OG_full[OG_full[[sp_col]] != "", c("Orthogroup", sp_col)]
  if (nrow(rows) == 0) return(NULL)
  genes <- strsplit(as.character(rows[[sp_col]]), ",\\s*")
  data.frame(
    Orthogroup = rep(rows$Orthogroup, lengths(genes)),
    MarkerGene = trimws(unlist(genes)),
    stringsAsFactors = FALSE
  )
}))
Ortho <- Ortho[!duplicated(Ortho$MarkerGene), ]
# Also create a normalized lookup (hyphens -> underscores) for flexible matching
Ortho$MarkerGene_norm <- gsub("-", "_", Ortho$MarkerGene)
cat("Built lookup:", nrow(Ortho), "gene-to-orthogroup mappings\n")

# The reference species/tissue combinations used in comparisons (exclude "Unknow")
ref_pairs <- MG_all %>%
  filter(tissue != "Unknow") %>%
  distinct(species, tissue)

cat("Found", nrow(ref_pairs), "reference species-tissue combinations:\n")
for (i in 1:nrow(ref_pairs)) {
  cat("  ", ref_pairs$species[i], " - ", ref_pairs$tissue[i], "\n")
}

# =============================================================================
# 2. Load & prepare user query data
# =============================================================================
cat("\nLoading user marker genes...\n")
Sample_MG <- read.csv(input_csv)

# Flexible column name matching
normalize_col <- function(df, target, aliases) {
  found <- intersect(tolower(colnames(df)), tolower(aliases))
  if (length(found) == 0) return(NULL)
  # Use the first match
  original_name <- colnames(df)[tolower(colnames(df)) == found[1]]
  if (original_name != target) {
    colnames(df)[colnames(df) == original_name] <- target
  }
  return(df)
}

# Normalize cluster column: accept cluster, clusterName, cluster_name, clustername, Cluster, etc.
cluster_aliases <- c("cluster", "clustername", "cluster_name", "clusterName", "Cluster", "ClusterName", "Cluster_Name")
Sample_MG <- normalize_col(Sample_MG, "cluster", cluster_aliases)
if (is.null(Sample_MG) || !"cluster" %in% colnames(Sample_MG))
  stop("No cluster column found. Accepted names: ", paste(cluster_aliases, collapse = ", "))

# Normalize gene column
gene_aliases <- c("gene", "Gene", "gene_name", "geneName", "gene_id", "geneID")
Sample_MG <- normalize_col(Sample_MG, "gene", gene_aliases)
if (!"gene" %in% colnames(Sample_MG))
  stop("No gene column found. Accepted names: ", paste(gene_aliases, collapse = ", "))

# Normalize avg_log2FC column
fc_aliases <- c("avg_log2fc", "avg_log2FC", "avglog2fc", "log2fc", "log2FC", "avg_logFC", "logfc")
Sample_MG <- normalize_col(Sample_MG, "avg_log2FC", fc_aliases)
if (!"avg_log2FC" %in% colnames(Sample_MG))
  stop("No fold-change column found. Accepted names: ", paste(fc_aliases, collapse = ", "))

cat("Detected columns: gene=", colnames(Sample_MG)[colnames(Sample_MG) == "gene"],
    ", cluster=", colnames(Sample_MG)[colnames(Sample_MG) == "cluster"],
    ", avg_log2FC=", colnames(Sample_MG)[colnames(Sample_MG) == "avg_log2FC"], "\n")

# Select top N marker genes per cluster BEFORE orthogroup mapping
Sample_MG <- Sample_MG %>%
  arrange(cluster, desc(avg_log2FC)) %>%
  group_by(cluster) %>%
  slice_head(n = top_n_genes) %>%
  ungroup()

cat("Selected top", top_n_genes, "genes per cluster:",
    nrow(Sample_MG), "total genes across",
    n_distinct(Sample_MG$cluster), "clusters\n")

# Merge with orthogroups — try exact match first
Sample_OMG <- merge(Sample_MG, Ortho[, c("Orthogroup", "MarkerGene")],
                    by.x = "gene", by.y = "MarkerGene", all.x = TRUE)
n_exact <- sum(!is.na(Sample_OMG$Orthogroup))

# If few exact matches, retry with normalized gene names (hyphens <-> underscores)
if (n_exact < nrow(Sample_MG) * 0.5) {
  Sample_MG$gene_norm <- gsub("-", "_", Sample_MG$gene)
  Sample_OMG <- merge(Sample_MG, Ortho[, c("Orthogroup", "MarkerGene_norm")],
                      by.x = "gene_norm", by.y = "MarkerGene_norm", all.x = TRUE)
  n_norm <- sum(!is.na(Sample_OMG$Orthogroup))
  if (n_norm > n_exact) {
    cat("Gene name normalization (hyphen/underscore): matched", n_norm,
        "genes (vs", n_exact, "exact matches)\n")
  } else {
    # Exact was better or equal, redo exact merge
    Sample_OMG <- merge(Sample_MG, Ortho[, c("Orthogroup", "MarkerGene")],
                        by.x = "gene", by.y = "MarkerGene", all.x = TRUE)
  }
  Sample_OMG$gene_norm <- NULL
}

Sample_OMG <- Sample_OMG %>%
  filter(!is.na(Orthogroup))

cat("User query:", n_distinct(Sample_OMG$cluster), "clusters,",
    n_distinct(Sample_OMG$Orthogroup), "unique orthogroups\n\n")

# =============================================================================
# Helper functions
# =============================================================================

# Count common orthogroups between two species (returns matrix)
count_com_OG <- function(Species1, Species2, cl_col1 = "cluster", cl_col2 = "cluster") {
  clusters_S1 <- sort(unique(Species1[[cl_col1]]))
  clusters_S2 <- sort(unique(Species2[[cl_col2]]))
  mat <- matrix(0L, nrow = length(clusters_S1), ncol = length(clusters_S2))
  rownames(mat) <- clusters_S1
  colnames(mat) <- clusters_S2
  for (i in seq_along(clusters_S1)) {
    og1 <- Species1$Orthogroup[Species1[[cl_col1]] == clusters_S1[i]]
    for (j in seq_along(clusters_S2)) {
      og2 <- Species2$Orthogroup[Species2[[cl_col2]] == clusters_S2[j]]
      mat[i, j] <- length(intersect(og1, og2))
    }
  }
  return(mat)
}

# Fisher exact test with BH correction (returns melted data.frame)
test_significant <- function(Species1, Species2, fdr_thr,
                             cl_col1 = "cluster", cl_col2 = "cluster") {
  mat <- count_com_OG(Species1, Species2, cl_col1, cl_col2)

  row_sums <- rowSums(mat)
  col_sums <- colSums(mat)
  grand_total <- sum(mat)

  # Fisher exact test per cell
  p_value_mat <- matrix(NA_real_, nrow = nrow(mat), ncol = ncol(mat),
                         dimnames = dimnames(mat))
  for (i in rownames(mat)) {
    for (j in colnames(mat)) {
      a <- mat[i, j]
      b <- col_sums[j] - a
      c <- row_sums[i] - a
      d <- grand_total - a - b - c
      contingency <- matrix(c(a, b, c, d), nrow = 2)
      p_value_mat[i, j] <- fisher.test(contingency, alternative = "greater")$p.value
    }
  }

  # BH correction
  adj_pval <- p_value_mat
  adj_pval[] <- p.adjust(as.vector(p_value_mat), method = "BH")
  conclusion <- ifelse(adj_pval < fdr_thr, "Reject", "Fail")

  # Melt to long format
  melted <- melt(mat)
  melted$FDR <- melt(adj_pval)$value
  melted$test <- melt(conclusion)$value
  colnames(melted) <- c("reference_cluster", "query_cluster", "num_common_OMG", "FDR", "test")
  return(melted)
}

# Extract common gene names for significant pairs
extract_common_genes <- function(Species1, Species2, cl_col1 = "cluster", cl_col2 = "cluster") {
  clusters_S1 <- sort(unique(Species1[[cl_col1]]))
  clusters_S2 <- sort(unique(Species2[[cl_col2]]))
  results <- list()
  for (i in seq_along(clusters_S1)) {
    og1 <- Species1[Species1[[cl_col1]] == clusters_S1[i], ]
    for (j in seq_along(clusters_S2)) {
      og2 <- Species2[Species2[[cl_col2]] == clusters_S2[j], ]
      common_og <- intersect(og1$Orthogroup, og2$Orthogroup)
      if (length(common_og) > 0) {
        ref_genes <- og1 %>% filter(Orthogroup %in% common_og) %>%
          select(gene, Orthogroup, any_of(c("avg_log2FC", cl_col1))) %>%
          mutate(source = "reference")
        query_genes <- og2 %>% filter(Orthogroup %in% common_og) %>%
          select(gene, Orthogroup, any_of(c("avg_log2FC", cl_col2))) %>%
          mutate(source = "query")
        colnames(ref_genes)[colnames(ref_genes) == cl_col1] <- "cluster_name"
        colnames(query_genes)[colnames(query_genes) == cl_col2] <- "cluster_name"
        ref_genes$cluster_name <- as.character(ref_genes$cluster_name)
        query_genes$cluster_name <- as.character(query_genes$cluster_name)
        combined <- bind_rows(ref_genes, query_genes) %>%
          mutate(ref_cluster = clusters_S1[i], query_cluster = clusters_S2[j])
        results[[length(results) + 1]] <- combined
      }
    }
  }
  if (length(results) > 0) bind_rows(results) else data.frame()
}

# Generate heatmap plot
generate_heatmap <- function(test_result, title_x, title_y) {
  # Sort query clusters numerically
  qc_levels <- unique(as.character(test_result$query_cluster))
  qc_levels <- qc_levels[order(as.numeric(qc_levels))]
  test_result$query_cluster <- factor(test_result$query_cluster, levels = qc_levels)

  # Reorder reference clusters so significant pairs form a diagonal
  ref_order <- c()
  for (qc in qc_levels) {
    qc_sig <- test_result %>%
      filter(query_cluster == qc, test == "Reject") %>%
      arrange(FDR)
    new_refs <- setdiff(as.character(qc_sig$reference_cluster), ref_order)
    ref_order <- c(ref_order, new_refs)
  }
  all_refs <- sort(unique(as.character(test_result$reference_cluster)))
  ref_order <- c(ref_order, setdiff(all_refs, ref_order))
  test_result$reference_cluster <- factor(test_result$reference_cluster, levels = ref_order)

  rejected <- subset(test_result, test == "Reject")
  if (nrow(rejected) > 0) {
    rejected$xmin <- as.numeric(rejected$reference_cluster) - 0.5
    rejected$xmax <- as.numeric(rejected$reference_cluster) + 0.5
    rejected$ymin <- as.numeric(rejected$query_cluster) - 0.5
    rejected$ymax <- as.numeric(rejected$query_cluster) + 0.5
  }

  p <- ggplot(test_result, aes(x = reference_cluster, y = query_cluster, fill = num_common_OMG)) +
    geom_tile() +
    geom_text(aes(label = num_common_OMG), color = "black", size = 3) +
    scale_fill_gradient(low = "white", high = "darkgreen", name = "# OMGs") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(color = "black", face = "bold", size = 10,
                                 angle = 90, vjust = 0.4, hjust = 1),
      axis.text.y = element_text(color = "black", face = "bold", size = 10),
      axis.title.x = element_text(size = 14, face = "bold", colour = "darkgreen"),
      axis.title.y = element_text(size = 14, face = "bold", colour = "darkgreen"),
      plot.title = element_text(size = 12, face = "bold", hjust = 0.5)
    ) +
    labs(x = title_x, y = "Query clusters", title = paste(title_x, "vs Query"))

  if (nrow(rejected) > 0) {
    p <- p + geom_rect(data = rejected,
                       aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
                       fill = NA, color = "red", linewidth = 0.8, inherit.aes = FALSE)
  }
  return(p)
}

# =============================================================================
# 3. Pairwise comparisons: user query vs each species-tissue
# =============================================================================
cat("=== Running pairwise comparisons ===\n")

all_pairwise_results <- list()
all_pairwise_genes   <- list()

pairwise_dir <- file.path(out_dir, "pairwise")
dir.create(pairwise_dir, showWarnings = FALSE)

for (idx in 1:nrow(ref_pairs)) {
  sp  <- ref_pairs$species[idx]
  tis <- ref_pairs$tissue[idx]
  label <- paste0(sp, "_", gsub("[^A-Za-z0-9]", "_", tis))
  cat("  [", idx, "/", nrow(ref_pairs), "] ", sp, " - ", tis, "...\n")

  # Subset reference for this species-tissue
  ref_sub <- MG_all %>%
    filter(species == sp, tissue == tis)

  # Merge with orthogroups, top 200 per cluster
  ref_OMG <- merge(ref_sub, Ortho, by.x = "gene", by.y = "MarkerGene", all.x = TRUE) %>%
    filter(!is.na(Orthogroup)) %>%
    arrange(clusterName, desc(avg_log2FC)) %>%
    group_by(clusterName) %>%
    slice_head(n = 200) %>%
    ungroup()

  if (nrow(ref_OMG) == 0 || n_distinct(ref_OMG$clusterName) == 0) {
    cat("    Skipped (no orthogroup matches in reference)\n")
    next
  }
  if (nrow(Sample_OMG) == 0 || n_distinct(Sample_OMG$cluster) == 0) {
    cat("    Skipped (no orthogroup matches in query)\n")
    next
  }

  # Run significance test
  result <- test_significant(ref_OMG, Sample_OMG, FDR_threshold,
                             cl_col1 = "clusterName", cl_col2 = "cluster")
  result$species <- sp
  result$tissue  <- tis
  all_pairwise_results[[label]] <- result

  # Save significance table CSV
  write.csv(result, file.path(pairwise_dir, paste0(label, "_significance.csv")),
            row.names = FALSE)

  # Extract common gene names for significant pairs
  sig_pairs <- result %>% filter(test == "Reject")
  if (nrow(sig_pairs) > 0) {
    genes <- extract_common_genes(ref_OMG, Sample_OMG,
                                  cl_col1 = "clusterName", cl_col2 = "cluster")
    if (nrow(genes) > 0) {
      genes$species <- sp
      genes$tissue  <- tis
      # Keep only genes from significant pairs
      genes <- genes %>%
        inner_join(sig_pairs %>% select(reference_cluster, query_cluster),
                   by = c("ref_cluster" = "reference_cluster",
                          "query_cluster" = "query_cluster"))
      all_pairwise_genes[[label]] <- genes
      write.csv(genes, file.path(pairwise_dir, paste0(label, "_significant_genes.csv")),
                row.names = FALSE)
    }
  }

  # Generate heatmap
  p <- generate_heatmap(result, label)
  n_ref <- n_distinct(result$reference_cluster)
  n_qry <- n_distinct(result$query_cluster)
  plot_w <- max(8, n_ref * 0.6 + 3)
  plot_h <- max(6, n_qry * 0.5 + 3)
  ggsave(file.path(pairwise_dir, paste0(label, "_heatmap.pdf")),
         p, width = plot_w, height = plot_h, limitsize = FALSE)
}

# =============================================================================
# 4. Combined pairwise summary
# =============================================================================
if (length(all_pairwise_results) > 0) {
  combined_pairwise <- bind_rows(all_pairwise_results)
  cat("\nSaved combined pairwise results:",
      nrow(combined_pairwise), "comparisons\n")

  sig_summary <- combined_pairwise %>%
    filter(test == "Reject") %>%
    arrange(FDR)
  cat("Significant pairs (FDR <", FDR_threshold, "):",
      nrow(sig_summary), "\n")
}

if (length(all_pairwise_genes) > 0) {
  combined_genes <- bind_rows(all_pairwise_genes)
  cat("Total significant gene entries:", nrow(combined_genes), "\n")
}

# =============================================================================
# 5. All 15 species comparison at once (same as compare15species.R logic)
# =============================================================================
cat("\n=== Running combined 15-species comparison ===\n")

# Use the same reference filtering as the original compare15species.R
query_15 <- MG_all %>%
  filter(
    (species == "arabidopsis_thaliana" & tissue %in% c("Flower", "Inflorescence", "Leaf", "Root", "Seed", "Shoot axis apex", "Stem")) |
    (species == "oryza_sativa" & tissue %in% c("Inflorescence", "Leaf", "Pistil", "Root", "Root;Leaf")) |
    (species == "zea_mays" & tissue %in% c("Ear inflorescence", "Inflorescence", "Leaf", "Root", "Shoot axis apex")) |
    (species == "brassica_rapa" & tissue == "Leaf") |
    (species == "catharanthus_roseus" & tissue == "Leaf") |
    (species == "fragaria_vesca" & tissue == "Leaf") |
    (species == "glycine_max" & tissue %in% c("Flower", "Leaf", "Root", "Seed", "Stem")) |
    (species == "gossypium_bickii" & tissue %in% c("Leaf", "Seed")) |
    (species == "gossypium_hirsutum" & tissue %in% c("Flower", "Leaf")) |
    (species == "manihot_esculenta" & tissue %in% c("Leaf", "Root", "Tuberous root")) |
    (species == "medicago_truncatula" & tissue == "Root") |
    (species == "nicotiana_attenuata" & tissue %in% c("Corolla", "Flower")) |
    (species == "populus_alba_var_pyramidalis" & tissue == "Stem") |
    (species == "solanum_lycopersicum" & tissue %in% c("Leaf", "Root", "Shoot axis apex")) |
    (species == "sorghum_bicolor" & tissue == "Root")
  )

query_15$clusterName <- paste(query_15$species, query_15$tissue, query_15$clusterName)

query_15_OMG <- merge(query_15, Ortho, by.x = "gene", by.y = "MarkerGene", all.x = TRUE) %>%
  filter(!is.na(Orthogroup)) %>%
  group_by(species) %>%
  arrange(clusterName, desc(avg_log2FC)) %>%
  group_by(clusterName) %>%
  slice_head(n = 200) %>%
  ungroup()

# Use the Heatmap_with_count_comOMG logic from the original code
# (Fisher test on orthogroup overlap lists, not just counts)
# --- Inlined from compare15species.R ---
count_com_OG_number <- function(Species1, Species2) {
  clusters_S1 <- unique(Species1$clusterName)
  clusters_S2 <- unique(Species2$clusterName)
  two_plants <- matrix(nrow = length(clusters_S1), ncol = length(clusters_S2))
  for (i in 1:length(clusters_S1)) {
    for (j in 1:length(clusters_S2)) {
      list_overlap <- intersect(
        Species1[Species1$clusterName == clusters_S1[i], ]$Orthogroup,
        Species2[Species2$clusterName == clusters_S2[j], ]$Orthogroup
      )
      two_plants[i, j] <- length(list_overlap)
    }
  }
  rownames(two_plants) <- unique(Species1$clusterName)
  colnames(two_plants) <- unique(Species2$clusterName)
  return(two_plants[order(rownames(two_plants)), order(colnames(two_plants))])
}

count_com_OG_15sp <- function(Species1, Species2) {
  clusters_S1 <- unique(Species1$clusterName)
  clusters_S2 <- unique(Species2$clusterName)
  two_plants <- data.frame(matrix(vector("list", length(clusters_S1) * length(clusters_S2)),
                                  nrow = length(clusters_S1),
                                  ncol = length(clusters_S2)))
  for (i in 1:length(clusters_S1)) {
    for (j in 1:length(clusters_S2)) {
      list_overlap <- intersect(
        Species1[Species1$clusterName == clusters_S1[i], ]$Orthogroup,
        Species2[Species2$clusterName == clusters_S2[j], ]$Orthogroup
      )
      two_plants[[i, j]] <- list(list_overlap)
    }
  }
  rownames(two_plants) <- unique(Species1$clusterName)
  colnames(two_plants) <- unique(Species2$clusterName)
  return(two_plants[order(rownames(two_plants)), order(colnames(two_plants))])
}

Heatmap_with_count_comOMG <- function(reference_data, query_data) {
  df <- count_com_OG_15sp(reference_data, query_data)
  colnames(df) <- gsub(" ", "_", colnames(df))

  # Add row that sums for each column
  new_row <- list()
  for (col in names(df)) {
    concatenated_list <- unlist(df[[col]], use.names = FALSE)
    new_row[[col]] <- list(concatenated_list)
  }
  new_row_df <- t(enframe(new_row, name = NULL)) %>% as_tibble()
  colnames(new_row_df) <- colnames(df)
  row.names(new_row_df) <- "new_row"
  df <- rbind(df, new_row_df)

  # Add column that sums for each row
  new_col <- vector("list", nrow(df))
  for (i in 1:nrow(df)) {
    concatenated_list <- list()
    for (col in names(df)) {
      concatenated_list <- c(concatenated_list, unlist(df[[col]][[i]]))
    }
    new_col[[i]] <- unlist(concatenated_list, recursive = TRUE)
  }
  df$new_col <- new_col

  p_value_dataframe <- df[1:(nrow(df) - 1), 1:(ncol(df) - 1)]
  for (i in rownames(p_value_dataframe)) {
    for (j in colnames(p_value_dataframe)) {
      frame <- df[c(i, "new_row"), c(j, "new_col")]
      frame_int <- data.frame(Column1 = c(0, 0), Column2 = c(0, 0))
      frame_int[2, 2] <- length(setdiff(setdiff(setdiff(frame[2, 2][[1]], frame[1, 2][[1]]), frame[2, 1][[1]][[1]]), frame[1, 1][[1]][[1]]))
      frame_int[1, 2] <- length(setdiff(frame[1, 2][[1]], frame[1, 1][[1]][[1]]))
      frame_int[2, 1] <- length(setdiff(frame[2, 1][[1]][[1]], frame[1, 1][[1]][[1]]))
      frame_int[1, 1] <- length(frame[1, 1][[1]][[1]])
      p_value_dataframe[i, j] <- fisher.test(frame_int, alternative = "greater")$p.value
    }
  }
  adjusted_pvalue_dataframe <- p_value_dataframe
  adjusted_pvalue_dataframe[] <- p.adjust(unlist(p_value_dataframe), method = "BH")

  frame <- melt(count_com_OG_number(reference_data, query_data))
  frame$p_value <- melt(adjusted_pvalue_dataframe)$value
  frame$transform <- -log10(frame$p_value)
  return(frame)
}
# --- End inlined compare15species.R ---

frame_15sp <- Heatmap_with_count_comOMG(query_15_OMG, Sample_OMG %>% mutate(clusterName = as.character(cluster)))
colnames(frame_15sp) <- c("reference_cell_types", "query_clusters",
                           "number_OMGs", "p_value", "negative_log10_p_value")
write.csv(frame_15sp, file.path(out_dir, "compare_15species_all.csv"), row.names = FALSE)

# Create heatmap for the 15-species comparison
frame_15sp$significant <- ifelse(frame_15sp$p_value < FDR_threshold, "Yes", "No")

qc15_levels <- unique(as.character(frame_15sp$query_clusters))
qc15_levels <- qc15_levels[order(as.numeric(qc15_levels))]
frame_15sp$query_clusters <- factor(frame_15sp$query_clusters, levels = qc15_levels)

# Reorder reference cell types so significant pairs form a diagonal
# For each query cluster (0, 1, 2, ...), find its best-matching ref cell types
# and place them in that order on the x-axis
ref_order <- c()
for (qc in qc15_levels) {
  qc_data <- frame_15sp %>%
    filter(query_clusters == qc, significant == "Yes") %>%
    arrange(p_value)
  new_refs <- setdiff(as.character(qc_data$reference_cell_types), ref_order)
  ref_order <- c(ref_order, new_refs)
}
# Append any remaining ref cell types not yet placed
all_refs <- sort(unique(as.character(frame_15sp$reference_cell_types)))
ref_order <- c(ref_order, setdiff(all_refs, ref_order))
frame_15sp$reference_cell_types <- factor(frame_15sp$reference_cell_types, levels = ref_order)

p15 <- ggplot(frame_15sp, aes(x = reference_cell_types, y = query_clusters,
                               fill = negative_log10_p_value)) +
  geom_tile() +
  geom_text(aes(label = number_OMGs), color = "black", size = 2) +
  scale_fill_gradient(low = "white", high = "darkgreen", name = "-log10(p)") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(color = "black", size = 6, angle = 90, vjust = 0.4, hjust = 1),
    axis.text.y = element_text(color = "black", size = 10),
    axis.title.x = element_text(size = 12, face = "bold", colour = "darkgreen"),
    axis.title.y = element_text(size = 12, face = "bold", colour = "darkgreen"),
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5)
  ) +
  labs(x = "Reference (15 species)", y = "Query clusters",
       title = "Query vs All 15 Reference Species")

# Mark significant cells
sig_cells <- subset(frame_15sp, significant == "Yes")
if (nrow(sig_cells) > 0) {
  sig_cells$xn <- as.numeric(sig_cells$reference_cell_types)
  sig_cells$yn <- as.numeric(sig_cells$query_clusters)
  p15 <- p15 + geom_rect(data = sig_cells,
                          aes(xmin = xn - 0.5, xmax = xn + 0.5,
                              ymin = yn - 0.5, ymax = yn + 0.5),
                          fill = NA, color = "red", linewidth = 0.6, inherit.aes = FALSE)
}

n_ref15 <- n_distinct(frame_15sp$reference_cell_types)
n_qry15 <- n_distinct(frame_15sp$query_clusters)
ggsave(file.path(out_dir, "compare_15species_heatmap.pdf"),
       p15, width = max(20, n_ref15 * 0.3 + 4), height = max(6, n_qry15 * 0.5 + 3),
       limitsize = FALSE)

# Significant table from 15-species with gene names
cat("Extracting significant genes from 15-species comparison...\n")
sig_15 <- frame_15sp %>% filter(significant == "Yes")

cat("Significant pairs in multi-species comparison:", nrow(sig_15), "\n")

# =============================================================================
# 6. Cell type prediction summary
# =============================================================================
cat("\n=== Generating cell type predictions ===\n")

if (length(all_pairwise_results) > 0) {
  combined_pairwise <- bind_rows(all_pairwise_results)

  # --- Build extract_table style summary across all species ---
  # Number of OMGs per reference cluster (per species-tissue)
  ref_numOMG_list <- list()
  for (idx in 1:nrow(ref_pairs)) {
    sp  <- ref_pairs$species[idx]
    tis <- ref_pairs$tissue[idx]
    ref_sub <- MG_all %>% filter(species == sp, tissue == tis)
    ref_OMG <- merge(ref_sub, Ortho, by.x = "gene", by.y = "MarkerGene", all.x = TRUE) %>%
      filter(!is.na(Orthogroup)) %>%
      group_by(clusterName) %>%
      slice_head(n = 200) %>%
      ungroup()
    if (nrow(ref_OMG) > 0) {
      ref_numOMG_list[[paste0(sp, "_", tis)]] <- ref_OMG %>%
        group_by(clusterName) %>%
        summarise(ref_numOMG = n_distinct(Orthogroup), .groups = "drop") %>%
        mutate(species = sp, tissue = tis)
    }
  }
  ref_numOMG_all <- bind_rows(ref_numOMG_list)

  # Number of OMGs per query cluster
  query_numOMG <- Sample_OMG %>%
    group_by(cluster) %>%
    summarise(query_numOMG = n_distinct(Orthogroup), .groups = "drop")

  # Build the full extract table
  extract_all <- combined_pairwise %>%
    left_join(ref_numOMG_all,
              by = c("reference_cluster" = "clusterName", "species", "tissue")) %>%
    left_join(query_numOMG,
              by = c("query_cluster" = "cluster")) %>%
    select(reference_cluster, species, tissue, ref_numOMG,
           query_cluster, query_numOMG,
           num_common_OMG, FDR, test) %>%
    arrange(FDR)

  extract_sig <- extract_all %>% filter(test == "Reject")

  # --- Print the significant table to console ---
  if (nrow(extract_sig) > 0) {
    cat("\n")
    cat(sprintf("%-30s %8s  %-25s %8s  %6s  %12s  %s\n",
                "Reference_cluster", "Ref_OMG", "Query_cluster", "Qry_OMG",
                "comOMG", "FDR", "Test"))
    cat(paste(rep("-", 120), collapse = ""), "\n")
    n_print <- min(nrow(extract_sig), 50)
    for (r in 1:n_print) {
      cat(sprintf("%-30s %8d  %-25s %8d  %6d  %12s  %s\n",
                  substr(paste0(extract_sig$reference_cluster[r],
                                " (", extract_sig$species[r], ")"), 1, 30),
                  extract_sig$ref_numOMG[r],
                  substr(as.character(extract_sig$query_cluster[r]), 1, 25),
                  extract_sig$query_numOMG[r],
                  extract_sig$num_common_OMG[r],
                  formatC(extract_sig$FDR[r], format = "e", digits = 2),
                  extract_sig$test[r]))
    }
    if (nrow(extract_sig) > 50) {
      cat(sprintf("  ... and %d more rows\n",
                  nrow(extract_sig) - 50))
    }
    cat("\n")
  }

  # --- Cell type predictions: all significant matches per query cluster ---
  sig_hits <- extract_sig %>%
    mutate(query_cluster = as.character(query_cluster),
           reference_cluster = as.character(reference_cluster))

  if (nrow(sig_hits) > 0) {
    # --- Cell type grouping map: all 123 reference cell types -> broad labels ---
    cell_type_groups <- list(
      # --- Epidermis ---
      "Epidermis"      = c("Epidermis", "Leaf epidermis", "Shoot system epidermis",
                           "Abaxial epidermis", "Adaxial epidermis",
                           "Leaf pavement cell", "Pavement cell-A", "Pavement cell-N",
                           "Leaf guard cell", "Guard cell", "Leaf subsidiary cell",
                           "Meristem epidermis", "Outer cell layer",
                           "Cork cambium", "Root epidermis",
                           "Root epidermis/Root cortex"),
      "Root hair cell"  = c("Root hair", "Trichoblast"),
      "Non-hair cell"   = c("Non-hair", "Atrichoblast", "Non-hair root epidermal cell"),
      "Trichome"       = c("Trichome", "Idioblast cell", "Pigment gland"),
      # --- Mesophyll ---
      "Mesophyll"      = c("Mesophyll", "Palisade mesophyll", "Spongy mesophyll",
                           "Chlorenchyma", "Photosynthetic cell",
                           "Mesophyll/Cortex", "Bundle sheath", "Mestome sheath"),
      # --- Phloem ---
      "Phloem"         = c("Phloem", "Companion cell", "Sieve element",
                           "Protophloem", "Phloem parenchyma",
                           "Phloem pole pericycle", "Phloem/Pericycle",
                           "Phloem_Procamb"),
      # --- Xylem ---
      "Xylem"          = c("Xylem", "Metaxylem", "Protoxylem",
                           "Xylem parenchyma", "Xylem precursor cell",
                           "Xylem pole pericycle", "Fiber cell", "Sclerenchyma"),
      # --- Vascular (general) ---
      "Vascular"       = c("Vascular tissue", "Vascular cambium", "Vascular bundle",
                           "Procambium", "Root procambium", "Provascular",
                           "Explant vasculature and callus founder cell",
                           "Hydathodes"),
      # --- Cortex ---
      "Cortex"         = c("Cortex", "Root cortex", "Inner cortex", "Outer cortex",
                           "Cortex/Endodermis", "Cortex/Endodermis initial",
                           "Shoot system cortex", "Ear cortex"),
      # --- Endodermis ---
      "Endodermis"     = c("Endodermis", "Root endodermis",
                           "Shoot system endodermis", "Mature endodermis",
                           "Mesophyll precursor/Root endodermis"),
      # --- Exodermis ---
      "Exodermis"      = c("Exodermis", "Root exodermis",
                           "Root exodermis/Sclerenchym"),
      # --- Pericycle ---
      "Pericycle"      = c("Pericycle", "Pericycle cell",
                           "Xylem pole pericycle cell", "Phloem pole pericycle"),
      # --- Stele ---
      "Stele"          = c("Root stele", "Mature stele"),
      # --- Root cap ---
      "Root cap"       = c("Root cap", "Columella root cap", "Lateral root cap",
                           "Lateral root cap like cell", "Root cap junction",
                           "Columella root cap cell",
                           "RC_Col_QC", "Root broder"),
      # --- Root (other) ---
      "Root"           = c("Root pith", "Lateral root primordia"),
      # --- Meristematic ---
      "Meristematic"   = c("Meristematic cell", "Stem cell niche", "Initial cell",
                           "Initials", "Root initial cell",
                           "Shoot apical meristem", "Inflorescence meristem",
                           "Root meristem", "Root apical meristem",
                           "Shoot meristem initial cell", "Root system ground meristem",
                           "Meristem base", "Meristem boundary",
                           "Leaf primordium", "Leaf rim",
                           "Leaf marginal meristem", "Determinate lateral organ",
                           "Rib Zone", "MZ", "Proximal meristem",
                           "Adaxial meristem periphery",
                           "Branch meristem", "Flower meristem",
                           "Transitory",
                           "Proliferating cell", "Dividing protoderm"),
      # --- Cell cycle (subset of meristematic — more specific labels) ---
      "Cell cycle"     = c("S phase", "S-phase cell", "G2/M phase", "G2/M-phase cell", "G1/G0 phase"),
      # --- Pollen/Spore ---
      "Pollen/Spore"   = c("Pollen", "Sperm nuclei", "Generative nuclei",
                           "Microspore nuclei", "Vegetative nuclei"),
      # --- Parenchyma ---
      "Parenchyma"     = c("Parenchyma", "Pith", "Ray parenchyma",
                           "Inner cell layer"),
      # --- Inflorescence/Spikelet ---
      "Inflorescence"  = c("Inflorescence axis", "Spikelet",
                           "Spikelet floret", "Spikelet meristem"),
      # --- Ovule/Seed ---
      "Ovule/Seed"     = c("Ovule", "CT-Ovule", "S-Ovule", "Nucellus",
                           "Integument", "Inner ovary wall/Outer integument",
                           "Outer ovary wall", "Ovary wall"),
      # --- Flower ---
      "Flower"         = c("Stigma", "Style", "Stamen residue"),
      # --- Nodule/Symbiosis ---
      "Nodule"         = c("Infected cell", "Uninfected cell"),
      # --- Other/Artifact ---
      "Other"          = c("Contaminating nuclei", "Contamination",
                           "Stress response", "Outer pigment layer",
                           "Middle cell layer")
    )

    # Lookup function: map a cell type to its broad group
    get_broad_group <- function(ct) {
      for (group_name in names(cell_type_groups)) {
        if (ct %in% cell_type_groups[[group_name]]) return(group_name)
      }
      return(ct)  # keep original if no group match
    }

    # Build summary table: one row per query cluster
    # Count frequency of each cell type across species, sort high to low
    predictions_summary <- sig_hits %>%
      group_by(query_cluster) %>%
      summarise(
        cell_type_prediction_frequency = {
          # Count frequency of each reference cell type
          ct_freq <- sort(table(reference_cluster), decreasing = TRUE)
          paste(paste0(names(ct_freq), " (", ct_freq, ")"), collapse = "; ")
        },
        consolidated_cell_type = {
          # Map each hit to its broad group, then pick the most frequent group
          broad <- sapply(reference_cluster, get_broad_group)

          # Reassign "Epidermis" based on Root hair vs Non-hair counts
          rh_count <- sum(broad == "Root hair cell")
          nh_count <- sum(broad == "Non-hair cell")
          if (sum(broad == "Epidermis") > 0 && (rh_count > 0 || nh_count > 0)) {
            if (rh_count > nh_count) {
              broad[broad == "Epidermis"] <- "Root hair cell"
            } else if (nh_count > rh_count) {
              broad[broad == "Epidermis"] <- "Non-hair cell"
            }
          }

          # Reassign "Epidermis" if Cortex/Exodermis/Endodermis are more dominant
          epi_count <- sum(broad == "Epidermis")
          ground_counts <- c(Cortex = sum(broad == "Cortex"),
                             Exodermis = sum(broad == "Exodermis"),
                             Endodermis = sum(broad == "Endodermis"))
          total_ground <- sum(ground_counts)
          if (epi_count > 0 && total_ground > epi_count) {
            top_ground <- names(which.max(ground_counts))
            broad[broad == "Epidermis"] <- top_ground
          }

          # --- Major layer grouping to reduce ties ---
          major_layer_map <- c(
            "Root cap"       = "Outer", "Epidermis"      = "Outer",
            "Root hair cell" = "Outer", "Non-hair cell"  = "Outer",
            "Trichome"       = "Outer",
            "Exodermis"      = "Middle", "Cortex"        = "Middle",
            "Endodermis"     = "Middle",
            "Pericycle"      = "Inner", "Stele"          = "Inner",
            "Phloem"         = "Inner", "Xylem"          = "Inner",
            "Vascular"       = "Inner"
          )

          # Map broad groups to major layers
          major <- major_layer_map[broad]

          # Count hits per major layer (ignore NA = groups not in the 3 layers)
          layer_counts <- table(major[!is.na(major)])

          group_freq <- sort(table(broad), decreasing = TRUE)
          top_count <- group_freq[1]
          tied <- names(group_freq[group_freq == top_count])

          # If multiple tied groups, use major layer to reduce
          if (length(tied) > 1 && length(layer_counts) > 0) {
            # Find dominant major layer
            top_layer <- names(sort(layer_counts, decreasing = TRUE))[1]
            # Keep only tied groups that belong to the dominant layer
            tied_in_layer <- tied[tied %in% names(major_layer_map[major_layer_map == top_layer])]
            if (length(tied_in_layer) > 0) tied <- tied_in_layer
          }

          # Override Meristematic with cell cycle phase if G2/M or S phase counts are high
          # Compare against top individual cell type within Meristematic, not the
          # inflated broad group total (which lumps many cell types together)
          if ("Meristematic" %in% tied) {
            g2m_count <- sum(reference_cluster == "G2/M phase" |
                             reference_cluster == "G2/M-phase cell")
            s_count   <- sum(reference_cluster == "S phase" |
                             reference_cluster == "S-phase cell")
            # Get the count of the single most common individual cell type in Meristematic
            merist_types <- cell_type_groups[["Meristematic"]]
            merist_individual_counts <- table(reference_cluster[reference_cluster %in% merist_types])
            top_merist_individual <- if (length(merist_individual_counts) > 0) max(merist_individual_counts) else 0
            # If cell cycle phase hits are at least half of the top individual meristematic type
            if (g2m_count >= top_merist_individual / 2 || s_count >= top_merist_individual / 2) {
              tied <- setdiff(tied, "Meristematic")
              if (g2m_count >= s_count && g2m_count > 0) {
                tied <- c(tied, "G2/M phase")
              } else if (s_count > 0) {
                tied <- c(tied, "S phase")
              }
              if (length(tied) == 0) {
                if (g2m_count >= s_count) tied <- "G2/M phase" else tied <- "S phase"
              }
            }
          }

          # Override Phloem with more specific Sieve element or Companion cell
          if ("Phloem" %in% tied) {
            sieve_count <- sum(grepl("^[Ss]ieve element", reference_cluster))
            companion_count <- sum(grepl("^[Cc]ompanion cell", reference_cluster))
            if (sieve_count > 0 || companion_count > 0) {
              if (sieve_count >= companion_count) {
                tied[tied == "Phloem"] <- "Sieve element"
              } else {
                tied[tied == "Phloem"] <- "Companion cell"
              }
            }
          }

          # Compute prediction confidence: hits for predicted group(s) / total hits
          # S phase/G2M phase are specific labels for meristematic cells,
          # so Meristematic + Cell cycle hits all support the prediction
          count_groups <- tied
          # Meristematic and Cell cycle are related (dividing cells)
          if (any(tied %in% c("S phase", "G2/M phase", "Meristematic"))) {
            count_groups <- unique(c(count_groups, "Meristematic", "Cell cycle"))
          }
          # Sieve element/Companion cell are specific Phloem subtypes
          if (any(tied %in% c("Sieve element", "Companion cell"))) {
            count_groups <- unique(c(count_groups, "Phloem"))
          }
          total_hits <- length(broad)
          # If tied, use the count of one group (they're equal), not all combined
          if (length(tied) > 1) {
            hits_supporting <- as.integer(group_freq[1])
          } else {
            hits_supporting <- sum(broad %in% count_groups)
          }
          conf <- round(hits_supporting / total_hits, 3)

          paste(paste(tied, collapse = "/"), conf, sep = "|||")
        },
        .groups = "drop"
      ) %>%
      # Split consolidated_cell_type into cell type and confidence
      separate(consolidated_cell_type, into = c("consolidated_cell_type", "prediction_confidence"),
               sep = "\\|\\|\\|", convert = TRUE) %>%
      arrange(as.numeric(as.character(query_cluster)))

    write.csv(predictions_summary,
              file.path(out_dir, "cell_type_predictions.csv"),
              row.names = FALSE)

    # Print summary to console
    cat("\n=== Cell Type Predictions (one row per query cluster) ===\n")
    for (r in 1:nrow(predictions_summary)) {
      cat(sprintf("\n--- %s ---\n", predictions_summary$query_cluster[r]))
      cat("  Consolidated cell type:         ", predictions_summary$consolidated_cell_type[r], "\n")
      cat("  Prediction confidence:          ", predictions_summary$prediction_confidence[r], "\n")
      cat("  Cell type prediction frequency: ", predictions_summary$cell_type_prediction_frequency[r], "\n")
    }

    # Also list clusters with no significant match
    all_clusters <- sort(unique(as.character(Sample_OMG$cluster)))
    no_match <- setdiff(all_clusters, sort(unique(predictions_summary$query_cluster)))
    if (length(no_match) > 0) {
      cat("\n--- No significant match ---\n")
      for (cl in no_match) cat("  ", cl, "\n")
    }

    # Print explanation of the consolidation logic
    cat("\n===================================================================================\n")
    cat("HOW TO READ THE CONSOLIDATED CELL TYPE PREDICTIONS\n")
    cat("===================================================================================\n")
    cat("\nCOLUMNS:\n")
    cat("  query_cluster                  - Your input cluster ID\n")
    cat("  cell_type_prediction_frequency - All significant reference cell types matched to\n")
    cat("                                   this cluster with hit counts across species.\n")
    cat("                                   Format: CellType (count); CellType (count); ...\n")
    cat("                                   Higher count = more species/tissues agree.\n")
    cat("  consolidated_cell_type         - Final predicted cell type after grouping and\n")
    cat("                                   tie-breaking (explained below).\n")
    cat("\nHOW CONSOLIDATED CELL TYPE IS DETERMINED:\n")
    cat("\n  Step 1: Broad grouping\n")
    cat("    Each reference cell type is mapped to a broad group. For example:\n")
    cat("      Root hair, Trichoblast                  -> Root hair cell\n")
    cat("      Non-hair, Atrichoblast                  -> Non-hair cell\n")
    cat("      Cortex, Root cortex, Inner/Outer cortex -> Cortex\n")
    cat("      Phloem, Companion cell, Sieve element   -> Phloem\n")
    cat("      Xylem, Metaxylem, Protoxylem            -> Xylem\n")
    cat("      (see cell_type_groups in run_all_comparisons.R for the full mapping)\n")
    cat("\n  Step 2: Epidermis reassignment\n")
    cat("    - If Root hair cell hits > Non-hair cell hits: Epidermis -> Root hair cell\n")
    cat("    - If Non-hair cell hits > Root hair cell hits: Epidermis -> Non-hair cell\n")
    cat("\n  Step 3: Major layer tie-breaking\n")
    cat("    Broad groups are classified into three root tissue layers:\n")
    cat("      Outer layer : Root cap, Epidermis, Root hair cell, Non-hair cell, Trichome\n")
    cat("      Middle layer: Exodermis, Cortex, Endodermis\n")
    cat("      Inner layer : Pericycle, Stele, Phloem, Xylem, Vascular\n")
    cat("    The most frequent broad group is selected. If multiple groups are tied,\n")
    cat("    only those belonging to the dominant major layer are kept.\n")
    cat("\n  Step 4: Cell cycle override\n")
    cat("    If the top group is Meristematic but G2/M phase or S phase hits are\n")
    cat("    at least half of the Meristematic count, the cluster is reclassified\n")
    cat("    as 'G2/M phase' or 'S phase' (whichever has more hits).\n")
    cat("\n  Step 5: Final output\n")
    cat("    The remaining group(s) are reported, separated by '/' if still tied.\n")
    cat("\nHOW TO FIND EVIDENCE:\n")
    cat("  1. Look at cell_type_prediction_frequency for which reference cell types\n")
    cat("     matched your cluster and how many species/tissues support each match.\n")
    cat("  2. Check extract_table_significant.csv for significant pairs with FDR values\n")
    cat("     - lower FDR = stronger evidence.\n")
    cat("  3. Check compare_15species_significant_genes.csv or all_significant_genes.csv\n")
    cat("     to see which orthologous marker genes are shared.\n")
    cat("  4. Check pairwise/*_significance.csv for per-species detail.\n")
    cat("===================================================================================\n")
    cat("\n")
  } else {
    cat("No significant matches found — cannot predict cell types.\n")
  }
}

# =============================================================================
# Done
# =============================================================================
cat("\n=== Pipeline complete! ===\n")
cat("All outputs saved to:", out_dir, "/\n")
cat("\nOutput files:\n")
cat("  ", out_dir, "/cell_type_predictions.csv           - Cell type predictions per cluster\n")
cat("  ", out_dir, "/compare_15species_all.csv          - Multi-species comparison table\n")
cat("  ", out_dir, "/compare_15species_heatmap.pdf      - Multi-species heatmap\n")
cat("  ", out_dir, "/pairwise/                         - Per species-tissue results\n")
cat("  ", out_dir, "/pairwise/*_significance.csv        - Significance tables\n")
cat("  ", out_dir, "/pairwise/*_significant_genes.csv   - Gene names for significant pairs\n")
cat("  ", out_dir, "/pairwise/*_heatmap.pdf             - Heatmap plots\n")
