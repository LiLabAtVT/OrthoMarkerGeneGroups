# =============================================================================
# Significance testing on orthogroup overlap (ported from Omg_annotation.R)
# =============================================================================

# Count common orthogroups between two species (returns integer matrix).
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
  mat
}

# Fisher exact test per cell with BH correction (returns melted data.frame).
test_significant <- function(Species1, Species2, fdr_thr,
                             cl_col1 = "cluster", cl_col2 = "cluster") {
  mat <- count_com_OG(Species1, Species2, cl_col1, cl_col2)

  row_sums <- rowSums(mat)
  col_sums <- colSums(mat)
  grand_total <- sum(mat)

  p_value_mat <- matrix(NA_real_, nrow = nrow(mat), ncol = ncol(mat),
                        dimnames = dimnames(mat))
  for (i in rownames(mat)) {
    for (j in colnames(mat)) {
      a <- mat[i, j]
      b <- col_sums[j] - a
      c <- row_sums[i] - a
      d <- grand_total - a - b - c
      contingency <- matrix(c(a, b, c, d), nrow = 2)
      p_value_mat[i, j] <- stats::fisher.test(contingency, alternative = "greater")$p.value
    }
  }

  adj_pval <- p_value_mat
  adj_pval[] <- stats::p.adjust(as.vector(p_value_mat), method = "BH")
  conclusion <- ifelse(adj_pval < fdr_thr, "Reject", "Fail")

  melted <- reshape2::melt(mat)
  melted$FDR <- reshape2::melt(adj_pval)$value
  melted$test <- reshape2::melt(conclusion)$value
  colnames(melted) <- c("reference_cluster", "query_cluster", "num_common_OMG", "FDR", "test")
  melted
}

# Extract common gene names for significant pairs.
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
        ref_genes <- og1 %>% dplyr::filter(Orthogroup %in% common_og) %>%
          dplyr::select(gene, Orthogroup, dplyr::any_of(c("avg_log2FC", cl_col1))) %>%
          dplyr::mutate(source = "reference")
        query_genes <- og2 %>% dplyr::filter(Orthogroup %in% common_og) %>%
          dplyr::select(gene, Orthogroup, dplyr::any_of(c("avg_log2FC", cl_col2))) %>%
          dplyr::mutate(source = "query")
        colnames(ref_genes)[colnames(ref_genes) == cl_col1] <- "cluster_name"
        colnames(query_genes)[colnames(query_genes) == cl_col2] <- "cluster_name"
        ref_genes$cluster_name <- as.character(ref_genes$cluster_name)
        query_genes$cluster_name <- as.character(query_genes$cluster_name)
        combined <- dplyr::bind_rows(ref_genes, query_genes) %>%
          dplyr::mutate(ref_cluster = clusters_S1[i], query_cluster = clusters_S2[j])
        results[[length(results) + 1]] <- combined
      }
    }
  }
  if (length(results) > 0) dplyr::bind_rows(results) else data.frame()
}

# --- Combined multi-species comparison (orthogroup-list Fisher test) ---------

count_com_OG_number <- function(Species1, Species2) {
  clusters_S1 <- unique(Species1$clusterName)
  clusters_S2 <- unique(Species2$clusterName)
  two_plants <- matrix(nrow = length(clusters_S1), ncol = length(clusters_S2))
  for (i in seq_along(clusters_S1)) {
    for (j in seq_along(clusters_S2)) {
      list_overlap <- intersect(
        Species1[Species1$clusterName == clusters_S1[i], ]$Orthogroup,
        Species2[Species2$clusterName == clusters_S2[j], ]$Orthogroup
      )
      two_plants[i, j] <- length(list_overlap)
    }
  }
  rownames(two_plants) <- unique(Species1$clusterName)
  colnames(two_plants) <- unique(Species2$clusterName)
  two_plants[order(rownames(two_plants)), order(colnames(two_plants)), drop = FALSE]
}

count_com_OG_15sp <- function(Species1, Species2) {
  clusters_S1 <- unique(Species1$clusterName)
  clusters_S2 <- unique(Species2$clusterName)
  two_plants <- data.frame(matrix(vector("list", length(clusters_S1) * length(clusters_S2)),
                                  nrow = length(clusters_S1),
                                  ncol = length(clusters_S2)))
  for (i in seq_along(clusters_S1)) {
    for (j in seq_along(clusters_S2)) {
      list_overlap <- intersect(
        Species1[Species1$clusterName == clusters_S1[i], ]$Orthogroup,
        Species2[Species2$clusterName == clusters_S2[j], ]$Orthogroup
      )
      two_plants[[i, j]] <- list(list_overlap)
    }
  }
  rownames(two_plants) <- unique(Species1$clusterName)
  colnames(two_plants) <- unique(Species2$clusterName)
  two_plants[order(rownames(two_plants)), order(colnames(two_plants))]
}

Heatmap_with_count_comOMG <- function(reference_data, query_data) {
  df <- count_com_OG_15sp(reference_data, query_data)
  colnames(df) <- gsub(" ", "_", colnames(df))

  new_row <- list()
  for (col in names(df)) {
    concatenated_list <- unlist(df[[col]], use.names = FALSE)
    new_row[[col]] <- list(concatenated_list)
  }
  new_row_df <- t(tibble::enframe(new_row, name = NULL)) %>%
    tibble::as_tibble(.name_repair = "minimal")
  colnames(new_row_df) <- colnames(df)
  row.names(new_row_df) <- "new_row"
  df <- rbind(df, new_row_df)

  new_col <- vector("list", nrow(df))
  for (i in seq_len(nrow(df))) {
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
      p_value_dataframe[i, j] <- stats::fisher.test(frame_int, alternative = "greater")$p.value
    }
  }
  adjusted_pvalue_dataframe <- p_value_dataframe
  adjusted_pvalue_dataframe[] <- stats::p.adjust(unlist(p_value_dataframe), method = "BH")

  frame <- reshape2::melt(count_com_OG_number(reference_data, query_data))
  frame$p_value <- reshape2::melt(adjusted_pvalue_dataframe)$value
  frame$transform <- -log10(frame$p_value)
  frame
}
