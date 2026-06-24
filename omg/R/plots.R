# =============================================================================
# Heatmap for a pairwise comparison result (ported from Omg_annotation.R)
# =============================================================================

generate_heatmap <- function(test_result, title_x) {
  qc_levels <- unique(as.character(test_result$query_cluster))
  qc_levels <- qc_levels[order(suppressWarnings(as.numeric(qc_levels)))]
  test_result$query_cluster <- factor(test_result$query_cluster, levels = qc_levels)

  # Reorder reference clusters so significant pairs form a diagonal
  ref_order <- c()
  for (qc in qc_levels) {
    qc_sig <- test_result %>%
      dplyr::filter(query_cluster == qc, test == "Reject") %>%
      dplyr::arrange(FDR)
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

  p <- ggplot2::ggplot(test_result,
                       ggplot2::aes(x = reference_cluster, y = query_cluster,
                                    fill = num_common_OMG)) +
    ggplot2::geom_tile() +
    ggplot2::geom_text(ggplot2::aes(label = num_common_OMG), color = "black", size = 3) +
    ggplot2::scale_fill_gradient(low = "white", high = "darkgreen", name = "# OMGs") +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(color = "black", face = "bold", size = 10,
                                          angle = 90, vjust = 0.4, hjust = 1),
      axis.text.y = ggplot2::element_text(color = "black", face = "bold", size = 10),
      axis.title.x = ggplot2::element_text(size = 14, face = "bold", colour = "darkgreen"),
      axis.title.y = ggplot2::element_text(size = 14, face = "bold", colour = "darkgreen"),
      plot.title = ggplot2::element_text(size = 12, face = "bold", hjust = 0.5)
    ) +
    ggplot2::labs(x = title_x, y = "Query clusters", title = paste(title_x, "vs Query"))

  if (nrow(rejected) > 0) {
    p <- p + ggplot2::geom_rect(data = rejected,
                                ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
                                fill = NA, color = "red", linewidth = 0.8, inherit.aes = FALSE)
  }
  p
}
