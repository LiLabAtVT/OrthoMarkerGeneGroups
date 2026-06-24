# Combined multi-species heatmap (ported from Omg_annotation.R).
write_combined_heatmap <- function(frame_multi, outdir) {
  qc_levels <- unique(as.character(frame_multi$query_clusters))
  qc_levels <- qc_levels[order(suppressWarnings(as.numeric(qc_levels)))]
  frame_multi$query_clusters <- factor(frame_multi$query_clusters, levels = qc_levels)

  ref_order <- c()
  for (qc in qc_levels) {
    qc_data <- frame_multi %>%
      dplyr::filter(query_clusters == qc, significant == "Yes") %>%
      dplyr::arrange(p_value)
    new_refs <- setdiff(as.character(qc_data$reference_cell_types), ref_order)
    ref_order <- c(ref_order, new_refs)
  }
  all_refs <- sort(unique(as.character(frame_multi$reference_cell_types)))
  ref_order <- c(ref_order, setdiff(all_refs, ref_order))
  frame_multi$reference_cell_types <- factor(frame_multi$reference_cell_types, levels = ref_order)

  p <- ggplot2::ggplot(frame_multi,
                       ggplot2::aes(x = reference_cell_types, y = query_clusters,
                                    fill = negative_log10_p_value)) +
    ggplot2::geom_tile() +
    ggplot2::geom_text(ggplot2::aes(label = number_OMGs), color = "black", size = 2) +
    ggplot2::scale_fill_gradient(low = "white", high = "darkgreen", name = "-log10(p)") +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(color = "black", size = 6, angle = 90, vjust = 0.4, hjust = 1),
      axis.text.y = ggplot2::element_text(color = "black", size = 10),
      axis.title.x = ggplot2::element_text(size = 12, face = "bold", colour = "darkgreen"),
      axis.title.y = ggplot2::element_text(size = 12, face = "bold", colour = "darkgreen"),
      plot.title = ggplot2::element_text(size = 12, face = "bold", hjust = 0.5)
    ) +
    ggplot2::labs(x = "Reference (all species)", y = "Query clusters",
                  title = "Query vs All Reference Species")

  sig_cells <- subset(frame_multi, significant == "Yes")
  if (nrow(sig_cells) > 0) {
    sig_cells$xn <- as.numeric(sig_cells$reference_cell_types)
    sig_cells$yn <- as.numeric(sig_cells$query_clusters)
    p <- p + ggplot2::geom_rect(data = sig_cells,
                                ggplot2::aes(xmin = xn - 0.5, xmax = xn + 0.5,
                                             ymin = yn - 0.5, ymax = yn + 0.5),
                                fill = NA, color = "red", linewidth = 0.6, inherit.aes = FALSE)
  }

  n_ref <- dplyr::n_distinct(frame_multi$reference_cell_types)
  n_qry <- dplyr::n_distinct(frame_multi$query_clusters)
  ggplot2::ggsave(file.path(outdir, "compare_15species_heatmap.pdf"), p,
                  width = max(20, n_ref * 0.3 + 4), height = max(6, n_qry * 0.5 + 3),
                  limitsize = FALSE)
}

# Print the per-cluster prediction summary to the console.
print_predictions <- function(predictions_summary, Sample_OMG) {
  cat("\n=== Cell Type Predictions (one row per query cluster) ===\n")
  for (r in seq_len(nrow(predictions_summary))) {
    cat(sprintf("\n--- %s ---\n", predictions_summary$query_cluster[r]))
    cat("  Consolidated cell type:         ", predictions_summary$consolidated_cell_type[r], "\n")
    cat("  Prediction confidence:          ", predictions_summary$prediction_confidence[r], "\n")
    cat("  Cell type prediction frequency: ", predictions_summary$cell_type_prediction_frequency[r], "\n")
  }
  all_clusters <- sort(unique(as.character(Sample_OMG$cluster)))
  no_match <- setdiff(all_clusters, sort(unique(predictions_summary$query_cluster)))
  if (length(no_match) > 0) {
    cat("\n--- No significant match ---\n")
    for (cl in no_match) cat("  ", cl, "\n")
  }
}
