#' omg: Ortholog Marker Gene cross-species cell type annotation
#'
#' The main entry point is \code{\link{omg}}, which reproduces the
#' \code{Omg_annotation.R} command-line pipeline as a single R call.
#' Helper functions \code{\link{omg_list_species}} and
#' \code{\link{omg_check_reference}} support adding new species.
#'
#' @keywords internal
#' @importFrom magrittr %>%
#' @importFrom stats fisher.test p.adjust
#' @importFrom utils read.csv
"_PACKAGE"

# Quiet R CMD check notes about NSE column names used by dplyr/ggplot2.
utils::globalVariables(c(
  "Orthogroup", "MarkerGene", "MarkerGene_norm", "gene", "avg_log2FC",
  "cluster", "clusterName", "species", "tissue", "test", "FDR",
  "num_common_OMG", "reference_cluster", "query_cluster", "query_clusters",
  "reference_cell_types", "number_OMGs", "negative_log10_p_value",
  "significant", "p_value", "ref_numOMG", "query_numOMG", "consolidated_cell_type",
  "value", "xmin", "xmax", "ymin", "ymax", "xn", "yn", "name"
))
