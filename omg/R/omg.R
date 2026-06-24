# Flexible column-name matching for the query CSV.
normalize_col <- function(df, target, aliases) {
  found <- intersect(tolower(colnames(df)), tolower(aliases))
  if (length(found) == 0) return(NULL)
  original_name <- colnames(df)[tolower(colnames(df)) == found[1]]
  if (original_name != target) {
    colnames(df)[colnames(df) == original_name] <- target
  }
  df
}

#' Annotate cell types by cross-species ortholog marker overlap
#'
#' Single-call equivalent of the \code{Omg_annotation.R} command-line pipeline.
#' Maps query marker genes to orthogroups, tests overlap against a reference
#' panel of species/tissues (Fisher exact test + BH FDR), runs a combined
#' all-species comparison, and consolidates the significant hits into one
#' cell-type prediction per query cluster.
#'
#' @param query Path to a marker CSV (e.g. Seurat \code{FindAllMarkers()}), or a
#'   data frame. Needs gene, cluster and avg_log2FC columns (common aliases are
#'   accepted automatically).
#' @param fdr FDR threshold for significance (default 0.01).
#' @param top_n Top N marker genes per cluster, ranked by avg_log2FC, kept before
#'   orthogroup mapping (default 200).
#' @param reference_markers Path to the reference markers TSV (default: bundled).
#' @param orthogroups Path to the wide OrthoFinder orthogroups TSV (default:
#'   bundled). Pass your own (e.g. an OrthoFinder run including a new query
#'   species) to annotate a species not in the bundled data.
#' @param reference_filter Optional data frame with columns `species` and
#'   `tissue` restricting the combined comparison to those combinations (e.g.
#'   \code{omg_paper_15species()}). Default \code{NULL} uses every species/tissue
#'   in the reference (excluding "Unknow"), so added species are included.
#' @param outdir Output directory. Default \code{NULL} writes a timestamped
#'   \code{output_YYYYMMDD_HHMMSS/} folder in the working directory.
#' @param cell_type_groups Named list mapping reference cell types to broad
#'   groups (default \code{\link{omg_cell_type_groups}}).
#' @param write_files Write CSVs/PDFs to \code{outdir} (default TRUE).
#' @param verbose Print progress (default TRUE).
#'
#' @return Invisibly, a list with `predictions`, `pairwise`, `combined`,
#'   `extract` and `out_dir`.
#' @export
#' @examples
#' \donttest{
#' ex <- system.file("extdata", "example_markers.csv", package = "omg")
#' res <- omg(ex, write_files = FALSE)
#' head(res$predictions)
#' }
omg <- function(query,
                fdr = 0.01,
                top_n = 200L,
                reference_markers = omg_default_markers(),
                orthogroups = omg_default_orthogroups(),
                reference_filter = NULL,
                outdir = NULL,
                cell_type_groups = omg_cell_type_groups(),
                write_files = TRUE,
                verbose = TRUE) {

  say <- function(...) if (verbose) cat(...)
  top_n <- as.integer(top_n)

  say("=== OMG Comparison Pipeline ===\n")
  say("FDR threshold :", fdr, "\n")
  say("Top N genes   :", top_n, "\n\n")

  if (write_files) {
    if (is.null(outdir)) {
      timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
      outdir <- paste0("output_", timestamp)
    }
    dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  }

  # --- 1. Reference data ---
  say("Loading reference data...\n")
  MG_all <- read_markers(reference_markers)
  OG_full <- read_orthogroups(orthogroups)

  say("Building gene-to-orthogroup mapping...\n")
  Ortho <- build_ortho_lookup(OG_full)
  say("Built lookup:", nrow(Ortho), "gene-to-orthogroup mappings\n")

  ref_pairs <- MG_all %>%
    dplyr::filter(tissue != "Unknow") %>%
    dplyr::distinct(species, tissue)
  say("Found", nrow(ref_pairs), "reference species-tissue combinations\n")

  # --- 2. Query data ---
  say("\nLoading user marker genes...\n")
  if (is.character(query)) {
    if (!file.exists(query)) stop("Input file not found: ", query)
    Sample_MG <- read.csv(query, stringsAsFactors = FALSE)
  } else {
    Sample_MG <- as.data.frame(query)
  }

  Sample_MG <- normalize_col(Sample_MG, "cluster",
                             c("cluster", "clustername", "cluster_name", "clusterName",
                               "Cluster", "ClusterName", "Cluster_Name"))
  if (is.null(Sample_MG) || !"cluster" %in% colnames(Sample_MG))
    stop("No cluster column found.")
  Sample_MG <- normalize_col(Sample_MG, "gene",
                             c("gene", "Gene", "gene_name", "geneName", "gene_id", "geneID"))
  if (!"gene" %in% colnames(Sample_MG)) stop("No gene column found.")
  Sample_MG <- normalize_col(Sample_MG, "avg_log2FC",
                             c("avg_log2fc", "avg_log2FC", "avglog2fc", "log2fc",
                               "log2FC", "avg_logFC", "logfc"))
  if (!"avg_log2FC" %in% colnames(Sample_MG)) stop("No fold-change column found.")

  Sample_MG <- Sample_MG %>%
    dplyr::arrange(cluster, dplyr::desc(avg_log2FC)) %>%
    dplyr::group_by(cluster) %>%
    dplyr::slice_head(n = top_n) %>%
    dplyr::ungroup()
  say("Selected top", top_n, "genes per cluster:", nrow(Sample_MG), "total genes across",
      dplyr::n_distinct(Sample_MG$cluster), "clusters\n")

  # Merge with orthogroups: exact match first, normalized fallback
  Sample_OMG <- merge(Sample_MG, Ortho[, c("Orthogroup", "MarkerGene")],
                      by.x = "gene", by.y = "MarkerGene", all.x = TRUE)
  n_exact <- sum(!is.na(Sample_OMG$Orthogroup))
  if (n_exact < nrow(Sample_MG) * 0.5) {
    Sample_MG$gene_norm <- gsub("-", "_", Sample_MG$gene)
    Sample_OMG <- merge(Sample_MG, Ortho[, c("Orthogroup", "MarkerGene_norm")],
                        by.x = "gene_norm", by.y = "MarkerGene_norm", all.x = TRUE)
    n_norm <- sum(!is.na(Sample_OMG$Orthogroup))
    if (n_norm > n_exact) {
      say("Gene-name normalization matched", n_norm, "genes (vs", n_exact, "exact)\n")
    } else {
      Sample_OMG <- merge(Sample_MG, Ortho[, c("Orthogroup", "MarkerGene")],
                          by.x = "gene", by.y = "MarkerGene", all.x = TRUE)
    }
    Sample_OMG$gene_norm <- NULL
  }
  Sample_OMG <- Sample_OMG %>% dplyr::filter(!is.na(Orthogroup))
  say("User query:", dplyr::n_distinct(Sample_OMG$cluster), "clusters,",
      dplyr::n_distinct(Sample_OMG$Orthogroup), "unique orthogroups\n\n")

  # --- 3. Pairwise comparisons ---
  say("=== Running pairwise comparisons ===\n")
  all_pairwise_results <- list()
  all_pairwise_genes   <- list()
  if (write_files) {
    pairwise_dir <- file.path(outdir, "pairwise")
    dir.create(pairwise_dir, showWarnings = FALSE)
  }

  for (idx in seq_len(nrow(ref_pairs))) {
    sp  <- ref_pairs$species[idx]
    tis <- ref_pairs$tissue[idx]
    label <- paste0(sp, "_", gsub("[^A-Za-z0-9]", "_", tis))
    say("  [", idx, "/", nrow(ref_pairs), "] ", sp, " - ", tis, "...\n")

    ref_OMG <- merge(MG_all %>% dplyr::filter(species == sp, tissue == tis),
                     Ortho, by.x = "gene", by.y = "MarkerGene", all.x = TRUE) %>%
      dplyr::filter(!is.na(Orthogroup)) %>%
      dplyr::arrange(clusterName, dplyr::desc(avg_log2FC)) %>%
      dplyr::group_by(clusterName) %>%
      dplyr::slice_head(n = 200) %>%
      dplyr::ungroup()

    if (nrow(ref_OMG) == 0 || dplyr::n_distinct(ref_OMG$clusterName) == 0) next
    if (nrow(Sample_OMG) == 0 || dplyr::n_distinct(Sample_OMG$cluster) == 0) next

    result <- test_significant(ref_OMG, Sample_OMG, fdr,
                               cl_col1 = "clusterName", cl_col2 = "cluster")
    result$species <- sp
    result$tissue  <- tis
    all_pairwise_results[[label]] <- result

    if (write_files) {
      utils::write.csv(result, file.path(pairwise_dir, paste0(label, "_significance.csv")),
                       row.names = FALSE)
    }

    sig_pairs <- result %>% dplyr::filter(test == "Reject")
    if (nrow(sig_pairs) > 0) {
      genes <- extract_common_genes(ref_OMG, Sample_OMG,
                                    cl_col1 = "clusterName", cl_col2 = "cluster")
      if (nrow(genes) > 0) {
        genes$species <- sp
        genes$tissue  <- tis
        genes <- genes %>%
          dplyr::inner_join(sig_pairs %>% dplyr::select(reference_cluster, query_cluster),
                            by = c("ref_cluster" = "reference_cluster",
                                   "query_cluster" = "query_cluster"))
        all_pairwise_genes[[label]] <- genes
        if (write_files) {
          utils::write.csv(genes, file.path(pairwise_dir, paste0(label, "_significant_genes.csv")),
                           row.names = FALSE)
        }
      }
    }

    if (write_files) {
      p <- generate_heatmap(result, label)
      n_ref <- dplyr::n_distinct(result$reference_cluster)
      n_qry <- dplyr::n_distinct(result$query_cluster)
      ggplot2::ggsave(file.path(pairwise_dir, paste0(label, "_heatmap.pdf")), p,
                      width = max(8, n_ref * 0.6 + 3), height = max(6, n_qry * 0.5 + 3),
                      limitsize = FALSE)
    }
  }

  combined_pairwise <- if (length(all_pairwise_results) > 0)
    dplyr::bind_rows(all_pairwise_results) else data.frame()

  # --- 4. Combined all-species comparison (dynamic; replaces hardcoded list) ---
  say("\n=== Running combined multi-species comparison ===\n")
  query_multi <- MG_all %>% dplyr::filter(tissue != "Unknow")
  if (!is.null(reference_filter)) {
    query_multi <- query_multi %>%
      dplyr::semi_join(reference_filter, by = c("species", "tissue"))
  }
  query_multi$clusterName <- paste(query_multi$species, query_multi$tissue, query_multi$clusterName)

  query_multi_OMG <- merge(query_multi, Ortho, by.x = "gene", by.y = "MarkerGene", all.x = TRUE) %>%
    dplyr::filter(!is.na(Orthogroup)) %>%
    dplyr::group_by(species) %>%
    dplyr::arrange(clusterName, dplyr::desc(avg_log2FC)) %>%
    dplyr::group_by(clusterName) %>%
    dplyr::slice_head(n = 200) %>%
    dplyr::ungroup()

  frame_multi <- Heatmap_with_count_comOMG(
    query_multi_OMG,
    Sample_OMG %>% dplyr::mutate(clusterName = as.character(cluster)))
  colnames(frame_multi) <- c("reference_cell_types", "query_clusters",
                             "number_OMGs", "p_value", "negative_log10_p_value")
  frame_multi$significant <- ifelse(frame_multi$p_value < fdr, "Yes", "No")

  if (write_files) {
    utils::write.csv(frame_multi, file.path(outdir, "compare_15species_all.csv"),
                     row.names = FALSE)
    write_combined_heatmap(frame_multi, outdir)
  }

  # --- 5. Cell type predictions ---
  say("\n=== Generating cell type predictions ===\n")
  predictions_summary <- NULL
  extract_all <- NULL
  if (nrow(combined_pairwise) > 0) {
    ref_numOMG_list <- list()
    for (idx in seq_len(nrow(ref_pairs))) {
      sp <- ref_pairs$species[idx]; tis <- ref_pairs$tissue[idx]
      ref_OMG <- merge(MG_all %>% dplyr::filter(species == sp, tissue == tis),
                       Ortho, by.x = "gene", by.y = "MarkerGene", all.x = TRUE) %>%
        dplyr::filter(!is.na(Orthogroup)) %>%
        dplyr::group_by(clusterName) %>% dplyr::slice_head(n = 200) %>% dplyr::ungroup()
      if (nrow(ref_OMG) > 0) {
        ref_numOMG_list[[paste0(sp, "_", tis)]] <- ref_OMG %>%
          dplyr::group_by(clusterName) %>%
          dplyr::summarise(ref_numOMG = dplyr::n_distinct(Orthogroup), .groups = "drop") %>%
          dplyr::mutate(species = sp, tissue = tis)
      }
    }
    ref_numOMG_all <- dplyr::bind_rows(ref_numOMG_list)
    query_numOMG <- Sample_OMG %>%
      dplyr::group_by(cluster) %>%
      dplyr::summarise(query_numOMG = dplyr::n_distinct(Orthogroup), .groups = "drop")

    extract_all <- combined_pairwise %>%
      dplyr::left_join(ref_numOMG_all,
                       by = c("reference_cluster" = "clusterName", "species", "tissue")) %>%
      dplyr::left_join(query_numOMG, by = c("query_cluster" = "cluster")) %>%
      dplyr::select(reference_cluster, species, tissue, ref_numOMG,
                    query_cluster, query_numOMG, num_common_OMG, FDR, test) %>%
      dplyr::arrange(FDR)

    sig_hits <- extract_all %>%
      dplyr::filter(test == "Reject") %>%
      dplyr::mutate(query_cluster = as.character(query_cluster),
                    reference_cluster = as.character(reference_cluster))

    if (nrow(sig_hits) > 0) {
      predictions_summary <- consolidate_cell_types(sig_hits, cell_type_groups)
      if (write_files) {
        utils::write.csv(predictions_summary,
                         file.path(outdir, "cell_type_predictions.csv"), row.names = FALSE)
        utils::write.csv(extract_all %>% dplyr::filter(test == "Reject"),
                         file.path(outdir, "extract_table_significant.csv"), row.names = FALSE)
      }
      if (verbose) print_predictions(predictions_summary, Sample_OMG)
    } else {
      say("No significant matches found - cannot predict cell types.\n")
    }
  }

  say("\n=== Pipeline complete! ===\n")
  if (write_files) say("All outputs saved to:", outdir, "/\n")

  invisible(list(
    predictions = predictions_summary,
    pairwise    = combined_pairwise,
    combined    = frame_multi,
    extract     = extract_all,
    out_dir     = if (write_files) outdir else NULL
  ))
}

#' Reference filter reproducing the published 15-species figure
#'
#' Pass to \code{omg(..., reference_filter = omg_paper_15species())} to restrict
#' the combined comparison to the species/tissues used in Chau et al. (2025).
#'
#' @return A data frame with `species` and `tissue` columns.
#' @export
omg_paper_15species <- function() {
  spec <- list(
    arabidopsis_thaliana = c("Flower", "Inflorescence", "Leaf", "Root", "Seed", "Shoot axis apex", "Stem"),
    oryza_sativa = c("Inflorescence", "Leaf", "Pistil", "Root", "Root;Leaf"),
    zea_mays = c("Ear inflorescence", "Inflorescence", "Leaf", "Root", "Shoot axis apex"),
    brassica_rapa = "Leaf",
    catharanthus_roseus = "Leaf",
    fragaria_vesca = "Leaf",
    glycine_max = c("Flower", "Leaf", "Root", "Seed", "Stem"),
    gossypium_bickii = c("Leaf", "Seed"),
    gossypium_hirsutum = c("Flower", "Leaf"),
    manihot_esculenta = c("Leaf", "Root", "Tuberous root"),
    medicago_truncatula = "Root",
    nicotiana_attenuata = c("Corolla", "Flower"),
    populus_alba_var_pyramidalis = "Stem",
    solanum_lycopersicum = c("Leaf", "Root", "Shoot axis apex"),
    sorghum_bicolor = "Root"
  )
  do.call(rbind, lapply(names(spec), function(s)
    data.frame(species = s, tissue = spec[[s]], stringsAsFactors = FALSE)))
}
