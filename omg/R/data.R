#' Paths to the bundled reference data
#'
#' The package ships a slimmed, gzipped reference marker table (all genes, the
#' five columns the pipeline uses) and the orthogroup table in the original wide
#' OrthoFinder format. \code{read.delim()} reads the gzipped files directly.
#'
#' @return A file path (character) to the bundled gzipped TSV.
#' @export
#' @examples
#' omg_default_markers()
#' omg_default_orthogroups()
omg_default_markers <- function() {
  system.file("extdata", "reference_markers.tsv.gz", package = "omg")
}

#' @rdname omg_default_markers
#' @export
omg_default_orthogroups <- function() {
  system.file("extdata", "orthogroups.tsv.gz", package = "omg")
}

# --- internal readers -------------------------------------------------------

# Read a reference marker table (plain or gzipped TSV).
read_markers <- function(path) {
  if (!file.exists(path)) stop("Reference markers file not found: ", path)
  read.delim(path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
}

# Read the wide OrthoFinder orthogroups table (plain or gzipped TSV).
read_orthogroups <- function(path) {
  if (!file.exists(path)) stop("Orthogroups file not found: ", path)
  read.delim(path, header = TRUE, sep = "\t", stringsAsFactors = FALSE,
             check.names = FALSE)
}

# Build a long gene -> orthogroup lookup from the wide orthogroups table.
# Column-agnostic: every column except `Orthogroup` is treated as a species,
# so an added species column is picked up automatically.
build_ortho_lookup <- function(OG_full) {
  if (!"Orthogroup" %in% colnames(OG_full))
    stop("Orthogroups table must contain an 'Orthogroup' column.")
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
  # normalized lookup (hyphens -> underscores) for flexible matching
  Ortho$MarkerGene_norm <- gsub("-", "_", Ortho$MarkerGene)
  Ortho
}

#' List the species present in the reference data
#'
#' Use this before attempting to add a new species: if your species is already
#' in the orthogroups table you only need to add reference markers (or, when it
#' is your query species, nothing at all).
#'
#' @param reference_markers Path to a reference markers TSV (default: bundled).
#' @param orthogroups Path to a wide OrthoFinder orthogroups TSV (default: bundled).
#' @return Invisibly, a list with `marker_species` and `orthogroup_species`.
#'   Also prints a summary.
#' @export
#' @examples
#' \donttest{ omg_list_species() }
omg_list_species <- function(reference_markers = omg_default_markers(),
                             orthogroups = omg_default_orthogroups()) {
  MG <- read_markers(reference_markers)
  OG <- read_orthogroups(orthogroups)
  marker_sp <- sort(unique(MG$species))
  og_sp <- setdiff(colnames(OG), "Orthogroup")
  cat("Reference marker species (", length(marker_sp), "):\n", sep = "")
  cat("  ", paste(marker_sp, collapse = ", "), "\n\n", sep = "")
  cat("Orthogroup species columns (", length(og_sp), "):\n", sep = "")
  cat("  ", paste(og_sp, collapse = ", "), "\n", sep = "")
  invisible(list(marker_species = marker_sp, orthogroup_species = og_sp))
}

#' Check how well a marker set maps to orthogroups
#'
#' Reports the fraction of genes in a marker CSV (query) or reference table that
#' map to an orthogroup. A low rate usually means the gene IDs do not match the
#' IDs used for that species in the orthogroups table.
#'
#' @param markers Path to a CSV/TSV with a `gene` column, or a data frame.
#' @param orthogroups Path to a wide OrthoFinder orthogroups TSV (default: bundled).
#' @param gene_col Name of the gene column (default `"gene"`).
#' @return Invisibly, a list with `n`, `mapped`, `rate`. Also prints a summary.
#' @export
#' @examples
#' \donttest{
#' omg_check_reference(system.file("extdata", "example_markers.csv", package = "omg"))
#' }
omg_check_reference <- function(markers, orthogroups = omg_default_orthogroups(),
                                gene_col = "gene") {
  if (is.character(markers)) {
    df <- if (grepl("\\.tsv(\\.gz)?$", markers)) {
      read.delim(markers, sep = "\t", stringsAsFactors = FALSE)
    } else {
      read.csv(markers, stringsAsFactors = FALSE)
    }
  } else {
    df <- as.data.frame(markers)
  }
  if (!gene_col %in% colnames(df))
    stop("No '", gene_col, "' column found in markers.")
  Ortho <- build_ortho_lookup(read_orthogroups(orthogroups))
  genes <- unique(df[[gene_col]])
  mapped <- sum(genes %in% Ortho$MarkerGene |
                  gsub("-", "_", genes) %in% Ortho$MarkerGene_norm)
  rate <- if (length(genes)) mapped / length(genes) else NA_real_
  cat(sprintf("Genes: %d unique | mapped to orthogroups: %d (%.1f%%)\n",
              length(genes), mapped, 100 * rate))
  if (!is.na(rate) && rate < 0.3)
    cat("  Low mapping rate - check that gene IDs match the orthogroups table.\n")
  invisible(list(n = length(genes), mapped = mapped, rate = rate))
}
