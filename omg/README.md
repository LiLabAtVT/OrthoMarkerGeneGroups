# omg

**Ortholog Marker Gene cross-species cell type annotation** — R package version of the
`Omg_annotation.R` pipeline ([Chau et al., *Nat Commun* 2025](https://doi.org/10.1038/s41467-024-55755-0)).

Annotate cell types in plant scRNA-seq/snRNA-seq data by testing marker-gene
ortholog overlap against a reference panel of plant species and tissues. Takes a
Seurat `FindAllMarkers()` CSV, returns cell-type predictions with confidence scores.

## Install

```r
install.packages(c("dplyr", "tidyr", "tibble", "ggplot2", "reshape2", "magrittr"))
# install.packages("devtools")
devtools::install_github("LiLabAtVT/OrthoMarkerGeneGroups", subdir = "omg")
```

## Quick start

```r
library(omg)

# Default: uses the bundled reference (all genes) and orthogroups
omg("my_cluster_markers.csv")

# Tune FDR and number of genes (same as the Rscript arguments)
omg("my_cluster_markers.csv", fdr = 0.05, top_n = 500)
```

This writes a timestamped `output_YYYYMMDD_HHMMSS/` folder containing
`cell_type_predictions.csv`, `compare_15species_all.csv`, heatmaps, and a
`pairwise/` folder — the same outputs as the command-line pipeline. The call
also returns the results as R objects:

```r
res <- omg("my_cluster_markers.csv")
res$predictions   # one row per query cluster
```

## Input

A CSV with cluster markers. Column names are matched flexibly:

| Needed | Accepted aliases |
|--------|------------------|
| `gene` | `Gene`, `gene_id`, `geneID`, ... |
| `cluster` | `clusterName`, `cluster_name`, ... |
| `avg_log2FC` | `log2FC`, `avg_logFC`, ... |

## Adding a new species

The bundled orthogroups table covers 34 species (more than the 17 with reference
markers). Check first:

```r
omg_list_species()        # which species are already covered
omg_check_reference("my_markers.csv")   # how many of your genes map
```

**If your query species is new** (not in the 34): run OrthoFinder so your
species' genes are placed into orthogroups, then point `omg()` at the result —
no reinstall, the reference panel is left unchanged:

```r
omg("my_markers.csv", orthogroups = "path/to/Orthogroups.tsv")
```

**To add a species into the reference panel** (advanced): also pass an updated
markers table:

```r
omg("my_markers.csv",
    orthogroups       = "path/to/Orthogroups.tsv",
    reference_markers = "path/to/reference_markers.tsv")
```

`reference_markers.tsv` needs 5 columns: `gene`, `avg_log2FC`, `clusterName`,
`species`, `tissue`. Plain or gzipped (`.tsv.gz`) both work.

## Reproducing the published figure

By default the combined comparison uses every species/tissue in the reference.
To restrict it to the exact set from the paper:

```r
omg("my_markers.csv", reference_filter = omg_paper_15species())
```

## Citation

> Chau, T.N., Timilsena, P.R., Bathala, S.P., Kundu, S., Bargmann, B.O.R. & Li, S.
> Orthologous marker groups reveal broad cell identity conservation across plant
> single-cell transcriptomes. *Nature Communications* 16, 201 (2025).
> https://doi.org/10.1038/s41467-024-55755-0
