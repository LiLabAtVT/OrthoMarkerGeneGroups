# Ortholog Marker Gene cross-species cell type annotation

Annotate cell types in plant scRNA-seq/snRNA-seq data by testing marker-gene ortholog overlap against a reference panel of plant species and tissues. 

## Install

```r
install.packages(c("dplyr", "tidyr", "tibble", "ggplot2", "reshape2", "magrittr"))
# install.packages("devtools")
devtools::install_github("LiLabAtVT/OrthoMarkerGeneGroups", subdir = "omg")
```

## Quick start

```r
library(omg)

omg("my_cluster_markers.csv")
omg("my_cluster_markers.csv", fdr = 0.05, top_n = 500)
```

This writes a timestamped `output_YYYYMMDD_HHMMSS/` folder containing
`cell_type_predictions.csv`, `compare_15species_all.csv`, heatmaps, and a
`pairwise/` folder — the same outputs as the command-line pipeline. The call
also returns the results as R objects:


## Input

A CSV with cluster markers. Column names are matched flexibly:

| Needed | Accepted aliases |
|--------|------------------|
| `gene` | `Gene`, `gene_id`, `geneID`, ... |
| `cluster` | `clusterName`, `cluster_name`, ... |
| `avg_log2FC` | `log2FC`, `avg_logFC`, ... |

## Adding a new species

**If your query species is new** (not in the 33): run OrthoFinder so your
species' genes are placed into orthogroups, then point `omg()` at the result —
no reinstall, the reference panel is left unchanged:

```r
omg("my_markers.csv", orthogroups = "path/to/Orthogroups.tsv")
```

**To add a species into the reference panel** (advanced): also pass an updated
markers table:

```r
omg("my_markers.csv", 
    orthogroups = "path/to/Orthogroups.tsv",
    reference_markers = "path/to/reference_markers.tsv")
```

`reference_markers.tsv` needs 5 columns: `gene`, `avg_log2FC`, `clusterName`,
`species`, `tissue`. Plain or gzipped (`.tsv.gz`) both work.


## Citation
> Chau, T.N., Timilsena, P.R., Bathala, S.P., Kundu, S., Bargmann, B.O.R. & Li, S.
> Orthologous marker groups reveal broad cell identity conservation across plant
> single-cell transcriptomes. *Nature Communications* 16, 201 (2025).
> https://doi.org/10.1038/s41467-024-55755-0
