# OMG: Ortholog Marker Gene Cross-Species Cell Type Annotation

Annotate cell types in plant scRNA-seq/snRNA-seq data by testing marker gene ortholog overlap against a reference panel of 15 species with different tissue types ([Chau et al., *Nat Commun* 2025](https://doi.org/10.1038/s41467-024-55755-0)). Takes a Seurat `FindAllMarkers()` CSV as input, returns cell type predictions with confidence scores.

## Quick Start

```bash
# Clone and set up
git clone https://github.com/LiLabAtVT/OrthoMarkerGeneGroups.git
cd OrthoMarkerGeneGroups/Omg_annotation
unzip "*.zip"

# Install R dependencies (one time)
Rscript -e 'install.packages(c("tidyverse", "reshape2", "pheatmap"))'

# Run the pipeline
Rscript Omg_annotation.R my_cluster_markers.csv 
```

## Requirements

- R (>= 4.0)
- R packages: `tidyverse`, `reshape2`, `pheatmap`

## Input

### User-provided: marker gene CSV

A CSV file with cluster markers (e.g., from Seurat `FindAllMarkers()`):

| Column | Description |
|--------|-------------|
| `gene` | Gene name or ID |
| `cluster` | Cluster ID |
| `avg_log2FC` | Log2 fold change |

Gene names with hyphens vs. underscores (e.g., `Gene-01` vs `Gene_01`) are handled automatically.

### Reference data (included)

| File | Description |
|------|-------------|
| `orthogroups.tsv.zip` | Orthogroup definitions across 34 species from OrthoFinder (unzip before running) |
| `reference_markers.tsv.zip` | Reference marker genes from 17 plant species, 41 tissue types (unzip before running) |

## Usage

```bash
Rscript Omg_annotation.R <marker_gene.csv> [FDR_threshold] [top_n_genes]
```

| Argument | Required | Default | Description |
|----------|----------|---------|-------------|
| `marker_gene.csv` | Yes | - | Path to marker gene CSV |
| `FDR_threshold` | No | 0.01 | FDR cutoff for significance |
| `top_n_genes` | No | 200 | Top N marker genes per cluster (ranked by avg_log2FC) |

## Output

Results are saved to a timestamped directory `output_YYYYMMDD_HHMMSS/`:

| File | Description |
|------|-------------|
| `cell_type_predictions.csv` | **Main result.** Cell type prediction with confidence score per cluster |
| `compare_15species_all.csv` | Multi-species pooled comparison (p-values and shared orthogroups) |
| `compare_15species_heatmap.pdf` | Heatmap visualization of multi-species comparison |
| `pairwise/` | Per species-tissue results (significance tables, gene lists, heatmaps) |

## OMG Annotation Guide

### Step 1: Start with `cell_type_predictions.csv`

This is the main output. Each row contains:
- **consolidated_cell_type** — the predicted cell type
- **prediction_confidence** — proportion of evidence supporting the prediction (0 to 1)
- **cell_type_prediction_frequency** — all matched reference cell types with hit counts

Use the **prediction confidence** to guide how much to trust the result:

| Confidence | Level | Action |
|-----------|-------|--------|
| 0.7 -- 1.0 | High | Use directly. One cell type clearly dominates across species. |
| 0.4 -- 0.7 | Moderate | Likely correct. Confirm with the heatmap (Step 2). |
| 0.2 -- 0.4 | Low | Multiple competing cell types. Inspect `compare_15species_all.csv` (Step 3). |
| < 0.2 | Very low | Too ambiguous for automated prediction. Manual annotation recommended. |

Low confidence does not necessarily mean the prediction is wrong — it often indicates transitional cells, mixed identity clusters, or cell types with few matches in the reference panel.

### Step 2: Inspect the multi-species heatmap

Open `compare_15species_heatmap.pdf` for a visual overview of orthogroup overlap significance between your query clusters and all reference cell types. Look for strong, consistent dark cells to confirm or refine ambiguous predictions.

### Step 3: Dig into `compare_15species_all.csv`

For low-confidence clusters, open this file and filter by your query cluster.

Sort by:

1. **p_value** (ascending) — most significant matches first

2. **number_OMGs** (descending) — most shared orthogroups first

The combination of low p-value and high number of common orthogroups gives the strongest evidence for a cell type assignment.

### Step 4: Check pairwise results

The `pairwise/` folder contains per species-tissue comparisons. Use these to see which specific species or tissue drives a match, or to inspect shared marker genes for a particular pair.

> **Tip:** *Arabidopsis thaliana* root provides the most informative cell type predictions due to its well-characterized cell type atlas and extensive data availability. When resolving ambiguous clusters, the Arabidopsis root pairwise result is often the most reliable reference to check first.

### Adjusting parameters

If predictions have **too many competing cell types** (low confidence, noisy results), try lowering the FDR threshold to increase stringency:

```bash
Rscript Omg_annotation.R markers.csv 0.001 200   # stricter FDR
```

If certain clusters have **no significant match**, try relaxing the FDR threshold and/or increasing the number of marker genes:

```bash
Rscript Omg_annotation.R markers.csv 0.05 500    # more lenient FDR + more genes
```

<details>
<summary><strong>How consolidated_cell_type and prediction_confidence are calculated (click to expand)</strong></summary>

### Consolidated Cell Type

Each significant reference cell type match is mapped to a **broad group** (e.g., Root hair + Trichoblast → "Root hair cell"; Companion cell + Sieve element + Protophloem → "Phloem"). The full mapping of ~120 reference cell types to broad groups is defined in `cell_type_groups` inside `Omg_annotation.R`.

The broad group with the most hits becomes the prediction. When there are ties, the following rules apply in order:

1. **Epidermis reassignment**: If Epidermis hits coexist with Root hair cell or Non-hair cell hits, Epidermis is reassigned to whichever subtype has more hits. If Cortex/Exodermis/Endodermis outnumber Epidermis, Epidermis is reassigned to the dominant ground tissue type.

2. **Major tissue layer tie-breaking**: Broad groups belong to one of three root tissue layers:
   - **Outer**: Root cap, Epidermis, Root hair cell, Non-hair cell, Trichome
   - **Middle**: Exodermis, Cortex, Endodermis
   - **Inner**: Pericycle, Stele, Phloem, Xylem, Vascular

   When multiple broad groups are tied, only those in the dominant layer (the layer with the most total hits) are kept.

3. **Cell cycle override**: S phase, G2/M phase, and Proliferating cell are all dividing meristematic cells — S phase and G2/M phase are simply more specific labels. If the top broad group is Meristematic but G2/M phase or S phase hits are at least half the count of the top individual meristematic cell type, the prediction is overridden to the dominant cell cycle phase.

4. **Phloem refinement**: If the prediction is Phloem but Sieve element or Companion cell hits outnumber generic Phloem hits, the prediction is refined to the more specific subtype.

### Prediction Confidence

The confidence score measures how much of the evidence supports the predicted cell type:

```
confidence = hits for predicted broad group / total significant hits
```

**Key details:**

- **Single winner** (e.g., "Xylem"): all hits mapping to the Xylem broad group are counted as supporting. Example: Xylem (9) + Metaxylem (2) + Protoxylem (1) = 12 supporting out of 15 total → confidence = 0.80.

- **Tied predictions** (e.g., "Cortex/Endodermis/Exodermis"): confidence reflects **one** group's count, not the combined total. Example: each group has 1 hit out of 5 total → confidence = 0.20, not 0.60. This accurately represents the ambiguity.

- **Cell cycle phases** (S phase, G2/M phase): since these are specific labels for meristematic cells, all Meristematic + Cell cycle hits are counted together as supporting evidence. Example: Meristematic (16) + Cell cycle (7) = 23 supporting out of 53 total → confidence = 0.54 for an "S phase" prediction.

- **Phloem subtypes** (Sieve element, Companion cell): all Phloem-group hits count as supporting evidence. Example: Companion cell (8) + Phloem (8) + Pholem (1) = 16 supporting out of 19 total → confidence = 0.84 for a "Companion cell" prediction.

</details>

## Reference Species

| Species | Tissues |
|---------|---------|
| *Arabidopsis thaliana* | Flower, Inflorescence, Leaf, Root, Seed, Shoot axis apex, Stem |
| *Oryza sativa* | Inflorescence, Leaf, Pistil, Root |
| *Zea mays* | Ear inflorescence, Inflorescence, Leaf, Root, Shoot axis apex |
| *Glycine max* | Flower, Leaf, Root, Seed, Stem |
| *Solanum lycopersicum* | Leaf, Root, Shoot axis apex |
| *Gossypium hirsutum* | Flower, Leaf |
| *Gossypium bickii* | Leaf, Seed |
| *Manihot esculenta* | Leaf, Root, Tuberous root |
| *Nicotiana attenuata* | Corolla, Flower |
| *Brassica rapa* | Leaf |
| *Catharanthus roseus* | Leaf |
| *Fragaria vesca* | Leaf |
| *Medicago truncatula* | Root |
| *Populus alba* var. *pyramidalis* | Stem |
| *Sorghum bicolor* | Root |
| *Triticum aestivum* | Root |

## Statistical Method

For each query cluster vs. reference cell type pair:

1. Map marker genes to orthogroups using the OrthoFinder orthogroup table
2. Count overlapping orthogroups between query cluster and reference cell type
3. Construct a 2x2 contingency table of shared vs. unique orthogroups
4. Perform one-sided Fisher's exact test (alternative = "greater")
5. Apply Benjamini-Hochberg FDR correction across all comparisons

## Citation

If you use this pipeline, please cite:

> Chau, T.N., Timilsena, P.R., Bathala, S.P., Kundu, S., Bargmann, B.O.R. & Li, S. Orthologous marker groups reveal broad cell identity conservation across plant single-cell transcriptomes. *Nature Communications* 16, 201 (2025). https://doi.org/10.1038/s41467-024-55755-0

## Repository Structure

```
.
├── Omg_annotation.R              # Main pipeline script
├── orthogroups.tsv.zip           # Orthogroup definitions (34 species) — unzip before running
├── reference_markers.tsv.zip     # Reference marker genes (17 species, 41 tissues) — unzip before running
└── README.md
```
