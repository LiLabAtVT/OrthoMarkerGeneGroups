#' Default cell-type grouping map
#'
#' Maps the reference cell-type labels to broad groups used when consolidating
#' predictions. Pass a modified copy to \code{\link{omg}} via
#' \code{cell_type_groups=} to extend it for new species' labels.
#'
#' @return A named list: group name -> character vector of reference cell types.
#' @export
omg_cell_type_groups <- function() {
  list(
    "Epidermis"      = c("Epidermis", "Leaf epidermis", "Shoot system epidermis",
                         "Abaxial epidermis", "Adaxial epidermis",
                         "Leaf pavement cell", "Pavement cell-A", "Pavement cell-N",
                         "Leaf guard cell", "Guard cell", "Leaf subsidiary cell",
                         "Meristem epidermis", "Outer cell layer",
                         "Cork cambium", "Root epidermis",
                         "Root epidermis/Root cortex"),
    "Root hair cell" = c("Root hair", "Trichoblast"),
    "Non-hair cell"  = c("Non-hair", "Atrichoblast", "Non-hair root epidermal cell"),
    "Trichome"       = c("Trichome", "Idioblast cell", "Pigment gland"),
    "Mesophyll"      = c("Mesophyll", "Palisade mesophyll", "Spongy mesophyll",
                         "Chlorenchyma", "Photosynthetic cell",
                         "Mesophyll/Cortex", "Bundle sheath", "Mestome sheath"),
    "Phloem"         = c("Phloem", "Companion cell", "Sieve element",
                         "Protophloem", "Phloem parenchyma",
                         "Phloem pole pericycle", "Phloem/Pericycle",
                         "Phloem_Procamb"),
    "Xylem"          = c("Xylem", "Metaxylem", "Protoxylem",
                         "Xylem parenchyma", "Xylem precursor cell",
                         "Xylem pole pericycle", "Fiber cell", "Sclerenchyma"),
    "Vascular"       = c("Vascular tissue", "Vascular cambium", "Vascular bundle",
                         "Procambium", "Root procambium", "Provascular",
                         "Explant vasculature and callus founder cell",
                         "Hydathodes"),
    "Cortex"         = c("Cortex", "Root cortex", "Inner cortex", "Outer cortex",
                         "Cortex/Endodermis", "Cortex/Endodermis initial",
                         "Shoot system cortex", "Ear cortex"),
    "Endodermis"     = c("Endodermis", "Root endodermis",
                         "Shoot system endodermis", "Mature endodermis",
                         "Mesophyll precursor/Root endodermis"),
    "Exodermis"      = c("Exodermis", "Root exodermis",
                         "Root exodermis/Sclerenchym"),
    "Pericycle"      = c("Pericycle", "Pericycle cell",
                         "Xylem pole pericycle cell", "Phloem pole pericycle"),
    "Stele"          = c("Root stele", "Mature stele"),
    "Root cap"       = c("Root cap", "Columella root cap", "Lateral root cap",
                         "Lateral root cap like cell", "Root cap junction",
                         "Columella root cap cell",
                         "RC_Col_QC", "Root broder"),
    "Root"           = c("Root pith", "Lateral root primordia"),
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
    "Cell cycle"     = c("S phase", "S-phase cell", "G2/M phase", "G2/M-phase cell", "G1/G0 phase"),
    "Pollen/Spore"   = c("Pollen", "Sperm nuclei", "Generative nuclei",
                         "Microspore nuclei", "Vegetative nuclei"),
    "Parenchyma"     = c("Parenchyma", "Pith", "Ray parenchyma",
                         "Inner cell layer"),
    "Inflorescence"  = c("Inflorescence axis", "Spikelet",
                         "Spikelet floret", "Spikelet meristem"),
    "Ovule/Seed"     = c("Ovule", "CT-Ovule", "S-Ovule", "Nucellus",
                         "Integument", "Inner ovary wall/Outer integument",
                         "Outer ovary wall", "Ovary wall"),
    "Flower"         = c("Stigma", "Style", "Stamen residue"),
    "Nodule"         = c("Infected cell", "Uninfected cell"),
    "Other"          = c("Contaminating nuclei", "Contamination",
                         "Stress response", "Outer pigment layer",
                         "Middle cell layer")
  )
}

# Consolidate significant hits into one prediction per query cluster.
# `sig_hits` must have columns: query_cluster, reference_cluster (characters).
consolidate_cell_types <- function(sig_hits, cell_type_groups = omg_cell_type_groups()) {

  get_broad_group <- function(ct) {
    for (group_name in names(cell_type_groups)) {
      if (ct %in% cell_type_groups[[group_name]]) return(group_name)
    }
    ct
  }

  major_layer_map <- c(
    "Root cap" = "Outer", "Epidermis" = "Outer",
    "Root hair cell" = "Outer", "Non-hair cell" = "Outer",
    "Trichome" = "Outer",
    "Exodermis" = "Middle", "Cortex" = "Middle",
    "Endodermis" = "Middle",
    "Pericycle" = "Inner", "Stele" = "Inner",
    "Phloem" = "Inner", "Xylem" = "Inner",
    "Vascular" = "Inner"
  )

  predictions_summary <- sig_hits %>%
    dplyr::group_by(query_cluster) %>%
    dplyr::summarise(
      cell_type_prediction_frequency = {
        ct_freq <- sort(table(reference_cluster), decreasing = TRUE)
        paste(paste0(names(ct_freq), " (", ct_freq, ")"), collapse = "; ")
      },
      consolidated_cell_type = {
        broad <- sapply(reference_cluster, get_broad_group)

        rh_count <- sum(broad == "Root hair cell")
        nh_count <- sum(broad == "Non-hair cell")
        if (sum(broad == "Epidermis") > 0 && (rh_count > 0 || nh_count > 0)) {
          if (rh_count > nh_count) {
            broad[broad == "Epidermis"] <- "Root hair cell"
          } else if (nh_count > rh_count) {
            broad[broad == "Epidermis"] <- "Non-hair cell"
          }
        }

        epi_count <- sum(broad == "Epidermis")
        ground_counts <- c(Cortex = sum(broad == "Cortex"),
                           Exodermis = sum(broad == "Exodermis"),
                           Endodermis = sum(broad == "Endodermis"))
        total_ground <- sum(ground_counts)
        if (epi_count > 0 && total_ground > epi_count) {
          top_ground <- names(which.max(ground_counts))
          broad[broad == "Epidermis"] <- top_ground
        }

        major <- major_layer_map[broad]
        layer_counts <- table(major[!is.na(major)])

        group_freq <- sort(table(broad), decreasing = TRUE)
        top_count <- group_freq[1]
        tied <- names(group_freq[group_freq == top_count])

        if (length(tied) > 1 && length(layer_counts) > 0) {
          top_layer <- names(sort(layer_counts, decreasing = TRUE))[1]
          tied_in_layer <- tied[tied %in% names(major_layer_map[major_layer_map == top_layer])]
          if (length(tied_in_layer) > 0) tied <- tied_in_layer
        }

        if ("Meristematic" %in% tied) {
          g2m_count <- sum(reference_cluster == "G2/M phase" |
                             reference_cluster == "G2/M-phase cell")
          s_count   <- sum(reference_cluster == "S phase" |
                             reference_cluster == "S-phase cell")
          merist_types <- cell_type_groups[["Meristematic"]]
          merist_individual_counts <- table(reference_cluster[reference_cluster %in% merist_types])
          top_merist_individual <- if (length(merist_individual_counts) > 0) max(merist_individual_counts) else 0
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

        count_groups <- tied
        if (any(tied %in% c("S phase", "G2/M phase", "Meristematic"))) {
          count_groups <- unique(c(count_groups, "Meristematic", "Cell cycle"))
        }
        if (any(tied %in% c("Sieve element", "Companion cell"))) {
          count_groups <- unique(c(count_groups, "Phloem"))
        }
        total_hits <- length(broad)
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
    tidyr::separate(consolidated_cell_type,
                    into = c("consolidated_cell_type", "prediction_confidence"),
                    sep = "\\|\\|\\|", convert = TRUE) %>%
    dplyr::arrange(suppressWarnings(as.numeric(as.character(query_cluster))))

  predictions_summary
}
