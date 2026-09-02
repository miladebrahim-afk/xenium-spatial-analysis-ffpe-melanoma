# ==============================================================================
# 01_distance_analysis.R
#
# Figure 3A: NF1Mut vs NF1WT differential spatial proximity.
# Uses per-core nearest-neighbor distances with FNN::get.knnx(k = 1).
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(FNN)
  library(ggplot2)
})

xlist_file <- file.path(
  "data", "processed",
  "xenium_list_with_meta_niche_published_reference.rds"
)

out_dir <- file.path(
  "results", "spatial_cd8", "distance"
)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

CELLTYPE_ORDER <- c(
  "Melanoma Melanocytic",
  "Melanoma Proliferative",
  "Melanoma Intermediate",
  "Melanoma NC",
  "Melanoma Mesenchymal",
  "CD8+ T Cytotoxic",
  "CD8+ T(CXCL13)Helper",
  "Tregs",
  "B cells",
  "Plasma cells",
  "Macrophages IFN-γ–activated",
  "Macrophages M2",
  "Activated Myeloid cells",
  "CAFS Inflammatory",
  "CAFS Myofibroblast",
  "Pericytes",
  "Endothelial cells",
  "Epithelial"
)

DISPLAY_LABEL <- c(
  "Melanoma Melanocytic" = "Melanocytic Mel",
  "Melanoma Proliferative" = "Proliferative Mel",
  "Melanoma Intermediate" = "Intermediate Mel",
  "Melanoma NC" = "NC-like Mel",
  "Melanoma Mesenchymal" = "Mesenchymal Mel",
  "CD8+ T Cytotoxic" = "Cytotoxic CD8+",
  "CD8+ T(CXCL13)Helper" = "CXCL13+ CD8+",
  "Tregs" = "Tregs",
  "B cells" = "B cells",
  "Plasma cells" = "Plasma cells",
  "Macrophages IFN-γ–activated" = "MΦ IFN-activated",
  "Macrophages M2" = "MΦ M2",
  "Activated Myeloid cells" = "Activated Myeloid",
  "CAFS Inflammatory" = "iCAFs",
  "CAFS Myofibroblast" = "myCAFs",
  "Pericytes" = "Pericytes",
  "Endothelial cells" = "Endothelial cells",
  "Epithelial" = "Epithelial"
)

MIN_CELLS_PER_GROUP <- 3L

normalize_nf1 <- function(x) {
  x <- as.character(x)
  out <- dplyr::case_when(
    x %in% c("Mut", "NF1Mut", "NF1_Mut", "NF1 Mut") ~ "Mut",
    x %in% c("WT", "NF1WT", "NF1_WT", "NF1 WT") ~ "WT",
    TRUE ~ NA_character_
  )
  if (anyNA(out)) {
    stop("Unexpected NF1 value(s): ",
         paste(sort(unique(x[is.na(out)])), collapse = ", "))
  }
  out
}

single_value <- function(x, field, sample_index) {
  x <- unique(as.character(x[!is.na(x)]))
  if (length(x) != 1L) {
    stop("Expected one ", field, " value in sample_index ", sample_index,
         "; found: ", paste(x, collapse = ", "))
  }
  x
}

safe_wilcox <- function(x, g) {
  g <- factor(g)
  if (length(unique(g[!is.na(g)])) != 2L) return(NA_real_)
  suppressWarnings(
    tryCatch(
      wilcox.test(x ~ g, exact = FALSE)$p.value,
      error = function(e) NA_real_
    )
  )
}

if (!file.exists(xlist_file)) stop("Missing input: ", xlist_file)

cat("Loading Xenium tissue-core objects...\n")
xenium_list <- readRDS(xlist_file)

if (length(xenium_list) != 39L) {
  stop("Expected 39 analyzed tissue cores; found ", length(xenium_list), ".")
}

cat("Loaded ", length(xenium_list), " tissue cores.\n", sep = "")
cat("\nComputing source-target nearest-neighbor distances...\n")

core_results <- vector("list", length(xenium_list))
core_qc <- vector("list", length(xenium_list))

for (i in seq_along(xenium_list)) {
  obj <- xenium_list[[i]]
  md <- obj@meta.data

  if (!all(c("patient_id", "NF1", "celltype") %in% colnames(md))) {
    stop("Missing required metadata in sample_index ", i)
  }

  sample_id <- single_value(md$patient_id, "patient_id", i)
  nf1 <- single_value(normalize_nf1(md$NF1), "NF1", i)

  coords <- GetTissueCoordinates(obj)

  if (!all(c("x", "y", "cell") %in% colnames(coords))) {
    stop("Coordinates in sample_index ", i, " must contain x, y, and cell.")
  }

  coords$celltype <- md$celltype[match(coords$cell, rownames(md))]

  coords <- coords %>%
    filter(
      !is.na(x),
      !is.na(y),
      !is.na(celltype),
      celltype %in% CELLTYPE_ORDER
    )

  coord_by_type <- lapply(
    CELLTYPE_ORDER,
    function(ct) {
      as.matrix(
        coords[
          coords$celltype == ct,
          c("x", "y"),
          drop = FALSE
        ]
      )
    }
  )
  names(coord_by_type) <- CELLTYPE_ORDER

  n_by_type <- vapply(coord_by_type, nrow, integer(1))

  pair_rows <- list()
  k <- 0L

  for (src in CELLTYPE_ORDER) {
    src_xy <- coord_by_type[[src]]
    if (nrow(src_xy) < MIN_CELLS_PER_GROUP) next

    for (trg in CELLTYPE_ORDER) {
      trg_xy <- coord_by_type[[trg]]
      if (nrow(trg_xy) < MIN_CELLS_PER_GROUP) next

      nn_dist <- FNN::get.knnx(
        data = trg_xy,
        query = src_xy,
        k = 1
      )$nn.dist[, 1]

      k <- k + 1L
      pair_rows[[k]] <- data.frame(
        sample_index = i,
        sample_id = sample_id,
        NF1 = nf1,
        source = src,
        target = trg,
        mean_distance = mean(nn_dist, na.rm = TRUE),
        median_distance = median(nn_dist, na.rm = TRUE),
        n_source = nrow(src_xy),
        n_target = nrow(trg_xy),
        stringsAsFactors = FALSE
      )
    }
  }

  core_results[[i]] <- bind_rows(pair_rows)

  core_qc[[i]] <- data.frame(
    sample_index = i,
    sample_id = sample_id,
    NF1 = nf1,
    annotated_cells_with_coordinates = nrow(coords),
    n_celltypes_present = sum(n_by_type > 0),
    n_celltypes_ge3 = sum(n_by_type >= MIN_CELLS_PER_GROUP),
    n_source_target_pairs = nrow(core_results[[i]]),
    stringsAsFactors = FALSE
  )

  cat(sprintf(
    "  core %02d | %-10s | NF1 %-3s | annotated=%7d | pairs=%3d\n",
    i, sample_id, nf1, nrow(coords), nrow(core_results[[i]])
  ))

  rm(obj, md, coords, coord_by_type, n_by_type, pair_rows)
  gc(verbose = FALSE)
}

distance_df <- bind_rows(core_results)
core_qc_df <- bind_rows(core_qc)

write.csv(
  distance_df,
  file.path(out_dir, "core_level_nearest_neighbor_distances.csv"),
  row.names = FALSE
)

write.csv(
  core_qc_df,
  file.path(out_dir, "distance_analysis_QC.csv"),
  row.names = FALSE
)

cat("\nCalculating NF1Mut vs NF1WT distance statistics...\n")

distance_stats <- distance_df %>%
  group_by(source, target) %>%
  summarise(
    n_mut_cores = sum(NF1 == "Mut"),
    n_wt_cores = sum(NF1 == "WT"),
    mean_distance_mut = mean(mean_distance[NF1 == "Mut"], na.rm = TRUE),
    mean_distance_wt = mean(mean_distance[NF1 == "WT"], na.rm = TRUE),
    median_core_distance_mut = median(mean_distance[NF1 == "Mut"], na.rm = TRUE),
    median_core_distance_wt = median(mean_distance[NF1 == "WT"], na.rm = TRUE),
    delta_mut_wt = mean_distance_mut - mean_distance_wt,
    wilcox_p = safe_wilcox(mean_distance, NF1),
    .groups = "drop"
  ) %>%
  mutate(
    # Published Figure 3A uses a relative distance effect, not an absolute
    # micrometer difference. Positive = farther in NF1Mut; negative = closer.
    # Contact effect for Figure 3A:
    # positive = closer / greater contact in NF1Mut
    # negative = farther / reduced contact in NF1Mut
    contact_effect = ifelse(
      source == target,
      0,
      log2(mean_distance_wt / mean_distance_mut)
    ),
    wilcox_p_BH = p.adjust(wilcox_p, method = "BH"),
    source_display = unname(DISPLAY_LABEL[source]),
    target_display = unname(DISPLAY_LABEL[target])
  )

write.csv(
  distance_stats,
  file.path(out_dir, "figure3A_distance_stats.csv"),
  row.names = FALSE
)

heatmap_df <- distance_stats %>%
  mutate(
    source_display = factor(
      source_display,
      levels = rev(unname(DISPLAY_LABEL[CELLTYPE_ORDER]))
    ),
    target_display = factor(
      target_display,
      levels = unname(DISPLAY_LABEL[CELLTYPE_ORDER])
    )
  )

# Final assembled Figure 3A uses a fixed relative-effect scale of -0.2 to 0.2.
# Values outside this range are saturated at the corresponding endpoint.
FIG3A_COLOR_LIMIT <- 0.2

p_heatmap <- ggplot(
  heatmap_df,
  aes(
    x = target_display,
    y = source_display,
    fill = contact_effect
  )
) +
  geom_tile(color = "white", linewidth = 0.15) +
  scale_fill_gradient2(
    # Contact effect:
    # negative = farther in NF1Mut = blue
    # positive = closer in NF1Mut = red
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0,
    limits = c(-FIG3A_COLOR_LIMIT, FIG3A_COLOR_LIMIT),
    oob = scales::squish,
    breaks = c(-0.2, -0.1, 0, 0.1, 0.2),
    name = expression(log[2]("mean distance NF1WT / NF1Mut"))
  ) +
  labs(
    x = "Target",
    y = "Source"
  ) +
  theme_classic(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1),
    axis.text.y = element_text(size = 9),
    axis.ticks = element_blank(),
    panel.border = element_rect(fill = NA, linewidth = 0.4)
  )

ggsave(
  file.path(out_dir, "Figure3A_differential_contact_heatmap.pdf"),
  p_heatmap,
  width = 11,
  height = 9,
  useDingbats = FALSE
)

ggsave(
  file.path(out_dir, "Figure3A_differential_contact_heatmap.png"),
  p_heatmap,
  width = 11,
  height = 9,
  dpi = 300
)

MELANOMA_STATES <- c(
  "Melanoma Melanocytic",
  "Melanoma Proliferative",
  "Melanoma Intermediate",
  "Melanoma NC",
  "Melanoma Mesenchymal"
)

KEY_TARGETS <- c(
  "CAFS Inflammatory",
  "CAFS Myofibroblast",
  "Tregs",
  "Macrophages M2",
  "Pericytes",
  "CD8+ T Cytotoxic",
  "CD8+ T(CXCL13)Helper",
  "Endothelial cells"
)

manuscript_qc <- distance_stats %>%
  filter(
    source %in% MELANOMA_STATES,
    target %in% KEY_TARGETS
  ) %>%
  arrange(source, target)

write.csv(
  manuscript_qc,
  file.path(out_dir, "figure3A_melanoma_key_target_direction_QC.csv"),
  row.names = FALSE
)

cat("\nMelanoma-state -> key-target distance changes (Mut - WT):\n")
print(
  manuscript_qc %>%
    select(
      source,
      target,
      delta_mut_wt,
      contact_effect,
      wilcox_p,
      wilcox_p_BH,
      n_mut_cores,
      n_wt_cores
    ),
  n = Inf,
  width = Inf
)

cat("\nDistance-analysis QC:\n")
cat("  tissue cores: ", nrow(core_qc_df), "\n", sep = "")
cat("  NF1Mut cores: ", sum(core_qc_df$NF1 == "Mut"), "\n", sep = "")
cat("  NF1WT cores: ", sum(core_qc_df$NF1 == "WT"), "\n", sep = "")
cat(
  "  directed core-level source-target measurements: ",
  format(nrow(distance_df), big.mark = ","),
  "\n",
  sep = ""
)

cat("\nDONE: Figure 3A spatial-distance analysis complete.\n")
