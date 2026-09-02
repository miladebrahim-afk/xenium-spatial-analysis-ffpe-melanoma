# ==============================================================================
# 03_cd8_abundance_and_spatial.R
#
# Figure 3C and source spatial panels for Figure 3B/E.
#
# Figure 3C:
#   - Cytotoxic CD8+ fraction / total annotated cells
#   - CXCL13+ CD8+ fraction / total annotated cells
#   - Treg fraction / total annotated cells
#   - Treg / total CD8 ratio, where total CD8 =
#       CD8+ T Cytotoxic + CD8+ T(CXCL13)Helper
#   - tissue-core-level two-sided Wilcoxon comparisons by NF1 genotype
#
# Spatial outputs:
#   - major-lineage maps for all 39 analyzed cores
#   - CD8-focused maps for all 39 analyzed cores
#
# NOTE:
# `patient_id` in the saved Xenium objects behaves operationally as a
# tissue-core/specimen identifier (39 unique analyzed IDs). Outputs therefore
# use `sample_id` terminology.
#
# Input:
#   data/processed/xenium_list_with_meta_niche_published_reference.rds
#
# Outputs:
#   results/spatial_cd8/cd8_abundance/
#   results/spatial_cd8/spatial_maps/
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
})

# ------------------------------------------------------------------------------
# 0. Paths
# ------------------------------------------------------------------------------

xlist_file <- file.path(
  "data",
  "processed",
  "xenium_list_with_meta_niche_published_reference.rds"
)

abundance_dir <- file.path(
  "results",
  "spatial_cd8",
  "cd8_abundance"
)

spatial_root <- file.path(
  "results",
  "spatial_cd8",
  "spatial_maps"
)

major_dir <- file.path(
  spatial_root,
  "major_lineages"
)

cd8_dir <- file.path(
  spatial_root,
  "cd8_focus"
)

dir.create(
  abundance_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  major_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  cd8_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# ------------------------------------------------------------------------------
# 1. Constants
# ------------------------------------------------------------------------------

CYTOTOXIC <- "CD8+ T Cytotoxic"
HELPER <- "CD8+ T(CXCL13)Helper"
TREG <- "Tregs"

MAJOR_LINEAGE_MAP <- c(
  "Melanoma Melanocytic" = "Melanoma",
  "Melanoma Proliferative" = "Melanoma",
  "Melanoma Intermediate" = "Melanoma",
  "Melanoma NC" = "Melanoma",
  "Melanoma Mesenchymal" = "Melanoma",

  "CD8+ T Cytotoxic" = "T cells",
  "CD8+ T(CXCL13)Helper" = "T cells",
  "Tregs" = "T cells",

  "Macrophages M2" = "Myeloid",
  "Macrophages IFN-γ–activated" = "Myeloid",
  "Activated Myeloid cells" = "Myeloid",

  "B cells" = "B cells",
  "Plasma cells" = "Plasma cells",

  "CAFS Inflammatory" = "Stromal",
  "CAFS Myofibroblast" = "Stromal",
  "Pericytes" = "Stromal",
  "Endothelial cells" = "Stromal",

  "Epithelial" = "Epithelial"
)

MAJOR_COLORS <- c(
  "Melanoma" = "#E41A1C",
  "T cells" = "#377EB8",
  "Myeloid" = "#4DAF4A",
  "B cells" = "#984EA3",
  "Plasma cells" = "#FF7F00",
  "Stromal" = "#FFBF00",
  "Epithelial" = "#A65628"
)

CD8_FOCUS_COLORS <- c(
  "Other" = "#D9D9D9",
  "CD8+ T Cytotoxic" = "#2166AC",
  "CD8+ T(CXCL13)Helper" = "#B2182B",
  "Tregs" = "#762A83"
)

NF1_COLORS <- c(
  "Mut" = "#D73027",
  "WT" = "#4575B4"
)

# ------------------------------------------------------------------------------
# 2. Helpers
# ------------------------------------------------------------------------------

normalize_nf1 <- function(x) {
  x <- as.character(x)

  out <- dplyr::case_when(
    x %in% c("Mut", "NF1Mut", "NF1_Mut", "NF1 Mut") ~ "Mut",
    x %in% c("WT", "NF1WT", "NF1_WT", "NF1 WT") ~ "WT",
    TRUE ~ NA_character_
  )

  if (anyNA(out)) {
    stop(
      "Unexpected NF1 value(s): ",
      paste(
        sort(unique(x[is.na(out)])),
        collapse = ", "
      )
    )
  }

  out
}

single_value <- function(x, field, i) {
  x <- unique(
    as.character(
      x[!is.na(x)]
    )
  )

  if (length(x) != 1L) {
    stop(
      "Expected one unique ",
      field,
      " in sample_index ",
      i,
      "; found: ",
      paste(x, collapse = ", ")
    )
  }

  x
}

safe_wilcox <- function(x, g) {
  suppressWarnings(
    tryCatch(
      wilcox.test(
        x ~ factor(g),
        exact = FALSE
      )$p.value,
      error = function(e) NA_real_
    )
  )
}

format_p <- function(p) {
  if (is.na(p)) return("p = NA")
  if (p < 0.001) return("p < 0.001")
  paste0(
    "p = ",
    formatC(
      p,
      format = "f",
      digits = 3
    )
  )
}

sanitize_filename <- function(x) {
  x <- gsub(
    "[^A-Za-z0-9._-]+",
    "_",
    x
  )
  x <- gsub(
    "_+",
    "_",
    x
  )
  gsub(
    "^_|_$",
    "",
    x
  )
}

# ------------------------------------------------------------------------------
# 3. Load analyzed cores
# ------------------------------------------------------------------------------

if (!file.exists(xlist_file)) {
  stop(
    "Missing input: ",
    xlist_file
  )
}

cat(
  "Loading publication-reference Xenium tissue cores...\n"
)

xenium_list <- readRDS(
  xlist_file
)

if (length(xenium_list) != 39L) {
  stop(
    "Expected 39 analyzed cores; found ",
    length(xenium_list),
    "."
  )
}

# ------------------------------------------------------------------------------
# 4. Per-core CD8/Treg composition
#
# Denominator is all cells with a retained biological `celltype` annotation.
# This matches the filtered integrated cell set used for the historical
# composition analysis, while the per-core objects additionally retain filtered
# geometry-only cells with celltype = NA.
# ------------------------------------------------------------------------------

cat(
  "\nCalculating per-core CD8/Treg abundance...\n"
)

rows <- vector(
  "list",
  length(xenium_list)
)

for (i in seq_along(xenium_list)) {

  obj <- xenium_list[[i]]
  md <- obj@meta.data

  required <- c(
    "patient_id",
    "NF1",
    "celltype"
  )

  missing_cols <- setdiff(
    required,
    colnames(md)
  )

  if (length(missing_cols) > 0L) {
    stop(
      "Core ",
      i,
      " missing: ",
      paste(
        missing_cols,
        collapse = ", "
      )
    )
  }

  sample_id <- single_value(
    md$patient_id,
    "patient_id",
    i
  )

  nf1 <- single_value(
    normalize_nf1(md$NF1),
    "NF1",
    i
  )

  ct <- as.character(
    md$celltype
  )

  retained <- !is.na(ct)

  n_total <- sum(retained)

  n_cytotoxic <- sum(
    ct == CYTOTOXIC,
    na.rm = TRUE
  )

  n_helper <- sum(
    ct == HELPER,
    na.rm = TRUE
  )

  n_treg <- sum(
    ct == TREG,
    na.rm = TRUE
  )

  n_total_cd8 <-
    n_cytotoxic +
    n_helper

  rows[[i]] <- data.frame(
    sample_index = i,
    sample_id = sample_id,
    NF1 = nf1,

    n_total_annotated = n_total,
    n_cd8_cytotoxic = n_cytotoxic,
    n_cd8_cxcl13 = n_helper,
    n_treg = n_treg,
    n_total_cd8 = n_total_cd8,

    frac_cd8_cytotoxic =
      n_cytotoxic / n_total,

    frac_cd8_cxcl13 =
      n_helper / n_total,

    frac_treg =
      n_treg / n_total,

    treg_total_cd8_ratio =
      if (
        n_total_cd8 > 0
      ) {
        n_treg / n_total_cd8
      } else {
        NA_real_
      },

    stringsAsFactors = FALSE
  )
}

core_comp <- bind_rows(
  rows
) %>%
  mutate(
    NF1 = factor(
      NF1,
      levels = c(
        "Mut",
        "WT"
      )
    )
  )

write.csv(
  core_comp,
  file.path(
    abundance_dir,
    "figure3C_core_level_CD8_Treg_metrics.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------------------------
# 5. Figure 3C Wilcoxon statistics
# ------------------------------------------------------------------------------

metric_info <- data.frame(
  metric = c(
    "frac_cd8_cytotoxic",
    "frac_cd8_cxcl13",
    "frac_treg",
    "treg_total_cd8_ratio"
  ),
  label = c(
    "Cytotoxic CD8+ / total cells",
    "CXCL13+ CD8+ / total cells",
    "Tregs / total cells",
    "Treg / total CD8 ratio"
  ),
  stringsAsFactors = FALSE
)

stats_list <- lapply(
  seq_len(
    nrow(metric_info)
  ),
  function(j) {

    metric <- metric_info$metric[j]

    d <- core_comp %>%
      filter(
        !is.na(
          .data[[metric]]
        )
      )

    data.frame(
      metric = metric,
      label = metric_info$label[j],

      n_mut = sum(
        d$NF1 == "Mut"
      ),

      n_wt = sum(
        d$NF1 == "WT"
      ),

      median_mut = median(
        d[[metric]][
          d$NF1 == "Mut"
        ],
        na.rm = TRUE
      ),

      median_wt = median(
        d[[metric]][
          d$NF1 == "WT"
        ],
        na.rm = TRUE
      ),

      mean_mut = mean(
        d[[metric]][
          d$NF1 == "Mut"
        ],
        na.rm = TRUE
      ),

      mean_wt = mean(
        d[[metric]][
          d$NF1 == "WT"
        ],
        na.rm = TRUE
      ),

      wilcox_p = safe_wilcox(
        d[[metric]],
        d$NF1
      ),

      stringsAsFactors = FALSE
    )
  }
)

figure3c_stats <- bind_rows(
  stats_list
)

write.csv(
  figure3c_stats,
  file.path(
    abundance_dir,
    "figure3C_wilcoxon_stats.csv"
  ),
  row.names = FALSE
)

cat(
  "\nFigure 3C Wilcoxon statistics:\n"
)

print(
  figure3c_stats,
  row.names = FALSE
)

# ------------------------------------------------------------------------------
# 6. Figure 3C plot
# ------------------------------------------------------------------------------

plot_df <- core_comp %>%
  select(
    sample_index,
    sample_id,
    NF1,
    all_of(
      metric_info$metric
    )
  ) %>%
  pivot_longer(
    cols = all_of(
      metric_info$metric
    ),
    names_to = "metric",
    values_to = "value"
  ) %>%
  left_join(
    metric_info,
    by = "metric"
  ) %>%
  mutate(
    label = factor(
      label,
      levels = metric_info$label
    )
  )

plot_ann <- figure3c_stats %>%
  select(
    label,
    wilcox_p
  )

p3c <- ggplot(
  plot_df,
  aes(
    x = NF1,
    y = value,
    fill = NF1
  )
) +
  geom_boxplot(
    width = 0.55,
    outlier.shape = NA,
    alpha = 0.65
  ) +
  geom_jitter(
    width = 0.12,
    size = 1.7,
    alpha = 0.75
  ) +
  facet_wrap(
    ~ label,
    scales = "free_y",
    nrow = 1
  ) +
  scale_fill_manual(
    values = NF1_COLORS
  ) +
  labs(
    x = "NF1 status",
    y = "Fraction / ratio"
  ) +
  theme_classic(
    base_size = 12
  ) +
  theme(
    legend.position = "none",
    strip.background = element_blank(),
    strip.text = element_text(
      face = "bold"
    ),
    axis.text.x = element_text(
      face = "bold"
    )
  )

# Add p-value labels separately to each facet.
panel_y <- plot_df %>%
  group_by(
    label
  ) %>%
  summarise(
    ymax = max(
      value,
      na.rm = TRUE
    ),
    ymin = min(
      value,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  mutate(
    ypos =
      ymax +
      0.10 *
      pmax(
        ymax - ymin,
        abs(ymax) * 0.10,
        1e-6
      )
  ) %>%
  left_join(
    plot_ann,
    by = "label"
  ) %>%
  mutate(
    p_label = vapply(
      wilcox_p,
      format_p,
      character(1)
    )
  )

p3c <- p3c +
  geom_text(
    data = panel_y,
    aes(
      x = 1.5,
      y = ypos,
      label = p_label
    ),
    inherit.aes = FALSE,
    size = 3.2
  )

ggsave(
  file.path(
    abundance_dir,
    "Figure3C_CD8_Treg_abundance.pdf"
  ),
  p3c,
  width = 13,
  height = 4.2,
  useDingbats = FALSE
)

ggsave(
  file.path(
    abundance_dir,
    "Figure3C_CD8_Treg_abundance.png"
  ),
  p3c,
  width = 13,
  height = 4.2,
  dpi = 300
)

# ------------------------------------------------------------------------------
# 7. Spatial source plots for Figure 3B/E
#
# Save all 39 rather than hard-coding a representative core.
# This keeps representative-panel selection transparent and non-circular.
# ------------------------------------------------------------------------------

cat(
  "\nGenerating major-lineage and CD8-focused spatial source plots...\n"
)

map_qc <- vector(
  "list",
  length(xenium_list)
)

for (i in seq_along(xenium_list)) {

  obj <- xenium_list[[i]]
  md <- obj@meta.data

  sample_id <- single_value(
    md$patient_id,
    "patient_id",
    i
  )

  nf1 <- single_value(
    normalize_nf1(md$NF1),
    "NF1",
    i
  )

  ct <- as.character(
    md$celltype
  )

  # ---------------------------------------------------------------------------
  # Major lineages
  # ---------------------------------------------------------------------------

  major <- unname(
    MAJOR_LINEAGE_MAP[ct]
  )

  major[
    is.na(major)
  ] <- NA_character_

  obj$major_celltype_plot <- factor(
    major,
    levels = names(
      MAJOR_COLORS
    )
  )

  # ---------------------------------------------------------------------------
  # CD8/Treg focused source plot
  # ---------------------------------------------------------------------------

  focus <- rep(
    "Other",
    length(ct)
  )

  focus[
    ct == CYTOTOXIC
  ] <- CYTOTOXIC

  focus[
    ct == HELPER
  ] <- HELPER

  focus[
    ct == TREG
  ] <- TREG

  focus[
    is.na(ct)
  ] <- NA_character_

  obj$cd8_focus_plot <- factor(
    focus,
    levels = names(
      CD8_FOCUS_COLORS
    )
  )

  stub <- sanitize_filename(
    paste0(
      sprintf(
        "core%02d",
        i
      ),
      "_",
      sample_id,
      "_NF1",
      nf1
    )
  )

  major_ok <- TRUE
  cd8_ok <- TRUE
  major_msg <- NA_character_
  cd8_msg <- NA_character_

  tryCatch(
    {
      p_major <- ImageDimPlot(
        obj,
        group.by =
          "major_celltype_plot",
        cols = MAJOR_COLORS,
        size = 1
      ) +
        ggtitle(
          paste0(
            sample_id,
            " | NF1 ",
            nf1,
            " | major lineages"
          )
        ) +
        theme(
          legend.position = "right",
          plot.title = element_text(
            face = "bold",
            hjust = 0.5
          )
        )

      ggsave(
        file.path(
          major_dir,
          paste0(
            stub,
            "_major_lineages.png"
          )
        ),
        p_major,
        width = 8,
        height = 6,
        dpi = 400,
        bg = "white"
      )
    },
    error = function(e) {
      major_ok <<- FALSE
      major_msg <<-
        conditionMessage(e)
    }
  )

  tryCatch(
    {
      p_cd8 <- ImageDimPlot(
        obj,
        group.by =
          "cd8_focus_plot",
        cols =
          CD8_FOCUS_COLORS,
        size = 1
      ) +
        ggtitle(
          paste0(
            sample_id,
            " | NF1 ",
            nf1,
            " | CD8/Treg"
          )
        ) +
        theme(
          legend.position = "right",
          plot.title = element_text(
            face = "bold",
            hjust = 0.5
          )
        )

      ggsave(
        file.path(
          cd8_dir,
          paste0(
            stub,
            "_CD8_Treg_focus.png"
          )
        ),
        p_cd8,
        width = 8,
        height = 6,
        dpi = 400,
        bg = "white"
      )
    },
    error = function(e) {
      cd8_ok <<- FALSE
      cd8_msg <<-
        conditionMessage(e)
    }
  )

  map_qc[[i]] <- data.frame(
    sample_index = i,
    sample_id = sample_id,
    NF1 = nf1,
    major_map_ok = major_ok,
    major_message = major_msg,
    cd8_map_ok = cd8_ok,
    cd8_message = cd8_msg,
    stringsAsFactors = FALSE
  )

  cat(
    sprintf(
      "  core %02d | %-10s | NF1 %-3s | major=%s | CD8=%s\n",
      i,
      sample_id,
      nf1,
      ifelse(
        major_ok,
        "OK",
        "FAIL"
      ),
      ifelse(
        cd8_ok,
        "OK",
        "FAIL"
      )
    )
  )

  rm(
    obj
  )

  gc(
    verbose = FALSE
  )
}

map_qc_df <- bind_rows(
  map_qc
)

write.csv(
  map_qc_df,
  file.path(
    spatial_root,
    "figure3B_E_spatial_map_QC.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------------------------
# 8. Final QC
# ------------------------------------------------------------------------------

cat(
  "\nFigure 3C direction summary:\n"
)

print(
  figure3c_stats %>%
    mutate(
      Mut_minus_WT_mean =
        mean_mut -
        mean_wt
    ) %>%
    select(
      label,
      mean_mut,
      mean_wt,
      Mut_minus_WT_mean,
      wilcox_p
    ),
  row.names = FALSE
)

cat(
  "\nSpatial-map QC:\n"
)

cat(
  "  major-lineage maps successful: ",
  sum(
    map_qc_df$major_map_ok
  ),
  "/39\n",
  sep = ""
)

cat(
  "  CD8-focused maps successful: ",
  sum(
    map_qc_df$cd8_map_ok
  ),
  "/39\n",
  sep = ""
)

cat(
  "\nDONE: Figure 3C and spatial source plots complete.\n"
)
