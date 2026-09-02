# ==============================================================================
# 03_mn_composition_and_nf1_abundance.R
#
# Meta-niche composition, tissue-core abundance, occurrence, and NF1 association
#
# IMPORTANT
# ---------
# - The neighborhood table contains 18 biological cell-type fraction columns
#   plus QC/metadata columns. ONLY the 18 explicit cell-type columns below are
#   used for MN composition.
# - By default, this script applies the documented 140-cell publication-reference
#   correction table in memory so downstream outputs reproduce the archived final
#   manuscript MN assignments. The de novo RDS is NEVER overwritten.
# - Set USE_PUBLISHED_REFERENCE <- FALSE to analyze the independent de novo
#   k-means reconstruction instead.
# - "patient_id" in these Xenium objects operationally identifies the analyzed
#   tissue core/specimen. The public outputs therefore use "sample_id" and
#   "n_cores" terminology where appropriate.
#
# Required inputs
# ---------------
# results/meta_niches/cell_neighborhoods_40um_k12.rds
# data/reference/published_meta_niche_boundary_corrections.csv
#   (required only when USE_PUBLISHED_REFERENCE = TRUE)
#
# Main outputs
# ------------
# results/meta_niches/meta_niche_composition_<mode>.csv
# results/meta_niches/meta_niche_frequency_per_core_<mode>.csv
# results/meta_niches/meta_niche_frequency_complete_per_core_<mode>.csv
# results/meta_niches/meta_niche_abundance_NF1_wilcoxon_<mode>.csv
# results/meta_niches/meta_niche_prevalence_NF1_fisher_<mode>.csv
# results/meta_niches/figure1H_core_occurrence_<mode>.csv
# results/meta_niches/figure2A_meta_niche_NF1_stats_<mode>.csv
# plus PDF/PNG QC/figure-style plots
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
})

# ------------------------------------------------------------------------------
# 0. Configuration
# ------------------------------------------------------------------------------

USE_PUBLISHED_REFERENCE <- TRUE

mn_dir <- file.path("results", "meta_niches")
dir.create(mn_dir, recursive = TRUE, showWarnings = FALSE)

neighborhood_rds <- file.path(
  mn_dir,
  "cell_neighborhoods_40um_k12.rds"
)

correction_file <- file.path(
  "data",
  "reference",
  "published_meta_niche_boundary_corrections.csv"
)

analysis_mode <- if (USE_PUBLISHED_REFERENCE) {
  "published_reference"
} else {
  "de_novo"
}

cat("Meta-niche analysis mode:", analysis_mode, "\n")

# Final 18-lineage order used for the Xenium analysis.
CELLTYPE_ORDER_FINAL <- c(
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

MN_LEVELS <- paste0("MetaNiche_", 1:12)

# Final manuscript/figure labels.
MN_SHORT_MAP <- c(
  "MetaNiche_1"  = "MN1: Inter+M2",
  "MetaNiche_2"  = "MN2: NC+M2+CAF",
  "MetaNiche_3"  = "MN3: IFN-myeloid+T cell",
  "MetaNiche_4"  = "MN4: Vasc+perivascular",
  "MetaNiche_5"  = "MN5: Prol+Melan",
  "MetaNiche_6"  = "MN6: iCAF+M2",
  "MetaNiche_7"  = "MN7: M2+myCAFs",
  "MetaNiche_8"  = "MN8: myCAF+Mes",
  "MetaNiche_9"  = "MN9: Inter+Prol+Melan",
  "MetaNiche_10" = "MN10: Melan+Prol+NC",
  "MetaNiche_11" = "MN11: Treg+B cell",
  "MetaNiche_12" = "MN12: Melan+M2"
)

# Figure 1H archived final occurrence counts.
EXPECTED_FIG1H <- c(
  "MetaNiche_1"  = 37L,
  "MetaNiche_2"  = 39L,
  "MetaNiche_3"  = 39L,
  "MetaNiche_4"  = 39L,
  "MetaNiche_5"  = 33L,
  "MetaNiche_6"  = 20L,
  "MetaNiche_7"  = 39L,
  "MetaNiche_8"  = 39L,
  "MetaNiche_9"  = 26L,
  "MetaNiche_10" = 30L,
  "MetaNiche_11" = 33L,
  "MetaNiche_12" = 36L
)

# ------------------------------------------------------------------------------
# 1. Helpers
# ------------------------------------------------------------------------------

same_label <- function(a, b) {
  (is.na(a) & is.na(b)) |
    (!is.na(a) & !is.na(b) & a == b)
}

normalize_nf1 <- function(x) {
  x <- as.character(x)

  out <- case_when(
    x %in% c("Mut", "NF1Mut", "NF1_Mut", "NF1 Mut") ~ "Mut",
    x %in% c("WT", "NF1WT", "NF1_WT", "NF1 WT") ~ "WT",
    TRUE ~ NA_character_
  )

  if (anyNA(out)) {
    bad <- sort(unique(x[is.na(out)]))
    stop(
      "Unexpected NF1 value(s): ",
      paste(bad, collapse = ", ")
    )
  }

  factor(out, levels = c("Mut", "WT"))
}

safe_wilcox <- function(x, g) {
  g <- droplevels(factor(g))
  if (nlevels(g) != 2L) return(NA_real_)

  tryCatch(
    wilcox.test(x ~ g, exact = FALSE)$p.value,
    error = function(e) NA_real_
  )
}

save_plot_both <- function(plot_obj, stem, width, height) {
  ggsave(
    filename = file.path(mn_dir, paste0(stem, ".pdf")),
    plot = plot_obj,
    width = width,
    height = height,
    units = "in",
    useDingbats = FALSE
  )

  ggsave(
    filename = file.path(mn_dir, paste0(stem, ".png")),
    plot = plot_obj,
    width = width,
    height = height,
    units = "in",
    dpi = 300
  )
}

# ------------------------------------------------------------------------------
# 2. Load exact 40-um neighborhood table
# ------------------------------------------------------------------------------

if (!file.exists(neighborhood_rds)) {
  stop("Missing neighborhood RDS: ", neighborhood_rds)
}

cat("\nLoading neighborhood table...\n")
cell_neighborhood_df <- readRDS(neighborhood_rds)

required_meta <- c(
  "cell",
  "patient_id",
  "NF1",
  "sample_index",
  "meta_niche"
)

missing_meta <- setdiff(required_meta, colnames(cell_neighborhood_df))
if (length(missing_meta) > 0L) {
  stop(
    "Neighborhood table is missing required columns: ",
    paste(missing_meta, collapse = ", ")
  )
}

if (nrow(cell_neighborhood_df) != 1085053L) {
  warning(
    "Expected 1,085,053 valid neighborhoods; found ",
    nrow(cell_neighborhood_df),
    "."
  )
}

# CRITICAL FIX:
# Never infer cell-type columns using setdiff(metadata columns), because the
# neighborhood table also contains QC fields such as n_neighbors_40um.
celltype_cols <- intersect(
  CELLTYPE_ORDER_FINAL,
  colnames(cell_neighborhood_df)
)

missing_celltypes <- setdiff(
  CELLTYPE_ORDER_FINAL,
  celltype_cols
)

if (length(missing_celltypes) > 0L) {
  stop(
    "Missing expected cell-type fraction columns: ",
    paste(missing_celltypes, collapse = ", ")
  )
}

if (length(celltype_cols) != 18L) {
  stop(
    "Expected exactly 18 biological cell-type columns; found ",
    length(celltype_cols),
    "."
  )
}

cat("Using exactly", length(celltype_cols), "biological cell-type columns.\n")

# Normalize NF1 labels only for downstream analysis; do not alter source RDS.
cell_neighborhood_df$NF1 <- normalize_nf1(
  cell_neighborhood_df$NF1
)

# ------------------------------------------------------------------------------
# 3. Optionally apply the documented publication-reference corrections
# ------------------------------------------------------------------------------

n_reference_changes <- 0L

if (USE_PUBLISHED_REFERENCE) {

  if (!file.exists(correction_file)) {
    stop(
      "USE_PUBLISHED_REFERENCE = TRUE but correction file is missing: ",
      correction_file
    )
  }

  cat("\nApplying publication-reference MN boundary corrections in memory...\n")
  corr <- read.csv(
    correction_file,
    stringsAsFactors = FALSE
  )

  required_corr <- c(
    "sample_index",
    "cell",
    "de_novo_meta_niche",
    "published_meta_niche"
  )

  if (!all(required_corr %in% colnames(corr))) {
    stop(
      "Correction table must contain: ",
      paste(required_corr, collapse = ", ")
    )
  }

  corr$sample_index <- as.integer(corr$sample_index)

  if (nrow(corr) != 140L) {
    stop(
      "Expected the documented 140 boundary corrections; found ",
      nrow(corr),
      "."
    )
  }

  if (anyDuplicated(
    paste(corr$sample_index, corr$cell, sep = "__")
  )) {
    stop("Duplicate sample_index + cell key in correction table.")
  }

  neighborhood_key <- paste(
    cell_neighborhood_df$sample_index,
    cell_neighborhood_df$cell,
    sep = "__"
  )

  correction_key <- paste(
    corr$sample_index,
    corr$cell,
    sep = "__"
  )

  if (anyDuplicated(neighborhood_key)) {
    stop("Neighborhood sample_index + cell key is not unique.")
  }

  pos <- match(
    correction_key,
    neighborhood_key
  )

  if (anyNA(pos)) {
    stop(
      sum(is.na(pos)),
      " publication-reference cells were not found in the valid neighborhood table."
    )
  }

  current_labels <- as.character(
    cell_neighborhood_df$meta_niche[pos]
  )

  expected_de_novo <- corr$de_novo_meta_niche

  if (!all(same_label(current_labels, expected_de_novo))) {
    bad <- which(
      !same_label(current_labels, expected_de_novo)
    )

    stop(
      "Publication-reference validation failed for ",
      length(bad),
      " row(s). The correction table does not match this de novo neighborhood RDS."
    )
  }

  cell_neighborhood_df$meta_niche[pos] <-
    corr$published_meta_niche

  n_reference_changes <- length(pos)

  cat(
    "Applied",
    n_reference_changes,
    "documented boundary-cell corrections.\n"
  )

  rm(
    corr,
    neighborhood_key,
    correction_key,
    pos,
    current_labels,
    expected_de_novo
  )
  gc()
}

# Validate MN values.
observed_mn <- sort(
  unique(
    as.character(
      cell_neighborhood_df$meta_niche[
        !is.na(cell_neighborhood_df$meta_niche)
      ]
    )
  )
)

unexpected_mn <- setdiff(
  observed_mn,
  MN_LEVELS
)

if (length(unexpected_mn) > 0L) {
  stop(
    "Unexpected meta-niche label(s): ",
    paste(unexpected_mn, collapse = ", ")
  )
}

# ------------------------------------------------------------------------------
# 4. Core/specimen QC
# ------------------------------------------------------------------------------

sample_info <- cell_neighborhood_df %>%
  distinct(
    sample_index,
    patient_id,
    NF1
  ) %>%
  mutate(
    sample_id = as.character(patient_id)
  ) %>%
  arrange(sample_index)

if (nrow(sample_info) != 39L) {
  warning(
    "Expected 39 analyzed Xenium tissue cores; found ",
    nrow(sample_info),
    "."
  )
}

nf1_core_counts <- sample_info %>%
  count(
    NF1,
    name = "n_cores"
  )

cat("\nAnalyzed tissue cores by NF1 status:\n")
print(nf1_core_counts, row.names = FALSE)

write.csv(
  sample_info,
  file.path(
    mn_dir,
    paste0("meta_niche_sample_info_", analysis_mode, ".csv")
  ),
  row.names = FALSE
)

# ------------------------------------------------------------------------------
# 5. Figure 1G: stacked-bar composition of the 12 meta-niches
# ------------------------------------------------------------------------------
# The final manuscript Figure 1G is a STACKED BAR plot, not a heatmap.
# This block follows the surviving historical Niche_revision.R logic:
#   1. average the 18-dimensional 40-um neighborhood composition within each MN;
#   2. order MNs by broad biological class, then MN number;
#   3. show the mean cell-type composition as stacked bars;
#   4. add the broad MN class as a thin annotation strip above the bars.
#
# The exact historical `celltype_colors` object was stored in the interactive R
# session and is not defined in the surviving script. We therefore use a fixed,
# deterministic 18-color palette here so repeated runs are identical. This does
# not alter any numerical composition values.
# ------------------------------------------------------------------------------

cat("\nCalculating mean neighbor composition of each MN...\n")

niche_composition <- cell_neighborhood_df %>%
  filter(!is.na(meta_niche)) %>%
  mutate(
    meta_niche = factor(
      as.character(meta_niche),
      levels = MN_LEVELS
    )
  ) %>%
  group_by(meta_niche) %>%
  summarise(
    across(
      all_of(celltype_cols),
      ~ mean(.x, na.rm = TRUE)
    ),
    n_neighborhoods = n(),
    .groups = "drop"
  ) %>%
  arrange(meta_niche)

write.csv(
  niche_composition,
  file.path(
    mn_dir,
    paste0(
      "meta_niche_composition_",
      analysis_mode,
      ".csv"
    )
  ),
  row.names = FALSE
)

# Historical Figure 1G cell-type display order.
FIG1G_CELLTYPE_ORDER <- c(
  "Melanoma Melanocytic",
  "Melanoma Proliferative",
  "Melanoma NC",
  "Melanoma Intermediate",
  "Melanoma Mesenchymal",
  "Epithelial",
  "Endothelial cells",
  "Pericytes",
  "CAFS Myofibroblast",
  "CAFS Inflammatory",
  "CD8+ T Cytotoxic",
  "CD8+ T(CXCL13)Helper",
  "Tregs",
  "B cells",
  "Plasma cells",
  "Macrophages M2",
  "Macrophages IFN-γ–activated",
  "Activated Myeloid cells"
)

# Broad MN classes recovered from the final historical Figure 1G workflow.
MN_CLASS_MAP <- c(
  "MetaNiche_1"  = "Melanoma-rich",
  "MetaNiche_2"  = "Melanoma-rich",
  "MetaNiche_3"  = "Immune-rich",
  "MetaNiche_4"  = "Vascular/perivascular",
  "MetaNiche_5"  = "Melanoma-rich",
  "MetaNiche_6"  = "Stromal-rich",
  "MetaNiche_7"  = "Stromal-rich",
  "MetaNiche_8"  = "Stromal-rich",
  "MetaNiche_9"  = "Melanoma-rich",
  "MetaNiche_10" = "Melanoma-rich",
  "MetaNiche_11" = "Immune-rich",
  "MetaNiche_12" = "Melanoma-rich"
)

MN_CLASS_ORDER <- c(
  "Melanoma-rich",
  "Immune-rich",
  "Stromal-rich",
  "Vascular/perivascular"
)

# Historical class colors.
MN_CLASS_COLORS <- c(
  "Melanoma-rich" = "#D73027",
  "Immune-rich" = "#2C7BB6",
  "Stromal-rich" = "#FEC44F",
  "Vascular/perivascular" = "#1A9850"
)

mn_order_df <- data.frame(
  meta_niche = names(MN_CLASS_MAP),
  mn_class = unname(MN_CLASS_MAP),
  stringsAsFactors = FALSE
) %>%
  mutate(
    mn_class = factor(mn_class, levels = MN_CLASS_ORDER),
    mn_number = as.integer(sub("MetaNiche_", "", meta_niche)),
    meta_niche_label = unname(MN_SHORT_MAP[meta_niche])
  ) %>%
  arrange(mn_class, mn_number)

mn_order <- mn_order_df$meta_niche
mn_label_order <- mn_order_df$meta_niche_label

# IMPORTANT FOR FIGURE 1G DISPLAY
# --------------------------------
# A valid spatial neighborhood can contain >= 3 spatial neighbors but still have
# no ANNOTATED neighbors if all neighboring cells have NA celltype. In the exact
# historical neighborhood construction, such a neighborhood is represented by
# an all-zero 18-cell-type vector. Therefore, the raw mean vector for an MN can
# sum to < 1 when an MN contains many of these zero vectors (MN10 is the clearest
# example in the reproduced data).
#
# Figure 1G is a CELL-TYPE COMPOSITION plot. To prevent annotation coverage from
# being encoded as bar height, keep the raw MN means above for provenance, then
# normalize each MN's mean 18-cell-type vector to sum to 1 FOR PLOTTING ONLY.
# This does NOT change neighborhood vectors, k-means/MN assignments, Figure 1H,
# NF1 loss-associated abundance/prevalence, or any downstream analysis.
composition_long <- niche_composition %>%
  select(
    meta_niche,
    all_of(celltype_cols)
  ) %>%
  pivot_longer(
    cols = all_of(celltype_cols),
    names_to = "celltype",
    values_to = "raw_mean_fraction"
  ) %>%
  mutate(
    meta_niche_chr = as.character(meta_niche),
    raw_mean_fraction = tidyr::replace_na(raw_mean_fraction, 0)
  ) %>%
  group_by(meta_niche_chr) %>%
  mutate(
    raw_composition_sum = sum(raw_mean_fraction),
    fraction = ifelse(
      raw_composition_sum > 0,
      raw_mean_fraction / raw_composition_sum,
      NA_real_
    )
  ) %>%
  ungroup() %>%
  mutate(
    meta_niche_label = unname(MN_SHORT_MAP[meta_niche_chr]),
    mn_class = unname(MN_CLASS_MAP[meta_niche_chr]),
    meta_niche = factor(meta_niche_chr, levels = mn_order),
    meta_niche_label = factor(meta_niche_label, levels = mn_label_order),
    mn_class = factor(mn_class, levels = MN_CLASS_ORDER),
    celltype = factor(celltype, levels = FIG1G_CELLTYPE_ORDER)
  )

# Save the exact values used for the Figure 1G stacked bars, while preserving the
# raw mean-composition CSV written above.
write.csv(
  composition_long %>%
    transmute(
      meta_niche = meta_niche_chr,
      celltype = as.character(celltype),
      raw_mean_fraction = raw_mean_fraction,
      raw_composition_sum = raw_composition_sum,
      figure1G_fraction = fraction
    ),
  file.path(
    mn_dir,
    paste0(
      "figure1G_meta_niche_composition_normalized_",
      analysis_mode,
      ".csv"
    )
  ),
  row.names = FALSE
)

fig1g_raw_sum_qc <- composition_long %>%
  distinct(meta_niche_chr, raw_composition_sum) %>%
  arrange(match(meta_niche_chr, MN_LEVELS))
cat("\nRaw MN mean-composition sums before Figure 1G normalization:\n")
print(fig1g_raw_sum_qc, row.names = FALSE)

if (any(!is.finite(composition_long$fraction))) {
  stop(
    "At least one MN has zero total annotated cell-type composition and cannot ",
    "be normalized for Figure 1G."
  )
}

# Numerical QC: after display normalization every stacked bar must sum to 1.
fig1g_sum_qc <- composition_long %>%
  group_by(meta_niche) %>%
  summarise(
    composition_sum = sum(fraction),
    .groups = "drop"
  )

if (any(abs(fig1g_sum_qc$composition_sum - 1) > 1e-8)) {
  stop("Figure 1G normalization failed: one or more MN bars do not sum to 1.")
}

# Fixed reproducible cell-type palette. The historical session-level color
# vector was not recoverable from the surviving plotting script.
CELLTYPE_COLORS_FIG1G <- setNames(
  grDevices::hcl.colors(
    length(FIG1G_CELLTYPE_ORDER),
    palette = "Dynamic"
  ),
  FIG1G_CELLTYPE_ORDER
)

mn_annot_df <- mn_order_df %>%
  mutate(
    meta_niche_label = factor(meta_niche_label, levels = mn_label_order),
    mn_class = factor(mn_class, levels = MN_CLASS_ORDER)
  )

p_top <- ggplot(
  mn_annot_df,
  aes(
    x = meta_niche_label,
    y = 1,
    fill = mn_class
  )
) +
  geom_tile(
    color = "white",
    linewidth = 0.4
  ) +
  scale_fill_manual(
    values = MN_CLASS_COLORS,
    name = "MN class"
  ) +
  scale_y_continuous(
    expand = c(0, 0)
  ) +
  labs(
    x = NULL,
    y = NULL
  ) +
  theme_void(base_size = 12) +
  theme(
    axis.text.x = element_blank(),
    legend.position = "right",
    plot.margin = margin(0, 5, 0, 5)
  )

p_stack <- ggplot(
  composition_long,
  aes(
    x = meta_niche_label,
    y = fraction,
    fill = celltype
  )
) +
  geom_col(
    width = 0.85,
    color = NA
  ) +
  scale_fill_manual(
    values = CELLTYPE_COLORS_FIG1G,
    breaks = FIG1G_CELLTYPE_ORDER,
    drop = FALSE,
    name = "Cell type"
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.02)),
    limits = c(0, 1)
  ) +
  labs(
    x = "Meta-niche",
    y = "Mean neighborhood composition"
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1,
      face = "bold"
    ),
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold"),
    legend.position = "right"
  )

p_comp <- p_top / p_stack +
  patchwork::plot_layout(
    heights = c(0.08, 1),
    guides = "collect"
  ) &
  theme(legend.position = "right")

save_plot_both(
  p_comp,
  paste0(
    "figure1G_meta_niche_composition_",
    analysis_mode
  ),
  width = 13,
  height = 7
)

# ------------------------------------------------------------------------------
# 6. Per-core MN abundance
#
# Original workflow:
#   count neighborhoods within sample/core -> divide by total valid
#   neighborhoods for that sample/core.
#
# Here we explicitly complete all 39 x 12 combinations with zero so an MN that
# is absent from a core has abundance = 0.
# ------------------------------------------------------------------------------

cat("\nCalculating per-core MN abundance...\n")

meta_niche_freq_present <- cell_neighborhood_df %>%
  filter(!is.na(meta_niche)) %>%
  count(
    sample_index,
    patient_id,
    NF1,
    meta_niche,
    name = "count"
  ) %>%
  group_by(
    sample_index,
    patient_id,
    NF1
  ) %>%
  mutate(
    freq = count / sum(count)
  ) %>%
  ungroup() %>%
  mutate(
    sample_id = as.character(patient_id)
  )

# Keep a copy of the historical "present-only" intermediate for provenance.
write.csv(
  meta_niche_freq_present,
  file.path(
    mn_dir,
    paste0(
      "meta_niche_frequency_present_only_per_core_",
      analysis_mode,
      ".csv"
    )
  ),
  row.names = FALSE
)

meta_niche_freq_complete <- sample_info %>%
  select(
    sample_index,
    patient_id,
    sample_id,
    NF1
  ) %>%
  crossing(
    meta_niche = MN_LEVELS
  ) %>%
  left_join(
    meta_niche_freq_present %>%
      select(
        sample_index,
        meta_niche,
        count,
        freq
      ),
    by = c(
      "sample_index",
      "meta_niche"
    )
  ) %>%
  mutate(
    count = replace_na(count, 0L),
    freq = replace_na(freq, 0),
    meta_niche = factor(
      meta_niche,
      levels = MN_LEVELS
    ),
    meta_niche_name = unname(
      MN_SHORT_MAP[
        as.character(meta_niche)
      ]
    )
  )

# Verify fractions sum to 1 for every core.
freq_sum_qc <- meta_niche_freq_complete %>%
  group_by(
    sample_index,
    sample_id
  ) %>%
  summarise(
    sum_freq = sum(freq),
    .groups = "drop"
  )

if (any(abs(freq_sum_qc$sum_freq - 1) > 1e-8)) {
  stop(
    "At least one core's MN fractions do not sum to 1."
  )
}

write.csv(
  meta_niche_freq_complete,
  file.path(
    mn_dir,
    paste0(
      "meta_niche_frequency_complete_per_core_",
      analysis_mode,
      ".csv"
    )
  ),
  row.names = FALSE
)

# ------------------------------------------------------------------------------
# 7. Figure 1H: number of analyzed Xenium tissue cores containing each MN
#
# Presence definition used for the archived figure:
#   any MN-assigned neighborhood in the core (freq > 0).
# ------------------------------------------------------------------------------

cat("\nCalculating Figure 1H core occurrence...\n")

mn_core_occurrence <- meta_niche_freq_complete %>%
  mutate(
    present = freq > 0
  ) %>%
  group_by(meta_niche) %>%
  summarise(
    n_cores = sum(present),
    .groups = "drop"
  ) %>%
  mutate(
    meta_niche = factor(
      as.character(meta_niche),
      levels = MN_LEVELS
    ),
    meta_niche_name = unname(
      MN_SHORT_MAP[
        as.character(meta_niche)
      ]
    )
  ) %>%
  arrange(meta_niche)

write.csv(
  mn_core_occurrence,
  file.path(
    mn_dir,
    paste0(
      "figure1H_core_occurrence_",
      analysis_mode,
      ".csv"
    )
  ),
  row.names = FALSE
)

cat("\nFigure 1H occurrence:\n")
print(
  mn_core_occurrence %>%
    select(
      meta_niche,
      n_cores
    ),
  row.names = FALSE
)

if (USE_PUBLISHED_REFERENCE) {
  observed_fig1h <- setNames(
    mn_core_occurrence$n_cores,
    as.character(mn_core_occurrence$meta_niche)
  )[MN_LEVELS]

  if (!identical(
    as.integer(observed_fig1h),
    as.integer(EXPECTED_FIG1H[MN_LEVELS])
  )) {
    stop(
      "Published-reference Figure 1H verification failed. Observed: ",
      paste(observed_fig1h, collapse = ", ")
    )
  }

  cat(
    "Figure 1H matches the archived final analysis exactly.\n"
  )
}

p_occurrence <- ggplot(
  mn_core_occurrence,
  aes(
    x = n_cores,
    y = factor(
      meta_niche_name,
      levels = rev(meta_niche_name)
    )
  )
) +
  geom_col(
    width = 0.72
  ) +
  geom_text(
    aes(label = n_cores),
    hjust = -0.25,
    size = 3.6
  ) +
  scale_x_continuous(
    limits = c(0, 42),
    breaks = seq(0, 40, by = 10),
    expand = expansion(
      mult = c(0, 0.02)
    )
  ) +
  labs(
    x = "Number of analyzed Xenium tissue cores",
    y = NULL
  ) +
  theme_classic(base_size = 12)

save_plot_both(
  p_occurrence,
  paste0(
    "figure1H_meta_niche_core_occurrence_",
    analysis_mode
  ),
  width = 7,
  height = 6
)

# ------------------------------------------------------------------------------
# 8. NF1 abundance statistics
#
# Abundance effect:
#   mean fraction of valid neighborhoods per analyzed core
#   delta = NF1Mut - NF1WT
#
# Wilcoxon test includes zero-abundance cores.
# ------------------------------------------------------------------------------

cat("\nCalculating NF1 abundance statistics...\n")

mn_abundance_stats <- meta_niche_freq_complete %>%
  group_by(meta_niche) %>%
  summarise(
    mean_mut = mean(
      freq[NF1 == "Mut"]
    ),
    mean_wt = mean(
      freq[NF1 == "WT"]
    ),
    median_mut = median(
      freq[NF1 == "Mut"]
    ),
    median_wt = median(
      freq[NF1 == "WT"]
    ),
    n_mut = sum(
      NF1 == "Mut"
    ),
    n_wt = sum(
      NF1 == "WT"
    ),
    wilcox_p = safe_wilcox(
      freq,
      NF1
    ),
    .groups = "drop"
  ) %>%
  mutate(
    delta_mut_wt = mean_mut - mean_wt,
    wilcox_p_BH = p.adjust(
      wilcox_p,
      method = "BH"
    ),
    meta_niche_name = unname(
      MN_SHORT_MAP[
        as.character(meta_niche)
      ]
    )
  ) %>%
  arrange(
    factor(
      as.character(meta_niche),
      levels = MN_LEVELS
    )
  )

write.csv(
  mn_abundance_stats,
  file.path(
    mn_dir,
    paste0(
      "meta_niche_abundance_NF1_wilcoxon_",
      analysis_mode,
      ".csv"
    )
  ),
  row.names = FALSE
)

# ------------------------------------------------------------------------------
# 9. NF1 prevalence / Fisher exact statistics
#
# Final prevalence definition:
#   present if freq > 0 (ANY presence).
#
# Fisher P value is exact.
# A 0.5 pseudocount is used ONLY for the finite plotted odds ratio.
# ------------------------------------------------------------------------------

cat("\nCalculating NF1 prevalence/Fisher statistics...\n")

mn_prevalence_stats <- bind_rows(
  lapply(
    MN_LEVELS,
    function(mn) {
      dat <- meta_niche_freq_complete %>%
        filter(
          as.character(meta_niche) == mn
        ) %>%
        mutate(
          present = freq > 0
        )

      mut_present <- sum(
        dat$NF1 == "Mut" &
          dat$present
      )
      mut_absent <- sum(
        dat$NF1 == "Mut" &
          !dat$present
      )
      wt_present <- sum(
        dat$NF1 == "WT" &
          dat$present
      )
      wt_absent <- sum(
        dat$NF1 == "WT" &
          !dat$present
      )

      contingency <- matrix(
        c(
          mut_present,
          mut_absent,
          wt_present,
          wt_absent
        ),
        nrow = 2,
        byrow = TRUE,
        dimnames = list(
          NF1 = c(
            "Mut",
            "WT"
          ),
          MN = c(
            "Present",
            "Absent"
          )
        )
      )

      fisher_obj <- fisher.test(
        contingency
      )

      odds_ratio_plot <- (
        (mut_present + 0.5) *
          (wt_absent + 0.5)
      ) / (
        (mut_absent + 0.5) *
          (wt_present + 0.5)
      )

      data.frame(
        meta_niche = mn,
        meta_niche_name = unname(
          MN_SHORT_MAP[mn]
        ),
        mut_present = mut_present,
        mut_absent = mut_absent,
        wt_present = wt_present,
        wt_absent = wt_absent,
        prevalence_mut_percent =
          100 * mut_present /
          (mut_present + mut_absent),
        prevalence_wt_percent =
          100 * wt_present /
          (wt_present + wt_absent),
        fisher_odds_ratio =
          if (length(fisher_obj$estimate) == 1L) {
            unname(fisher_obj$estimate)
          } else {
            NA_real_
          },
        odds_ratio_plot =
          odds_ratio_plot,
        log2_or_plot =
          log2(odds_ratio_plot),
        fisher_p =
          fisher_obj$p.value,
        stringsAsFactors = FALSE
      )
    }
  )
) %>%
  mutate(
    fisher_p_BH = p.adjust(
      fisher_p,
      method = "BH"
    ),
    neglog10_fisher_p = -log10(
      pmax(
        fisher_p,
        .Machine$double.xmin
      )
    )
  )

write.csv(
  mn_prevalence_stats,
  file.path(
    mn_dir,
    paste0(
      "meta_niche_prevalence_NF1_fisher_",
      analysis_mode,
      ".csv"
    )
  ),
  row.names = FALSE
)

# ------------------------------------------------------------------------------
# 10. Figure 2A-style combined statistics
#
# x / fill = mean MN abundance difference (NF1Mut - NF1WT)
# size     = -log10(Fisher exact P) for presence/absence
# ------------------------------------------------------------------------------

figure2A_stats <- mn_abundance_stats %>%
  left_join(
    mn_prevalence_stats %>%
      select(
        meta_niche,
        mut_present,
        mut_absent,
        wt_present,
        wt_absent,
        prevalence_mut_percent,
        prevalence_wt_percent,
        fisher_odds_ratio,
        odds_ratio_plot,
        log2_or_plot,
        fisher_p,
        fisher_p_BH,
        neglog10_fisher_p
      ),
    by = "meta_niche"
  ) %>%
  mutate(
    meta_niche = factor(
      as.character(meta_niche),
      levels = MN_LEVELS
    ),
    meta_niche_name = factor(
      meta_niche_name,
      levels = rev(
        unname(
          MN_SHORT_MAP[MN_LEVELS]
        )
      )
    )
  ) %>%
  arrange(meta_niche)

write.csv(
  figure2A_stats,
  file.path(
    mn_dir,
    paste0(
      "figure2A_meta_niche_NF1_stats_",
      analysis_mode,
      ".csv"
    )
  ),
  row.names = FALSE
)

p_nf1 <- ggplot(
  figure2A_stats,
  aes(
    x = delta_mut_wt,
    y = meta_niche_name,
    size = neglog10_fisher_p,
    fill = delta_mut_wt
  )
) +
  geom_point(
    shape = 21,
    color = "black",
    stroke = 0.35
  ) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed"
  ) +
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "grey85",
    high = "#B2182B",
    midpoint = 0,
    name = "Mean abundance\nNF1Mut - WT"
  ) +
  scale_size_continuous(
    range = c(2, 8),
    name = expression(
      -log[10](Fisher~italic(p))
    )
  ) +
  labs(
    x = "Mean meta-niche abundance difference: NF1Mut - WT",
    y = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.y = element_text(
      face = "bold"
    )
  )

save_plot_both(
  p_nf1,
  paste0(
    "figure2A_meta_niche_NF1_dotplot_",
    analysis_mode
  ),
  width = 9,
  height = 6
)

# ------------------------------------------------------------------------------
# 11. Compact QC summary
# ------------------------------------------------------------------------------

qc <- data.frame(
  metric = c(
    "analysis_mode",
    "valid_neighborhoods",
    "biological_celltype_columns",
    "analyzed_tissue_cores",
    "NF1Mut_cores",
    "NF1WT_cores",
    "publication_reference_changes_applied",
    "figure1H_exact_archived_match"
  ),
  value = c(
    analysis_mode,
    nrow(cell_neighborhood_df),
    length(celltype_cols),
    nrow(sample_info),
    sum(sample_info$NF1 == "Mut"),
    sum(sample_info$NF1 == "WT"),
    n_reference_changes,
    if (USE_PUBLISHED_REFERENCE) {
      TRUE
    } else {
      NA
    }
  ),
  stringsAsFactors = FALSE
)

write.csv(
  qc,
  file.path(
    mn_dir,
    paste0(
      "meta_niche_composition_NF1_QC_",
      analysis_mode,
      ".csv"
    )
  ),
  row.names = FALSE
)

cat("\nQC summary:\n")
print(
  qc,
  row.names = FALSE
)

cat("\nTop NF1 abundance differences:\n")
print(
  figure2A_stats %>%
    select(
      meta_niche,
      meta_niche_name,
      mean_mut,
      mean_wt,
      delta_mut_wt,
      fisher_p,
      fisher_p_BH
    ) %>%
    arrange(
      desc(
        abs(delta_mut_wt)
      )
    ),
  row.names = FALSE
)

cat(
  "\nDONE: MN composition / occurrence / NF1 abundance analysis complete.\n"
)
