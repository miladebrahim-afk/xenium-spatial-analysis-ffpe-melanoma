# ==============================================================================
# 03_egfr_pathway_scores.R
#
# Figure 5D: NF1Mut vs NF1WT melanoma Hallmark pathway scores.
#
# PUBLICATION-REPRODUCTION WORKFLOW
# ---------------------------------
# This script uses the EXACT historical Hallmark gene modules recovered from:
#
#   /gpfs/scratch/ibrahm07/R_Studio/Xenium/Mapping_back/Sketch/
#   EGFR_Images/EGFR_mut_WT.RData
#
# Specifically, the historical object:
#
#   pathways_list_filtered
#
# was exported to:
#
#   data/reference/figure5D_historical_pathway_modules.rds
#
# The recovered object contains 50 Hallmark modules (1,161 pathway-gene rows).
#
# Historical Figure 5D logic:
#   1. subset melanoma cells;
#   2. set RNA assay and LogNormalize;
#   3. AddModuleScore() using ALL 50 historical modules;
#   4. rename module-score columns to Hallmark pathway names;
#   5. average each pathway score within each independent tissue core;
#   6. Wilcoxon NF1Mut vs NF1WT for EACH OF THE 50 pathways;
#   7. BH-adjust across ALL 50 tests;
#   8. subset the nine pathways shown in the assembled final Figure 5D;
#   9. plot delta = mean(NF1Mut) - mean(NF1WT).
#
# IMPORTANT
# ---------
# Do NOT regenerate the Figure 5D modules from a new fgsea run. The historical
# module definitions are the publication reference. Re-running fgsea can alter
# leading-edge membership and can change pathway-score direction (P53 was the
# diagnostic example).
#
# ANALYSIS UNIT
# -------------
# Each Xenium tissue core is an independent sample. `patient_id` is used
# operationally as the tissue-core/specimen ID; paired-looking IDs are NOT merged.
#
# INPUTS
# ------
# data/processed/merged_obj_annotated_with_meta_niche_published_reference.rds
# data/reference/figure5D_historical_pathway_modules.rds
#
# OUTPUTS
# -------
# results/egfr/figure5D_pathway_scores/
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

# ------------------------------------------------------------------------------
# 0. Paths and constants
# ------------------------------------------------------------------------------

in_file <- file.path(
  "data",
  "processed",
  "merged_obj_annotated_with_meta_niche_published_reference.rds"
)

module_file <- file.path(
  "data",
  "reference",
  "figure5D_historical_pathway_modules.rds"
)

out_dir <- file.path(
  "results",
  "egfr",
  "figure5D_pathway_scores"
)

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

MELANOMA_TYPES <- c(
  "Melanoma Melanocytic",
  "Melanoma Proliferative",
  "Melanoma Intermediate",
  "Melanoma NC",
  "Melanoma Mesenchymal"
)

# Nine pathways visible in the assembled final Figure 5D.
FIG5D_PATHWAYS <- c(
  "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
  "HALLMARK_HYPOXIA",
  "HALLMARK_TGF_BETA_SIGNALING",
  "HALLMARK_IL2_STAT5_SIGNALING",
  "HALLMARK_INFLAMMATORY_RESPONSE",
  "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  "HALLMARK_INTERFERON_ALPHA_RESPONSE",
  "HALLMARK_P53_PATHWAY",
  "HALLMARK_IL6_JAK_STAT3_SIGNALING"
)

FIG5D_LABELS <- c(
  "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION" = "EMT",
  "HALLMARK_HYPOXIA" = "Hypoxia",
  "HALLMARK_TGF_BETA_SIGNALING" = "TGF B Sig.",
  "HALLMARK_IL2_STAT5_SIGNALING" = "IL2 STAT5 Sig.",
  "HALLMARK_INFLAMMATORY_RESPONSE" = "Inflamm. Resp",
  "HALLMARK_INTERFERON_GAMMA_RESPONSE" = "IFN-gamma Resp",
  "HALLMARK_INTERFERON_ALPHA_RESPONSE" = "IFN-alpha Resp",
  "HALLMARK_P53_PATHWAY" = "P53 Pathway",
  "HALLMARK_IL6_JAK_STAT3_SIGNALING" = "IL6 JAK STAT3"
)

# Historical P53 values recovered directly from the old RData workspace.
# Used only as a transparent validation target; they are NOT injected into results.
HISTORICAL_P53 <- data.frame(
  pathway = "HALLMARK_P53_PATHWAY",
  historical_mean_mut = -0.00558,
  historical_mean_wt = 0.0598,
  historical_delta = -0.0654,
  historical_p_value = 0.000485,
  historical_p_adj = 0.00347,
  stringsAsFactors = FALSE
)

# ------------------------------------------------------------------------------
# 1. Helpers
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

safe_wilcox <- function(x, g) {

  keep <- is.finite(x) & !is.na(g)

  x <- x[keep]
  g <- factor(g[keep])

  if (length(unique(g)) != 2L) {
    return(NA_real_)
  }

  suppressWarnings(
    tryCatch(
      stats::wilcox.test(
        x ~ g,
        exact = FALSE
      )$p.value,
      error = function(e) NA_real_
    )
  )
}

# ------------------------------------------------------------------------------
# 2. Load the exact historical module definitions
# ------------------------------------------------------------------------------

if (!file.exists(module_file)) {
  stop(
    "Missing historical Figure 5D module reference: ",
    module_file
  )
}

cat("Loading historical Figure 5D Hallmark modules...\n")

historical_modules <- readRDS(module_file)

if (!is.list(historical_modules)) {
  stop("Historical Figure 5D module reference is not a list.")
}

if (is.null(names(historical_modules))) {
  stop("Historical Figure 5D modules are not named.")
}

if (anyDuplicated(names(historical_modules))) {
  stop("Historical Figure 5D module names are duplicated.")
}

cat(
  "Historical Hallmark modules: ",
  length(historical_modules),
  "\n",
  sep = ""
)

if (length(historical_modules) != 50L) {
  warning(
    "Expected 50 historical Hallmark modules; found ",
    length(historical_modules),
    "."
  )
}

missing_final_pathways <- setdiff(
  FIG5D_PATHWAYS,
  names(historical_modules)
)

if (length(missing_final_pathways) > 0L) {
  stop(
    "Final Figure 5D pathway(s) absent from historical module reference: ",
    paste(missing_final_pathways, collapse = ", ")
  )
}

module_table <- dplyr::bind_rows(
  lapply(
    names(historical_modules),
    function(pw) {
      data.frame(
        pathway = pw,
        gene = historical_modules[[pw]],
        stringsAsFactors = FALSE
      )
    }
  )
)

cat(
  "Historical pathway-gene rows: ",
  nrow(module_table),
  "\n",
  sep = ""
)

utils::write.csv(
  module_table,
  file.path(
    out_dir,
    "Figure5D_historical_pathway_modules_used.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------------------------
# 3. Load merged object and subset melanoma cells
# ------------------------------------------------------------------------------

if (!file.exists(in_file)) {
  stop("Missing input: ", in_file)
}

cat("\nLoading publication-reference merged object...\n")
obj <- readRDS(in_file)

required_meta <- c(
  "patient_id",
  "NF1",
  "celltype"
)

missing_meta <- setdiff(
  required_meta,
  colnames(obj@meta.data)
)

if (length(missing_meta) > 0L) {
  stop(
    "Missing required metadata: ",
    paste(missing_meta, collapse = ", ")
  )
}

mel_cells <- rownames(obj@meta.data)[
  as.character(obj@meta.data$celltype) %in%
    MELANOMA_TYPES
]

if (length(mel_cells) == 0L) {
  stop("No melanoma cells found.")
}

mel <- subset(
  obj,
  cells = mel_cells
)

rm(obj)
gc(verbose = FALSE)

mel$NF1 <- normalize_nf1(
  mel$NF1
)

cat(
  "Melanoma cells: ",
  format(ncol(mel), big.mark = ","),
  "\n",
  "NF1Mut melanoma cells: ",
  sum(mel$NF1 == "Mut"),
  "\n",
  "NF1WT melanoma cells: ",
  sum(mel$NF1 == "WT"),
  "\n",
  sep = ""
)

# ------------------------------------------------------------------------------
# 4. Historical RNA normalization
# ------------------------------------------------------------------------------

DefaultAssay(mel) <- "RNA"

cat(
  "RNA layers before normalization: ",
  paste(
    SeuratObject::Layers(mel[["RNA"]]),
    collapse = ", "
  ),
  "\n",
  sep = ""
)

mel <- Seurat::NormalizeData(
  object = mel,
  assay = "RNA",
  normalization.method = "LogNormalize",
  verbose = FALSE
)

# Filter every historical module to genes still present in this reconstructed
# object while preserving the historical pathway order.
historical_modules_present <- lapply(
  historical_modules,
  function(genes) {
    unique(
      genes[
        genes %in% rownames(mel)
      ]
    )
  }
)

module_qc <- data.frame(
  pathway = names(historical_modules),
  historical_gene_count = vapply(
    historical_modules,
    length,
    integer(1)
  ),
  genes_present = vapply(
    historical_modules_present,
    length,
    integer(1)
  ),
  stringsAsFactors = FALSE
)

module_qc$genes_missing <-
  module_qc$historical_gene_count -
  module_qc$genes_present

utils::write.csv(
  module_qc,
  file.path(
    out_dir,
    "Figure5D_historical_module_gene_presence_QC.csv"
  ),
  row.names = FALSE
)

if (any(module_qc$genes_present == 0L)) {
  stop(
    "At least one historical pathway module has zero genes present."
  )
}

cat(
  "Historical module genes retained in reconstructed object: ",
  sum(module_qc$genes_present),
  "/",
  sum(module_qc$historical_gene_count),
  "\n",
  sep = ""
)

# ------------------------------------------------------------------------------
# 5. AddModuleScore with ALL 50 historical pathways
#
# IMPORTANT:
# - do not reorder the list;
# - do not score only the final nine;
# - BH adjustment later is across all 50 historical pathway tests.
# ------------------------------------------------------------------------------

cat("\nScoring all historical Hallmark modules...\n")

set.seed(1)

mel <- Seurat::AddModuleScore(
  object = mel,
  features = historical_modules_present,
  name = "PathwayScore",
  assay = "RNA",
  seed = 1
)

score_cols <- paste0(
  "PathwayScore",
  seq_along(historical_modules_present)
)

missing_score_cols <- setdiff(
  score_cols,
  colnames(mel@meta.data)
)

if (length(missing_score_cols) > 0L) {
  stop(
    "Missing AddModuleScore output column(s): ",
    paste(missing_score_cols, collapse = ", ")
  )
}

score_map <- data.frame(
  score_column = score_cols,
  pathway = names(historical_modules_present),
  stringsAsFactors = FALSE
)

utils::write.csv(
  score_map,
  file.path(
    out_dir,
    "Figure5D_AddModuleScore_column_map.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------------------------
# 6. Aggregate per independent tissue core
# ------------------------------------------------------------------------------

core_scores <- mel@meta.data |>
  as.data.frame() |>
  dplyr::transmute(
    sample_id = as.character(patient_id),
    NF1 = as.character(NF1),
    dplyr::across(
      dplyr::all_of(score_cols)
    )
  ) |>
  dplyr::group_by(
    sample_id,
    NF1
  ) |>
  dplyr::summarise(
    dplyr::across(
      dplyr::all_of(score_cols),
      mean,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

cat(
  "\nCore-level QC:\n",
  "  independent tissue cores: ",
  dplyr::n_distinct(core_scores$sample_id),
  "\n",
  "  NF1Mut cores: ",
  sum(core_scores$NF1 == "Mut"),
  "\n",
  "  NF1WT cores: ",
  sum(core_scores$NF1 == "WT"),
  "\n",
  sep = ""
)

if (dplyr::n_distinct(core_scores$sample_id) != 39L) {
  warning(
    "Expected 39 independent tissue cores; found ",
    dplyr::n_distinct(core_scores$sample_id),
    "."
  )
}

utils::write.csv(
  core_scores,
  file.path(
    out_dir,
    "Figure5D_core_pathway_scores_wide.csv"
  ),
  row.names = FALSE
)

long_df <- tidyr::pivot_longer(
  core_scores,
  cols = dplyr::all_of(score_cols),
  names_to = "score_column",
  values_to = "score"
)

long_df <- dplyr::left_join(
  long_df,
  score_map,
  by = "score_column"
)

utils::write.csv(
  long_df,
  file.path(
    out_dir,
    "Figure5D_core_pathway_scores_all50_long.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------------------------
# 7. Wilcoxon statistics for ALL 50 pathways
# ------------------------------------------------------------------------------

stats_list <- lapply(
  names(historical_modules_present),
  function(pw) {

    z <- long_df[
      long_df$pathway == pw,
      ,
      drop = FALSE
    ]

    mean_mut <- mean(
      z$score[z$NF1 == "Mut"],
      na.rm = TRUE
    )

    mean_wt <- mean(
      z$score[z$NF1 == "WT"],
      na.rm = TRUE
    )

    data.frame(
      pathway = pw,
      n_mut = sum(z$NF1 == "Mut"),
      n_wt = sum(z$NF1 == "WT"),
      mean_mut = mean_mut,
      mean_wt = mean_wt,
      delta = mean_mut - mean_wt,
      p_val = safe_wilcox(
        z$score,
        z$NF1
      ),
      stringsAsFactors = FALSE
    )
  }
)

all_stats <- do.call(
  rbind,
  stats_list
)

# Historical adjustment was across all pathway tests before selecting the
# display panel.
all_stats$p_val_adj <- stats::p.adjust(
  all_stats$p_val,
  method = "BH"
)

all_stats$neglog10_padj <- -log10(
  pmax(
    all_stats$p_val_adj,
    .Machine$double.xmin
  )
)

utils::write.csv(
  all_stats,
  file.path(
    out_dir,
    "Figure5D_all50_pathway_Wilcoxon_stats.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------------------------
# 8. Final nine Figure 5D pathways
# ------------------------------------------------------------------------------

fig5d_stats <- all_stats[
  match(
    FIG5D_PATHWAYS,
    all_stats$pathway
  ),
  ,
  drop = FALSE
]

fig5d_stats$label <- unname(
  FIG5D_LABELS[
    fig5d_stats$pathway
  ]
)

utils::write.csv(
  fig5d_stats,
  file.path(
    out_dir,
    "Figure5D_final9_pathway_stats.csv"
  ),
  row.names = FALSE
)

cat("\nFigure 5D final nine pathway statistics:\n")
print(
  fig5d_stats,
  row.names = FALSE
)

# ------------------------------------------------------------------------------
# 9. Historical P53 validation
# ------------------------------------------------------------------------------

p53_now <- fig5d_stats[
  fig5d_stats$pathway == "HALLMARK_P53_PATHWAY",
  c(
    "pathway",
    "mean_mut",
    "mean_wt",
    "delta",
    "p_val",
    "p_val_adj"
  ),
  drop = FALSE
]

p53_validation <- cbind(
  HISTORICAL_P53,
  current_mean_mut = p53_now$mean_mut,
  current_mean_wt = p53_now$mean_wt,
  current_delta = p53_now$delta,
  current_p_value = p53_now$p_val,
  current_p_adj = p53_now$p_val_adj
)

p53_validation$delta_difference <-
  p53_validation$current_delta -
  p53_validation$historical_delta

p53_validation$p_difference <-
  p53_validation$current_p_value -
  p53_validation$historical_p_value

utils::write.csv(
  p53_validation,
  file.path(
    out_dir,
    "Figure5D_P53_historical_validation.csv"
  ),
  row.names = FALSE
)

cat("\nHistorical P53 validation:\n")
print(
  p53_validation,
  row.names = FALSE
)

# ------------------------------------------------------------------------------
# 10. Plot
# ------------------------------------------------------------------------------

fig5d_plot <- fig5d_stats

# Preserve the assembled Figure 5D pathway order.
fig5d_plot$label <- factor(
  fig5d_plot$label,
  levels = rev(
    unname(
      FIG5D_LABELS[
        FIG5D_PATHWAYS
      ]
    )
  )
)

p_fig5d <- ggplot2::ggplot(
  fig5d_plot,
  ggplot2::aes(
    x = delta,
    y = label,
    color = neglog10_padj
  )
) +
  ggplot2::geom_point(
    size = 4
  ) +
  ggplot2::geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.6
  ) +
  ggplot2::scale_color_gradient(
    low = "grey80",
    high = "darkred",
    name = expression(
      -log[10]("adj. p-value")
    )
  ) +
  ggplot2::labs(
    title = "EGFR Pathway Enrichment",
    x = "Delta Pathway Score (Mut - WT)",
    y = NULL
  ) +
  ggplot2::theme_minimal(
    base_size = 13
  ) +
  ggplot2::theme(
    panel.grid.minor = ggplot2::element_blank(),
    plot.title = ggplot2::element_text(
      hjust = 0.5
    )
  )

ggplot2::ggsave(
  file.path(
    out_dir,
    "Figure5D_EGFR_pathway_enrichment.pdf"
  ),
  p_fig5d,
  width = 6.2,
  height = 5.3,
  useDingbats = FALSE
)

ggplot2::ggsave(
  file.path(
    out_dir,
    "Figure5D_EGFR_pathway_enrichment.png"
  ),
  p_fig5d,
  width = 6.2,
  height = 5.3,
  dpi = 300
)

cat(
  "\nDONE: Figure 5D historical-module reproduction complete.\n",
  "Positive delta = higher pathway activity in NF1Mut.\n",
  "Negative delta = higher pathway activity in NF1WT.\n",
  sep = ""
)
