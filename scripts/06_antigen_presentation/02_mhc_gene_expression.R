# ==============================================================================
# 02_mhc_gene_expression.R
#
# Figure 4E: TAP1, TAP2, CANX, and CIITA expression in melanoma cells.
#
# FINAL PAPER ANALYSIS UNIT
# -------------------------
# Cell level.
#
# The paper's Figure 4E uses melanoma cells as the observations for the
# Wilcoxon comparison between NF1Mut and NF1WT. Tissue-core summaries are
# optionally saved as QC only and are NOT used for the Figure 4E P-values.
#
# Workflow
# --------
# 1. Load the publication-reference merged Seurat object.
# 2. Subset the five melanoma states.
# 3. Set RNA as the default assay and LogNormalize RNA counts.
# 4. Extract normalized expression for TAP1, TAP2, CANX, CIITA.
# 5. Perform two-sided cell-level Wilcoxon rank-sum tests (Mut vs WT).
# 6. Generate violin plots using the cell-level normalized expression.
#
# Input
# -----
# data/processed/merged_obj_annotated_with_meta_niche_published_reference.rds
#
# Outputs
# -------
# results/antigen_presentation/figure4E/
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})

in_file <- file.path(
  "data",
  "processed",
  "merged_obj_annotated_with_meta_niche_published_reference.rds"
)

out_dir <- file.path(
  "results",
  "antigen_presentation",
  "figure4E"
)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

MELANOMA_TYPES <- c(
  "Melanoma Melanocytic",
  "Melanoma Proliferative",
  "Melanoma Intermediate",
  "Melanoma NC",
  "Melanoma Mesenchymal"
)

GENES <- c(
  "TAP1",
  "TAP2",
  "CANX",
  "CIITA"
)

NF1_COLORS <- c(
  "Mut" = "#D73027",
  "WT"  = "#4575B4"
)

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
      paste(sort(unique(x[is.na(out)])), collapse = ", ")
    )
  }

  out
}

format_p <- function(p) {
  if (is.na(p)) return("p = NA")
  if (p < 2.2e-16) return("p < 2.2e-16")
  paste0("p = ", formatC(p, format = "g", digits = 3))
}

if (!file.exists(in_file)) {
  stop("Missing input: ", in_file)
}

cat("Loading publication-reference merged object...\n")
obj <- readRDS(in_file)

required_meta <- c("patient_id", "NF1", "celltype")
missing_meta <- setdiff(required_meta, colnames(obj@meta.data))

if (length(missing_meta) > 0L) {
  stop(
    "Missing required metadata: ",
    paste(missing_meta, collapse = ", ")
  )
}

mel_cells <- rownames(obj@meta.data)[
  as.character(obj@meta.data$celltype) %in% MELANOMA_TYPES
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

cat(
  "Melanoma cells selected: ",
  format(ncol(mel), big.mark = ","),
  "\n",
  sep = ""
)

DefaultAssay(mel) <- "RNA"

missing_genes <- setdiff(GENES, rownames(mel))
if (length(missing_genes) > 0L) {
  stop(
    "Missing Figure 4E gene(s): ",
    paste(missing_genes, collapse = ", ")
  )
}

cat(
  "RNA layers before normalization: ",
  paste(SeuratObject::Layers(mel[["RNA"]]), collapse = ", "),
  "\n",
  sep = ""
)

mel <- Seurat::NormalizeData(
  mel,
  assay = "RNA",
  normalization.method = "LogNormalize",
  scale.factor = 10000,
  verbose = FALSE
)

nf1 <- normalize_nf1(mel$NF1)

cat(
  "NF1Mut melanoma cells: ",
  sum(nf1 == "Mut"),
  "\n",
  "NF1WT melanoma cells: ",
  sum(nf1 == "WT"),
  "\n",
  sep = ""
)

# ------------------------------------------------------------------------------
# Cell-level Figure 4E statistics
# ------------------------------------------------------------------------------

stats_list <- vector("list", length(GENES))
names(stats_list) <- GENES

plot_list <- vector("list", length(GENES))
names(plot_list) <- GENES

for (g in GENES) {

  expr <- Seurat::FetchData(
    mel,
    vars = g,
    layer = "data"
  )[[g]]

  z <- data.frame(
    expression = expr,
    NF1 = factor(
      nf1,
      levels = c("Mut", "WT")
    ),
    stringsAsFactors = FALSE
  )

  wt <- stats::wilcox.test(
    expression ~ NF1,
    data = z,
    exact = FALSE
  )

  mean_mut <- mean(
    z$expression[z$NF1 == "Mut"],
    na.rm = TRUE
  )

  mean_wt <- mean(
    z$expression[z$NF1 == "WT"],
    na.rm = TRUE
  )

  median_mut <- stats::median(
    z$expression[z$NF1 == "Mut"],
    na.rm = TRUE
  )

  median_wt <- stats::median(
    z$expression[z$NF1 == "WT"],
    na.rm = TRUE
  )

  stats_list[[g]] <- data.frame(
    gene = g,
    analysis_unit = "cell",
    n_mut_cells = sum(z$NF1 == "Mut"),
    n_wt_cells = sum(z$NF1 == "WT"),
    mean_mut = mean_mut,
    mean_wt = mean_wt,
    median_mut = median_mut,
    median_wt = median_wt,
    delta_mut_wt = mean_mut - mean_wt,
    wilcox_W = unname(wt$statistic),
    wilcox_p = wt$p.value,
    stringsAsFactors = FALSE
  )

  ymax <- max(z$expression, na.rm = TRUE)

  plot_list[[g]] <- ggplot2::ggplot(
    z,
    ggplot2::aes(
      x = NF1,
      y = expression,
      fill = NF1
    )
  ) +
    ggplot2::geom_violin(
      scale = "width",
      trim = TRUE,
      linewidth = 0.25
    ) +
    ggplot2::stat_summary(
      fun = stats::median,
      geom = "point",
      shape = 21,
      size = 2.2,
      fill = "white"
    ) +
    ggplot2::annotate(
      "text",
      x = 1.5,
      y = ymax * 1.03,
      label = format_p(wt$p.value),
      size = 3.3
    ) +
    ggplot2::scale_fill_manual(
      values = NF1_COLORS
    ) +
    ggplot2::labs(
      title = g,
      x = NULL,
      y = "Log-normalized RNA expression"
    ) +
    ggplot2::theme_classic(base_size = 11) +
    ggplot2::theme(
      legend.position = "none",
      plot.title = ggplot2::element_text(
        face = "bold",
        hjust = 0.5
      ),
      axis.text.x = ggplot2::element_text(
        face = "bold"
      )
    )
}

fig4e_stats <- do.call(
  rbind,
  stats_list
)

fig4e_stats$wilcox_p_BH <- stats::p.adjust(
  fig4e_stats$wilcox_p,
  method = "BH"
)

utils::write.csv(
  fig4e_stats,
  file.path(
    out_dir,
    "Figure4E_cell_level_Wilcoxon_stats.csv"
  ),
  row.names = FALSE
)

cat("\nFigure 4E cell-level statistics:\n")
print(fig4e_stats, row.names = FALSE)

# ------------------------------------------------------------------------------
# Optional core-level sensitivity/QC
#
# These values are NOT used for the paper Figure 4E P-values.
# ------------------------------------------------------------------------------

expr_all <- Seurat::FetchData(
  mel,
  vars = c(
    "patient_id",
    GENES
  ),
  layer = "data"
)

expr_all$NF1 <- nf1

core_means <- expr_all |>
  dplyr::group_by(
    patient_id,
    NF1
  ) |>
  dplyr::summarise(
    dplyr::across(
      dplyr::all_of(GENES),
      mean,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

core_stats <- do.call(
  rbind,
  lapply(
    GENES,
    function(g) {

      x <- core_means[[g]]
      grp <- core_means$NF1

      wt <- suppressWarnings(
        stats::wilcox.test(
          x ~ factor(grp),
          exact = FALSE
        )
      )

      data.frame(
        gene = g,
        analysis_unit = "tissue_core_QC_only",
        n_mut = sum(grp == "Mut"),
        n_wt = sum(grp == "WT"),
        mean_mut = mean(x[grp == "Mut"], na.rm = TRUE),
        mean_wt = mean(x[grp == "WT"], na.rm = TRUE),
        delta_mut_wt =
          mean(x[grp == "Mut"], na.rm = TRUE) -
          mean(x[grp == "WT"], na.rm = TRUE),
        wilcox_p = wt$p.value,
        stringsAsFactors = FALSE
      )
    }
  )
)

utils::write.csv(
  core_stats,
  file.path(
    out_dir,
    "Figure4E_core_level_sensitivity_QC_NOT_FIGURE_STATS.csv"
  ),
  row.names = FALSE
)

utils::write.csv(
  core_means,
  file.path(
    out_dir,
    "Figure4E_core_mean_expression_QC.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------------------------
# Figure
# ------------------------------------------------------------------------------

p_fig4e <- patchwork::wrap_plots(
  plot_list,
  nrow = 1
)

ggplot2::ggsave(
  file.path(
    out_dir,
    "Figure4E_MHC_gene_expression_violin.pdf"
  ),
  p_fig4e,
  width = 10,
  height = 4.2,
  useDingbats = FALSE
)

ggplot2::ggsave(
  file.path(
    out_dir,
    "Figure4E_MHC_gene_expression_violin.png"
  ),
  p_fig4e,
  width = 10,
  height = 4.2,
  dpi = 300
)

cat(
  "\nDONE: Figure 4E cell-level melanoma analysis complete.\n"
)
