# ==============================================================================
# 04_egfr_gene_correlations.R
#
# Figure 5F: Spearman correlation of EGFR expression with immune-regulatory
# genes in NF1Mut melanoma tissue cores.
#
# HISTORICAL WORKFLOW RECOVERED FROM ORIGINAL ANALYSIS
# ---------------------------------------------------
# 1. Use the publication-reference merged Seurat object.
# 2. Set the default assay to SCT, matching the historical analysis state.
# 3. Restrict to the five final melanoma states.
# 4. Restrict to NF1Mut tissue.
# 5. Fetch SCT expression for EGFR and the exact 27 genes used in the
#    historical Figure 5F plot.
# 6. Average expression within each independent tissue core (`patient_id`).
# 7. Calculate Spearman correlation of EGFR with each target gene across
#    the 20 NF1Mut cores.
# 8. Plot the per-core correlation results as the historical lollipop plot.
#
# IMPORTANT
# ---------
# - This is NOT a cell-level correlation.
# - This is NOT core x melanoma-state correlation.
# - Do NOT convert SCT values with expm1().
# - `patient_id` is used operationally as the tissue-core/specimen ID.
#
# INPUT
# -----
# data/processed/merged_obj_annotated_with_meta_niche_published_reference.rds
#
# OUTPUTS
# -------
# results/egfr/figure5F_gene_correlations/
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
})

in_file <- file.path(
  "data",
  "processed",
  "merged_obj_annotated_with_meta_niche_published_reference.rds"
)

out_dir <- file.path(
  "results",
  "egfr",
  "figure5F_gene_correlations"
)

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

MELANOMA_TYPES <- c(
  "Melanoma Mesenchymal",
  "Melanoma NC",
  "Melanoma Intermediate",
  "Melanoma Melanocytic",
  "Melanoma Proliferative"
)

# Exact gene list from the historical Figure 5F script / original plot.
GENES_OF_INTEREST <- c(
  "EGFR",
  "OAS1",
  "GBP1",
  "CD274",
  "HAVCR2",
  "TIGIT",
  "VSIR",
  "CD80",
  "IL1B",
  "IL6",
  "IL17A",
  "IFNA1",
  "IFNB1",
  "IL10",
  "TGFB1",
  "TGFB3",
  "CXCL9",
  "CXCL13",
  "CXCR3",
  "CXCR5",
  "IL2RA",
  "CD68",
  "CD163",
  "MRC1",
  "TAP1",
  "TAP2",
  "CANX",
  "CIITA"
)

# ------------------------------------------------------------------------------
# 1. Load publication-reference object
# ------------------------------------------------------------------------------

if (!file.exists(in_file)) {
  stop("Missing input: ", in_file)
}

cat("Loading publication-reference merged object...\n")
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

if (!"SCT" %in% Assays(obj)) {
  stop("SCT assay is missing from the merged object.")
}

missing_genes <- setdiff(
  GENES_OF_INTEREST,
  rownames(obj[["SCT"]])
)

if (length(missing_genes) > 0L) {
  stop(
    "Figure 5F gene(s) missing from SCT assay: ",
    paste(missing_genes, collapse = ", ")
  )
}

# Historical analysis used SCT as the active assay.
DefaultAssay(obj) <- "SCT"
Idents(obj) <- "celltype"

# ------------------------------------------------------------------------------
# 2. Subset final melanoma states and NF1Mut tissue
# ------------------------------------------------------------------------------

melanoma_cells <- subset(
  obj,
  idents = MELANOMA_TYPES
)

melanoma_mut <- subset(
  melanoma_cells,
  subset = NF1 == "Mut"
)

rm(obj, melanoma_cells)
gc(verbose = FALSE)

cat(
  "NF1Mut melanoma cells: ",
  format(ncol(melanoma_mut), big.mark = ","),
  "\n",
  "NF1Mut tissue cores: ",
  dplyr::n_distinct(as.character(melanoma_mut$patient_id)),
  "\n",
  "Default assay: ",
  DefaultAssay(melanoma_mut),
  "\n",
  sep = ""
)

if (
  dplyr::n_distinct(as.character(melanoma_mut$patient_id)) != 20L
) {
  warning(
    "Expected 20 NF1Mut tissue cores; found ",
    dplyr::n_distinct(as.character(melanoma_mut$patient_id)),
    "."
  )
}

# ------------------------------------------------------------------------------
# 3. Fetch SCT expression and average by independent tissue core
# ------------------------------------------------------------------------------

expr_data <- Seurat::FetchData(
  melanoma_mut,
  vars = c(
    GENES_OF_INTEREST,
    "patient_id"
  ),
  layer = "data"
) |>
  as.data.frame()

avg_expr <- expr_data |>
  dplyr::group_by(
    patient_id
  ) |>
  dplyr::summarise(
    dplyr::across(
      dplyr::all_of(GENES_OF_INTEREST),
      mean,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

if (nrow(avg_expr) != 20L) {
  warning(
    "Expected 20 NF1Mut core-level rows; found ",
    nrow(avg_expr),
    "."
  )
}

utils::write.csv(
  avg_expr,
  file.path(
    out_dir,
    "Figure5F_NF1Mut_melanoma_SCT_core_means.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------------------------
# 4. Spearman correlations across NF1Mut tissue cores
# ------------------------------------------------------------------------------

cor_results <- dplyr::bind_rows(
  lapply(
    GENES_OF_INTEREST[-1],
    function(gene) {

      keep <- is.finite(avg_expr$EGFR) &
        is.finite(avg_expr[[gene]])

      test <- suppressWarnings(
        stats::cor.test(
          avg_expr$EGFR[keep],
          avg_expr[[gene]][keep],
          method = "spearman",
          exact = FALSE
        )
      )

      data.frame(
        Gene = gene,
        n = sum(keep),
        Correlation = unname(test$estimate),
        P_value = test$p.value,
        stringsAsFactors = FALSE
      )
    }
  )
)

cor_results$logP <- -log10(
  pmax(
    cor_results$P_value,
    .Machine$double.xmin
  )
)

cor_results <- cor_results |>
  dplyr::arrange(
    dplyr::desc(Correlation)
  )

utils::write.csv(
  cor_results,
  file.path(
    out_dir,
    "Figure5F_EGFR_gene_Spearman_stats.csv"
  ),
  row.names = FALSE
)

cat("\nFigure 5F Spearman correlations:\n")
print(
  cor_results,
  row.names = FALSE
)

# ------------------------------------------------------------------------------
# 5. Direction QC for the relationships highlighted in the manuscript
# ------------------------------------------------------------------------------

expected_negative <- c(
  "TAP1",
  "TAP2",
  "CANX",
  "CIITA",
  "CXCL13",
  "CXCL9"
)

expected_positive <- c(
  "TGFB1",
  "TGFB3",
  "IL1B",
  "IL10",
  "CD274",
  "HAVCR2"
)

direction_qc <- data.frame(
  Gene = c(
    expected_negative,
    expected_positive
  ),
  expected_direction = c(
    rep(
      "negative",
      length(expected_negative)
    ),
    rep(
      "positive",
      length(expected_positive)
    )
  ),
  stringsAsFactors = FALSE
)

direction_qc <- dplyr::left_join(
  direction_qc,
  cor_results |>
    dplyr::select(
      Gene,
      Correlation,
      P_value
    ),
  by = "Gene"
)

direction_qc$direction_matches <- with(
  direction_qc,
  ifelse(
    expected_direction == "negative",
    Correlation < 0,
    Correlation > 0
  )
)

utils::write.csv(
  direction_qc,
  file.path(
    out_dir,
    "Figure5F_direction_QC.csv"
  ),
  row.names = FALSE
)

cat("\nDirection QC:\n")
print(
  direction_qc,
  row.names = FALSE
)

# ------------------------------------------------------------------------------
# 6. Historical lollipop plot
# ------------------------------------------------------------------------------

plot_data <- cor_results |>
  dplyr::arrange(Correlation)

p_fig5f <- ggplot2::ggplot(
  plot_data,
  ggplot2::aes(
    x = stats::reorder(Gene, Correlation),
    y = Correlation
  )
) +
  ggplot2::geom_segment(
    ggplot2::aes(
      xend = Gene,
      y = 0,
      yend = Correlation
    ),
    color = "gray80",
    linewidth = 0.5
  ) +
  ggplot2::geom_point(
    ggplot2::aes(
      color = logP
    ),
    size = 3
  ) +
  ggplot2::scale_color_gradient(
    low = "lightblue",
    high = "darkred",
    name = expression(-log[10](p))
  ) +
  ggplot2::coord_flip() +
  ggplot2::theme_minimal(
    base_size = 12
  ) +
  ggplot2::labs(
    title = "Correlation with EGFR",
    x = NULL,
    y = "Spearman Correlation"
  ) +
  ggplot2::theme(
    axis.text.y = ggplot2::element_text(
      size = 8
    ),
    legend.position = "right"
  )

ggplot2::ggsave(
  file.path(
    out_dir,
    "Figure5F_Correlation_w_EGFR.pdf"
  ),
  p_fig5f,
  width = 5.5,
  height = 5.5,
  useDingbats = FALSE
)

ggplot2::ggsave(
  file.path(
    out_dir,
    "Figure5F_Correlation_w_EGFR.png"
  ),
  p_fig5f,
  width = 5.5,
  height = 5.5,
  dpi = 300
)

cat(
  "\nDONE: Figure 5F historical SCT core-level correlation workflow complete.\n"
)
