# ==============================================================================
# 02_egfr_metaniche_dotplots.R
#
# Figure 5B
# EGFR prevalence/intensity across spatial meta-niches by NF1 status.
#
# FINAL CLEAN WORKFLOW
# --------------------
# - EGFR source: raw RNA counts
# - EGFR+ definition: raw EGFR count > 0
# - dot size: mean tissue-core fraction of EGFR+ cells within each meta-niche
# - dot color: mean tissue-core EGFR count among EGFR+ cells
#
# This script does NOT use SCT.
# Figure 5C / S6A are generated separately from the cached Figure 5B
# core-level metrics by:
#   02b_egfr_response_dotplots_from_cached_fig5b.R
#
# Direct Assay5 access is intentional: it extracts only the EGFR row and
# avoids materializing the full million-cell RNA matrix through LayerData().
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
})

in_file <- file.path(
  "data", "processed",
  "merged_obj_annotated_with_meta_niche_published_reference.rds"
)

out_dir <- file.path(
  "results", "egfr", "figure5BC_egfr_metaniche"
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(in_file)) {
  stop("Missing input: ", in_file)
}

cat("Loading publication-reference merged object...\n")
obj <- readRDS(in_file)
cat("Loaded: ", format(ncol(obj), big.mark = ","), " cells\n", sep = "")

required_meta <- c("patient_id", "NF1", "meta_niche")
missing_meta <- setdiff(required_meta, colnames(obj@meta.data))
if (length(missing_meta) > 0L) {
  stop("Missing metadata: ", paste(missing_meta, collapse = ", "))
}

extract_direct_layer_row <- function(assay_obj, layer, feature, cells = NULL) {
  layer_names <- names(assay_obj@layers)
  if (!layer %in% layer_names) {
    stop("Missing layer: ", layer)
  }

  feature_names <- as.character(
    SeuratObject::Features(assay_obj, layer = layer)
  )
  cell_names <- as.character(
    SeuratObject::Cells(assay_obj, layer = layer)
  )

  feature_index <- match(feature, feature_names)
  if (is.na(feature_index)) {
    stop("Feature not found in ", layer, ": ", feature)
  }

  x <- assay_obj@layers[[layer]]

  if (nrow(x) != length(feature_names) ||
      ncol(x) != length(cell_names)) {
    stop(
      "Assay layer dimensions do not agree with Seurat feature/cell maps: ",
      paste(dim(x), collapse = " x "), "."
    )
  }

  if (is.null(cells)) {
    values <- as.numeric(x[feature_index, , drop = TRUE])
    names(values) <- cell_names
  } else {
    cells <- as.character(cells)
    cell_index <- match(cells, cell_names)
    if (anyNA(cell_index)) {
      stop(sum(is.na(cell_index)), " requested cell(s) not present in layer.")
    }
    values <- as.numeric(
      x[feature_index, cell_index, drop = TRUE]
    )
    names(values) <- cells
  }

  values
}

cat("Directly extracting the EGFR row from RNA counts...\n")

rna_assay <- obj[["RNA"]]

t0 <- Sys.time()
egfr_raw <- extract_direct_layer_row(
  assay_obj = rna_assay,
  layer = "counts",
  feature = "EGFR"
)
cat(
  "EGFR row extracted for ",
  format(length(egfr_raw), big.mark = ","),
  " cells in ",
  round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 2),
  " seconds.\n",
  sep = ""
)

md <- obj@meta.data[names(egfr_raw), , drop = FALSE] |>
  as.data.frame()

md$sample_id <- as.character(md$patient_id)

nf1_raw <- as.character(md$NF1)
md$NF1_clean <- dplyr::case_when(
  nf1_raw %in% c("Mut", "NF1Mut", "NF1_Mut", "NF1 Mut") ~ "Mut",
  nf1_raw %in% c("WT", "NF1WT", "NF1_WT", "NF1 WT") ~ "WT",
  TRUE ~ NA_character_
)

if (anyNA(md$NF1_clean)) {
  stop(
    "Unexpected NF1 values: ",
    paste(
      sort(unique(nf1_raw[is.na(md$NF1_clean)])),
      collapse = ", "
    )
  )
}

md$meta_niche_clean <- as.character(md$meta_niche)
md$EGFR_raw <- egfr_raw

MN_LEVELS <- paste0("MetaNiche_", 1:12)

MN_PLOT_LEVELS <- c(
  "MetaNiche_1",
  "MetaNiche_10",
  "MetaNiche_11",
  "MetaNiche_12",
  "MetaNiche_2",
  "MetaNiche_3",
  "MetaNiche_4",
  "MetaNiche_5",
  "MetaNiche_6",
  "MetaNiche_7",
  "MetaNiche_8",
  "MetaNiche_9"
)

MN_NUMBER_LABELS <- setNames(
  as.character(c(1,10,11,12,2,3,4,5,6,7,8,9)),
  MN_PLOT_LEVELS
)

fig5b_cells <- md |>
  dplyr::filter(
    meta_niche_clean %in% MN_LEVELS,
    !is.na(sample_id),
    !is.na(NF1_clean)
  ) |>
  dplyr::mutate(
    EGFR_pos = EGFR_raw > 0
  )

fig5b_core <- fig5b_cells |>
  dplyr::group_by(
    sample_id,
    NF1_clean,
    meta_niche_clean
  ) |>
  dplyr::summarise(
    n_cells = dplyr::n(),
    n_EGFR_pos = sum(EGFR_pos),
    frac_EGFR_pos = mean(EGFR_pos),
    mean_EGFR_pos = ifelse(
      n_EGFR_pos > 0,
      mean(EGFR_raw[EGFR_pos]),
      NA_real_
    ),
    .groups = "drop"
  )

fig5b_df <- fig5b_core |>
  dplyr::group_by(
    meta_niche_clean,
    NF1_clean
  ) |>
  dplyr::summarise(
    n_cores = dplyr::n(),
    frac_EGFR_pos = mean(frac_EGFR_pos, na.rm = TRUE),
    mean_EGFR_pos = mean(mean_EGFR_pos, na.rm = TRUE),
    .groups = "drop"
  ) |>
  dplyr::mutate(
    NF1_clean = factor(NF1_clean, levels = c("Mut", "WT")),
    meta_niche_clean = factor(
      meta_niche_clean,
      levels = MN_PLOT_LEVELS
    )
  )

utils::write.csv(
  fig5b_core,
  file.path(
    out_dir,
    "Figure5B_core_level_EGFR_raw_count_metrics.csv"
  ),
  row.names = FALSE
)

utils::write.csv(
  fig5b_df,
  file.path(
    out_dir,
    "Figure5B_EGFR_dotplot_values_FINAL.csv"
  ),
  row.names = FALSE
)

cat(
  "Figure 5B fraction range: ",
  paste(range(fig5b_df$frac_EGFR_pos, na.rm = TRUE), collapse = " to "),
  "\n",
  "Figure 5B EGFR intensity range: ",
  paste(range(fig5b_df$mean_EGFR_pos, na.rm = TRUE), collapse = " to "),
  "\n",
  sep = ""
)

p_fig5b <- ggplot2::ggplot(
  fig5b_df,
  ggplot2::aes(
    x = NF1_clean,
    y = meta_niche_clean
  )
) +
  ggplot2::geom_point(
    ggplot2::aes(
      size = frac_EGFR_pos,
      color = mean_EGFR_pos
    )
  ) +
  ggplot2::scale_y_discrete(
    labels = MN_NUMBER_LABELS
  ) +
  ggplot2::scale_x_discrete(
    labels = c(
      "Mut" = expression(NF1^Mut),
      "WT" = expression(NF1^WT)
    )
  ) +
  ggplot2::scale_size_continuous(
    range = c(1, 7),
    breaks = c(0.02, 0.04, 0.06),
    name = "EGFR+"
  ) +
  ggplot2::scale_color_gradient(
    low = "grey80",
    high = "midnightblue",
    breaks = c(1.00, 1.05, 1.10, 1.15, 1.20),
    name = "EGFR"
  ) +
  ggplot2::labs(
    x = NULL,
    y = "Meta-Niche"
  ) +
  ggplot2::theme_classic(base_size = 13) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(face = "bold"),
    axis.text.y = ggplot2::element_text(size = 11),
    axis.title.y = ggplot2::element_text(size = 12),
    panel.grid.major.y = ggplot2::element_line(
      color = "grey90",
      linewidth = 0.35
    ),
    panel.grid.minor = ggplot2::element_blank(),
    legend.position = "right"
  )

ggplot2::ggsave(
  file.path(
    out_dir,
    "Figure5B_EGFR_by_NF1_meta_niche_FINAL.pdf"
  ),
  p_fig5b,
  width = 5.5,
  height = 6.3,
  useDingbats = FALSE
)

ggplot2::ggsave(
  file.path(
    out_dir,
    "Figure5B_EGFR_by_NF1_meta_niche_FINAL.png"
  ),
  p_fig5b,
  width = 5.5,
  height = 6.3,
  dpi = 300
)

cat("\nDONE — Figure 5B EGFR meta-niche analysis complete.\n")
cat("EGFR source: raw RNA counts. SCT was not accessed.\n")
