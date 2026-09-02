# ============================================================================
# 01_time_signature_scores.R
# Figure 4A: NF1Mut vs NF1WT immune/TME signature scores
# Analysis unit: independent Xenium tissue core (39 total: 20 Mut, 19 WT)
# ============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

in_file <- file.path(
  "data", "processed",
  "merged_obj_annotated_with_meta_niche_published_reference.rds"
)

out_dir <- file.path(
  "results", "antigen_presentation", "figure4A"
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

NF1_COLORS <- c(Mut = "#D73027", WT = "#4575B4")

MELANOMA_TYPES <- c(
  "Melanoma Melanocytic",
  "Melanoma Proliferative",
  "Melanoma NC",
  "Melanoma Intermediate",
  "Melanoma Mesenchymal"
)

# Exact historical signature definitions
signatures <- list(
  ISGs = c("STAT1", "IRF1", "IRF7", "MX1", "OAS1", "GBP1"),
  Inhibitory_Checkpoints = c("PDCD1", "CD274", "CTLA4", "LAG3", "HAVCR2", "TIGIT", "VSIR"),
  Stimulatory_Checkpoints = c("CD80", "CD86", "ICOS", "CD28", "TNFRSF9"),
  Pro_Inflammatory_Cytokines = c("IL1B", "IL6", "IL12A", "IL12B", "IL17A", "IL18", "TNF", "IFNG", "IFNA1", "IFNB1"),
  Anti_Inflammatory_Cytokines = c("IL10", "IL13", "IL4", "TGFB1", "TGFB2", "TGFB3"),
  Chemokines = c("CXCL9", "CXCL10", "CXCL11", "CXCL13"),
  Chemokine_Receptors = c("CCR5", "CCR7", "CXCR3", "CXCR4", "CXCR5"),
  NK_Markers = c("KLRD1", "KLRK1"),
  T_Cell_Activation = c("IL2RA", "CD40LG"),
  Macrophage_Markers = c("CD68", "CD163", "MRC1"),
  TIS = c("CCL5", "CD27", "CD274", "CD276", "CD8A", "CMKLR1", "CXCL9", "CXCR6", "HLA-DQA1", "HLA-DRB1", "HLA-E", "IDO1", "LAG3", "NKG7", "PDCD1LG2", "PSMB10", "STAT1", "TIGIT"),
  CYT = c("GZMA", "PRF1")
)

normalize_nf1 <- function(x) {
  x <- as.character(x)
  out <- dplyr::case_when(
    x %in% c("Mut", "NF1Mut", "NF1_Mut", "NF1 Mut") ~ "Mut",
    x %in% c("WT", "NF1WT", "NF1_WT", "NF1 WT") ~ "WT",
    TRUE ~ NA_character_
  )
  if (anyNA(out)) {
    stop("Unexpected NF1 value(s): ", paste(sort(unique(x[is.na(out)])), collapse = ", "))
  }
  out
}

safe_wilcox <- function(x, g) {
  suppressWarnings(stats::wilcox.test(x ~ factor(g), exact = FALSE)$p.value)
}

if (!file.exists(in_file)) stop("Missing input: ", in_file)
cat("Loading publication-reference merged object...\n")
obj <- readRDS(in_file)

required_meta <- c("patient_id", "NF1", "celltype")
missing_meta <- setdiff(required_meta, colnames(obj@meta.data))
if (length(missing_meta) > 0L) {
  stop("Missing required metadata: ", paste(missing_meta, collapse = ", "))
}

cat("Default assay used for historical-style FetchData scoring: ", DefaultAssay(obj), "\n", sep = "")

non_melanoma_cells <- rownames(obj@meta.data)[
  !as.character(obj@meta.data$celltype) %in% MELANOMA_TYPES &
    !is.na(obj@meta.data$celltype)
]

immune_tme <- subset(obj, cells = non_melanoma_cells)
rm(obj)
gc(verbose = FALSE)

cat("Non-melanoma cells: ", format(ncol(immune_tme), big.mark = ","), "\n", sep = "")

gene_presence <- vector("list", length(signatures))
names(gene_presence) <- names(signatures)

for (sig_name in names(signatures)) {
  genes_requested <- signatures[[sig_name]]
  genes_present <- genes_requested[genes_requested %in% rownames(immune_tme)]

  gene_presence[[sig_name]] <- data.frame(
    Signature = sig_name,
    gene = genes_requested,
    present = genes_requested %in% genes_present,
    stringsAsFactors = FALSE
  )

  cat(sprintf("  %-32s %2d/%2d genes present\n", sig_name, length(genes_present), length(genes_requested)))
  if (length(genes_present) == 0L) stop("No genes present for signature: ", sig_name)

  expr <- Seurat::FetchData(immune_tme, vars = genes_present)
  immune_tme[[paste0(sig_name, "_score")]] <- rowMeans(expr, na.rm = TRUE)
}

utils::write.csv(
  do.call(rbind, gene_presence),
  file.path(out_dir, "Figure4A_signature_gene_presence.csv"),
  row.names = FALSE
)

score_cols <- paste0(names(signatures), "_score")

core_scores <- immune_tme@meta.data |>
  as.data.frame() |>
  dplyr::transmute(
    sample_id = as.character(patient_id),
    NF1 = normalize_nf1(NF1),
    dplyr::across(dplyr::all_of(score_cols))
  ) |>
  dplyr::group_by(sample_id, NF1) |>
  dplyr::summarise(
    dplyr::across(dplyr::all_of(score_cols), mean, na.rm = TRUE),
    .groups = "drop"
  )

cat("\nCore-level QC:\n")
cat("  independent tissue cores: ", dplyr::n_distinct(core_scores$sample_id), "\n", sep = "")
cat("  NF1Mut cores: ", sum(core_scores$NF1 == "Mut"), "\n", sep = "")
cat("  NF1WT cores: ", sum(core_scores$NF1 == "WT"), "\n", sep = "")

utils::write.csv(
  core_scores,
  file.path(out_dir, "Figure4A_core_signature_scores.csv"),
  row.names = FALSE
)

score_long <- tidyr::pivot_longer(
  core_scores,
  cols = dplyr::all_of(score_cols),
  names_to = "Signature",
  values_to = "Score"
)

stats_list <- lapply(score_cols, function(sig) {
  z <- score_long[score_long$Signature == sig, , drop = FALSE]
  data.frame(
    Signature = sig,
    n_mut = sum(z$NF1 == "Mut"),
    n_wt = sum(z$NF1 == "WT"),
    mean_mut = mean(z$Score[z$NF1 == "Mut"], na.rm = TRUE),
    mean_wt = mean(z$Score[z$NF1 == "WT"], na.rm = TRUE),
    delta_mut_wt = mean(z$Score[z$NF1 == "Mut"], na.rm = TRUE) - mean(z$Score[z$NF1 == "WT"], na.rm = TRUE),
    wilcox_p = safe_wilcox(z$Score, z$NF1),
    stringsAsFactors = FALSE
  )
})

fig4a_stats <- do.call(rbind, stats_list)
fig4a_stats$wilcox_p_BH <- stats::p.adjust(fig4a_stats$wilcox_p, method = "BH")
fig4a_stats <- fig4a_stats[order(fig4a_stats$wilcox_p), , drop = FALSE]

utils::write.csv(
  fig4a_stats,
  file.path(out_dir, "Figure4A_Wilcoxon_stats.csv"),
  row.names = FALSE
)

cat("\nFigure 4A statistics:\n")
print(fig4a_stats, row.names = FALSE)

DISPLAY_LABELS <- c(
  ISGs_score = "ISGs",
  Inhibitory_Checkpoints_score = "Inhibitory checkpoints",
  Stimulatory_Checkpoints_score = "Stimulatory checkpoints",
  Pro_Inflammatory_Cytokines_score = "Pro-inflammatory cytokines",
  Anti_Inflammatory_Cytokines_score = "Anti-inflammatory cytokines",
  Chemokines_score = "Chemokines",
  Chemokine_Receptors_score = "Chemokine receptors",
  NK_Markers_score = "NK markers",
  T_Cell_Activation_score = "T-cell activation",
  Macrophage_Markers_score = "Macrophage markers",
  TIS_score = "TIS",
  CYT_score = "CYT"
)

score_long$Signature_label <- unname(DISPLAY_LABELS[score_long$Signature])
score_long$Signature_label <- factor(
  score_long$Signature_label,
  levels = unname(DISPLAY_LABELS[score_cols])
)
score_long$NF1 <- factor(score_long$NF1, levels = c("Mut", "WT"))

p_fig4a <- ggplot2::ggplot(
  score_long,
  ggplot2::aes(x = NF1, y = Score, fill = NF1)
) +
  ggplot2::geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.7) +
  ggplot2::geom_jitter(width = 0.12, size = 1.5, alpha = 0.8) +
  ggplot2::facet_wrap(~ Signature_label, scales = "free_y", ncol = 4) +
  ggplot2::scale_fill_manual(values = NF1_COLORS) +
  ggplot2::labs(x = "NF1 status", y = "Mean signature score per tissue core") +
  ggplot2::theme_classic(base_size = 11) +
  ggplot2::theme(
    legend.position = "none",
    strip.background = ggplot2::element_blank(),
    strip.text = ggplot2::element_text(face = "bold")
  )

ggplot2::ggsave(
  file.path(out_dir, "Figure4A_TIME_signature_scores.pdf"),
  p_fig4a, width = 11, height = 8, useDingbats = FALSE
)

ggplot2::ggsave(
  file.path(out_dir, "Figure4A_TIME_signature_scores.png"),
  p_fig4a, width = 11, height = 8, dpi = 300
)

cat("\nDONE: Figure 4A immune/TME signature analysis complete.\n")
