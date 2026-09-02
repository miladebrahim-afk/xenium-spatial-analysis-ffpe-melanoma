# ==============================================================================
# 04_cd8_core_pseudobulk_de.R
#
# NF1 melanoma Xenium analysis
# Figure 3F-G: CD8 T-cell core-level pseudobulk differential expression
# and GO Biological Process enrichment.
#
# FINAL ANALYSIS UNIT
# -------------------
# Each analyzed Xenium tissue core is treated as an independent sample.
# Core IDs are NEVER merged (e.g. "16-275" and "16-275_2" remain separate).
#
# CD8 populations
# ---------------
#   - CD8+ T Cytotoxic
#   - CD8+ T(CXCL13)Helper
#
# Figure 3F
# ---------
# Raw RNA counts are summed within each tissue core, DESeq2 is run with
# design = ~ NF1, and the six manuscript genes are displayed using
# log2(DESeq2-normalized pseudobulk counts + 1).
# Displayed P-values are two-sided Wilcoxon tests across tissue cores.
#
# Figure 3G
# ---------
# Significant DE genes (DESeq2 padj < 0.05) are split into NF1Mut-up and
# NF1Mut-down sets. GO Biological Process enrichment is run separately
# with clusterProfiler::enrichGO(), matching the surviving historical CD8
# workflow. The full significant GO tables are saved; the plotted panel is
# a curated, non-redundant subset of biologically representative terms.
#
# Input
# -----
# data/processed/merged_obj_annotated_with_meta_niche_published_reference.rds
#
# Outputs
# -------
# results/spatial_cd8/cd8_pseudobulk/
#
# NOTE
# ----
# This script deliberately namespaces functions that are commonly masked by
# Bioconductor packages (e.g. dplyr::select, dplyr::filter, DESeq2::counts).
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(DESeq2)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(clusterProfiler)
  library(org.Hs.eg.db)
})

# ------------------------------------------------------------------------------
# 0. Paths and constants
# ------------------------------------------------------------------------------

merged_file <- file.path(
  "data",
  "processed",
  "merged_obj_annotated_with_meta_niche_published_reference.rds"
)

out_dir <- file.path(
  "results",
  "spatial_cd8",
  "cd8_pseudobulk"
)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

CD8_TYPES <- c(
  "CD8+ T Cytotoxic",
  "CD8+ T(CXCL13)Helper"
)

FIG3F_GENES <- c(
  "FYN",
  "GZMB",
  "IRF1",
  "TNFRSF1B",
  "TNFSF4",
  "TNFRSF9"
)

DE_PADJ_CUTOFF <- 0.05

NF1_COLORS <- c(
  "Mut" = "#D73027",
  "WT" = "#4575B4"
)

# Curated Figure 3G terms.
# Full significant GO tables are saved separately, so this only controls
# which representative terms appear in the compact manuscript-style panel.
FIG3G_DOWN_IDS <- c(
  "GO:0046651",  # lymphocyte proliferation
  "GO:0030098",  # lymphocyte differentiation
  "GO:0001819",  # positive regulation of cytokine production
  "GO:0002286",  # T cell activation involved in immune response
  "GO:0030217",  # T cell differentiation
  "GO:0001906",  # cell killing
  "GO:0042267",  # natural killer cell mediated cytotoxicity
  "GO:0002697",  # regulation of immune effector process
  "GO:0002443"   # leukocyte mediated immunity
)

FIG3G_UP_IDS <- c(
  "GO:0036293"   # response to decreased oxygen levels
)

FIG3G_LABELS <- c(
  "GO:0046651" = "Lymphocyte proliferation",
  "GO:0030098" = "Lymphocyte differentiation",
  "GO:0001819" = "Cytokine production",
  "GO:0002286" = "T-cell activation",
  "GO:0030217" = "T-cell differentiation",
  "GO:0001906" = "Cell killing",
  "GO:0042267" = "NK-cell cytotoxicity",
  "GO:0002697" = "Immune effector regulation",
  "GO:0002443" = "Leukocyte-mediated immunity",
  "GO:0036293" = "Response to decreased oxygen"
)

# ------------------------------------------------------------------------------
# 1. Helper functions
# ------------------------------------------------------------------------------

normalize_nf1 <- function(x) {
  x <- as.character(x)

  out <- ifelse(
    x %in% c("Mut", "NF1Mut", "NF1_Mut", "NF1 Mut"),
    "Mut",
    ifelse(
      x %in% c("WT", "NF1WT", "NF1_WT", "NF1 WT"),
      "WT",
      NA_character_
    )
  )

  if (anyNA(out)) {
    stop(
      "Unexpected NF1 value(s): ",
      paste(sort(unique(x[is.na(out)])), collapse = ", ")
    )
  }

  out
}

safe_wilcox <- function(values, groups) {
  suppressWarnings(
    tryCatch(
      stats::wilcox.test(
        values ~ factor(groups),
        exact = FALSE
      )$p.value,
      error = function(e) NA_real_
    )
  )
}

ratio_to_num <- function(x) {
  z <- strsplit(as.character(x), "/", fixed = TRUE)

  vapply(
    z,
    function(y) {
      if (length(y) != 2L) return(NA_real_)
      as.numeric(y[1]) / as.numeric(y[2])
    },
    numeric(1)
  )
}

convert_to_entrez <- function(genes) {
  genes <- unique(as.character(genes))

  if (length(genes) == 0L) {
    return(character(0))
  }

  map <- suppressMessages(
    clusterProfiler::bitr(
      genes,
      fromType = "SYMBOL",
      toType = "ENTREZID",
      OrgDb = org.Hs.eg.db
    )
  )

  unique(map$ENTREZID)
}

# ------------------------------------------------------------------------------
# 2. Load publication-reference merged object
# ------------------------------------------------------------------------------

if (!file.exists(merged_file)) {
  stop("Missing input: ", merged_file)
}

cat("Loading publication-reference merged object...\n")
obj <- readRDS(merged_file)

required_meta <- c("patient_id", "NF1", "celltype")
missing_meta <- setdiff(required_meta, colnames(obj@meta.data))

if (length(missing_meta) > 0L) {
  stop(
    "Missing required metadata: ",
    paste(missing_meta, collapse = ", ")
  )
}

# Select the two CD8 populations.
cd8_cell_names <- rownames(obj@meta.data)[
  as.character(obj@meta.data$celltype) %in% CD8_TYPES
]

if (length(cd8_cell_names) == 0L) {
  stop("No CD8 cells found.")
}

cat(
  "CD8 cells selected: ",
  format(length(cd8_cell_names), big.mark = ","),
  "\n",
  sep = ""
)

cd8 <- subset(obj, cells = cd8_cell_names)

rm(obj)
gc(verbose = FALSE)

DefaultAssay(cd8) <- "RNA"

# ------------------------------------------------------------------------------
# 3. Extract raw RNA counts
# ------------------------------------------------------------------------------

rna_layers <- SeuratObject::Layers(cd8[["RNA"]])

cat(
  "RNA layers: ",
  paste(rna_layers, collapse = ", "),
  "\n",
  sep = ""
)

if ("counts" %in% rna_layers) {

  raw_counts <- SeuratObject::LayerData(
    cd8,
    assay = "RNA",
    layer = "counts"
  )

} else {

  count_layers <- grep("^counts", rna_layers, value = TRUE)

  if (length(count_layers) == 0L) {
    stop("No RNA count layer found.")
  }

  cat(
    "Joining ",
    length(count_layers),
    " RNA count layers...\n",
    sep = ""
  )

  cd8[["RNA"]] <- SeuratObject::JoinLayers(
    cd8[["RNA"]],
    layers = count_layers,
    new = "counts"
  )

  raw_counts <- SeuratObject::LayerData(
    cd8,
    assay = "RNA",
    layer = "counts"
  )
}

meta <- data.frame(
  core_id = as.character(cd8@meta.data$patient_id),
  NF1 = normalize_nf1(cd8@meta.data$NF1),
  row.names = rownames(cd8@meta.data),
  stringsAsFactors = FALSE
)

if (!identical(colnames(raw_counts), rownames(meta))) {
  stop("Raw-count columns do not align with CD8 metadata rows.")
}

# Verify one genotype per independent core.
core_nf1_check <- unique(meta[, c("core_id", "NF1")])

nf1_per_core <- table(core_nf1_check$core_id)

if (any(nf1_per_core != 1L)) {
  stop("At least one core has multiple NF1 labels.")
}

# ------------------------------------------------------------------------------
# 4. Build independent core-level pseudobulks
#
# IMPORTANT: do not alter core IDs and do not combine "_2" cores.
# ------------------------------------------------------------------------------

cat("\nBuilding independent tissue-core pseudobulks...\n")

cells_by_core <- split(
  rownames(meta),
  meta$core_id
)

pseudobulk_list <- lapply(
  cells_by_core,
  function(cells) {
    Matrix::rowSums(
      raw_counts[, cells, drop = FALSE]
    )
  }
)

pseudobulk_counts <- do.call(cbind, pseudobulk_list)
colnames(pseudobulk_counts) <- names(pseudobulk_list)

core_meta <- unique(meta[, c("core_id", "NF1")])
rownames(core_meta) <- core_meta$core_id

core_meta <- core_meta[
  colnames(pseudobulk_counts),
  "NF1",
  drop = FALSE
]

core_meta$NF1 <- factor(
  core_meta$NF1,
  levels = c("WT", "Mut")
)

if (!identical(colnames(pseudobulk_counts), rownames(core_meta))) {
  stop("Pseudobulk columns do not align with core metadata.")
}

core_cd8_qc <- data.frame(
  core_id = names(cells_by_core),
  n_cd8_cells = lengths(cells_by_core),
  NF1 = as.character(
    core_meta[names(cells_by_core), "NF1"]
  ),
  stringsAsFactors = FALSE
)

utils::write.csv(
  core_cd8_qc,
  file.path(out_dir, "CD8_core_pseudobulk_sample_QC.csv"),
  row.names = FALSE
)

cat("\nCore-level pseudobulk QC:\n")
cat(
  "  independent tissue cores: ",
  ncol(pseudobulk_counts),
  "\n",
  sep = ""
)
cat(
  "  NF1Mut cores: ",
  sum(core_meta$NF1 == "Mut"),
  "\n",
  sep = ""
)
cat(
  "  NF1WT cores: ",
  sum(core_meta$NF1 == "WT"),
  "\n",
  sep = ""
)

# Remove genes with zero pseudobulk counts across every core.
keep_gene <- Matrix::rowSums(pseudobulk_counts) > 0

pseudobulk_counts <- pseudobulk_counts[
  keep_gene,
  ,
  drop = FALSE
]

cat(
  "  genes with non-zero pseudobulk counts: ",
  nrow(pseudobulk_counts),
  "\n",
  sep = ""
)

# ------------------------------------------------------------------------------
# 5. DESeq2: NF1Mut versus NF1WT
# ------------------------------------------------------------------------------

dds <- DESeq2::DESeqDataSetFromMatrix(
  countData = round(pseudobulk_counts),
  colData = core_meta,
  design = ~ NF1
)

dds <- DESeq2::DESeq(dds)

res <- DESeq2::results(
  dds,
  contrast = c("NF1", "Mut", "WT")
)

res_df <- as.data.frame(res)
res_df$gene <- rownames(res_df)
rownames(res_df) <- NULL

res_df$group_higher <- ifelse(
  res_df$log2FoldChange > 0,
  "Mut",
  ifelse(
    res_df$log2FoldChange < 0,
    "WT",
    "Equal"
  )
)

res_df <- res_df[
  order(res_df$pvalue),
  ,
  drop = FALSE
]

utils::write.csv(
  res_df,
  file.path(
    out_dir,
    "CD8_NF1Mut_vs_WT_DESeq2_core_pseudobulk.csv"
  ),
  row.names = FALSE
)

saveRDS(
  dds,
  file.path(
    out_dir,
    "CD8_NF1Mut_vs_WT_DESeq2_core_pseudobulk_dds.rds"
  ),
  compress = FALSE
)

# ------------------------------------------------------------------------------
# 6. Figure 3F: selected CD8 genes
# ------------------------------------------------------------------------------

normalized_counts <- DESeq2::counts(
  dds,
  normalized = TRUE
)

genes_present <- intersect(
  FIG3F_GENES,
  rownames(normalized_counts)
)

genes_missing <- setdiff(
  FIG3F_GENES,
  genes_present
)

if (length(genes_missing) > 0L) {
  warning(
    "Figure 3F genes missing from RNA counts: ",
    paste(genes_missing, collapse = ", ")
  )
}

fig3f_wide <- as.data.frame(
  t(normalized_counts[genes_present, , drop = FALSE])
)

fig3f_wide$core_id <- rownames(fig3f_wide)
rownames(fig3f_wide) <- NULL

fig3f_wide$NF1 <- as.character(
  core_meta[fig3f_wide$core_id, "NF1"]
)

fig3f_long <- tidyr::pivot_longer(
  fig3f_wide,
  cols = dplyr::all_of(genes_present),
  names_to = "gene",
  values_to = "normalized_count"
)

fig3f_long$log2_normalized_expression <- log2(
  fig3f_long$normalized_count + 1
)

fig3f_long$gene <- factor(
  fig3f_long$gene,
  levels = FIG3F_GENES
)

fig3f_long$NF1 <- factor(
  fig3f_long$NF1,
  levels = c("Mut", "WT")
)

fig3f_stats_list <- lapply(
  FIG3F_GENES,
  function(g) {

    x <- fig3f_long[
      as.character(fig3f_long$gene) == g,
      ,
      drop = FALSE
    ]

    data.frame(
      gene = g,
      n_mut = sum(x$NF1 == "Mut"),
      n_wt = sum(x$NF1 == "WT"),
      mean_mut = mean(
        x$log2_normalized_expression[x$NF1 == "Mut"],
        na.rm = TRUE
      ),
      mean_wt = mean(
        x$log2_normalized_expression[x$NF1 == "WT"],
        na.rm = TRUE
      ),
      median_mut = stats::median(
        x$log2_normalized_expression[x$NF1 == "Mut"],
        na.rm = TRUE
      ),
      median_wt = stats::median(
        x$log2_normalized_expression[x$NF1 == "WT"],
        na.rm = TRUE
      ),
      wilcox_p = safe_wilcox(
        x$log2_normalized_expression,
        x$NF1
      ),
      stringsAsFactors = FALSE
    )
  }
)

fig3f_stats <- do.call(
  rbind,
  fig3f_stats_list
)

fig3f_stats$Mut_minus_WT_mean <-
  fig3f_stats$mean_mut -
  fig3f_stats$mean_wt

utils::write.csv(
  fig3f_stats,
  file.path(
    out_dir,
    "Figure3F_selected_gene_Wilcoxon_stats.csv"
  ),
  row.names = FALSE
)

fig3f_deseq <- res_df[
  res_df$gene %in% FIG3F_GENES,
  ,
  drop = FALSE
]

fig3f_deseq <- fig3f_deseq[
  match(FIG3F_GENES, fig3f_deseq$gene),
  ,
  drop = FALSE
]

utils::write.csv(
  fig3f_deseq,
  file.path(
    out_dir,
    "Figure3F_selected_gene_DESeq2_stats.csv"
  ),
  row.names = FALSE
)

panel_y <- dplyr::summarise(
  dplyr::group_by(fig3f_long, gene),
  ymax = max(log2_normalized_expression, na.rm = TRUE),
  ymin = min(log2_normalized_expression, na.rm = TRUE),
  .groups = "drop"
)

panel_y <- dplyr::left_join(
  panel_y,
  fig3f_stats[, c("gene", "wilcox_p")],
  by = "gene"
)

panel_y$ypos <- panel_y$ymax +
  0.12 * pmax(panel_y$ymax - panel_y$ymin, 0.25)

panel_y$p_label <- ifelse(
  panel_y$wilcox_p < 0.001,
  "p < 0.001",
  paste0(
    "p = ",
    formatC(
      panel_y$wilcox_p,
      format = "f",
      digits = 3
    )
  )
)

p_fig3f <- ggplot2::ggplot(
  fig3f_long,
  ggplot2::aes(
    x = NF1,
    y = log2_normalized_expression,
    fill = NF1
  )
) +
  ggplot2::geom_boxplot(
    width = 0.55,
    outlier.shape = NA,
    alpha = 0.65
  ) +
  ggplot2::geom_jitter(
    width = 0.12,
    size = 1.7,
    alpha = 0.75
  ) +
  ggplot2::facet_wrap(
    ~ gene,
    scales = "free_y",
    nrow = 1
  ) +
  ggplot2::geom_text(
    data = panel_y,
    ggplot2::aes(
      x = 1.5,
      y = ypos,
      label = p_label
    ),
    inherit.aes = FALSE,
    size = 3
  ) +
  ggplot2::scale_fill_manual(
    values = NF1_COLORS
  ) +
  ggplot2::labs(
    x = "NF1 status",
    y = expression(
      log[2]("DESeq2-normalized pseudobulk count + 1")
    )
  ) +
  ggplot2::theme_classic(base_size = 11) +
  ggplot2::theme(
    legend.position = "none",
    strip.background = ggplot2::element_blank(),
    strip.text = ggplot2::element_text(face = "bold"),
    axis.text.x = ggplot2::element_text(face = "bold")
  )

ggplot2::ggsave(
  file.path(
    out_dir,
    "Figure3F_CD8_selected_genes_core_pseudobulk.pdf"
  ),
  p_fig3f,
  width = 15,
  height = 4.2,
  useDingbats = FALSE
)

ggplot2::ggsave(
  file.path(
    out_dir,
    "Figure3F_CD8_selected_genes_core_pseudobulk.png"
  ),
  p_fig3f,
  width = 15,
  height = 4.2,
  dpi = 300
)

cat("\nFigure 3F selected-gene Wilcoxon statistics:\n")
print(fig3f_stats, row.names = FALSE)

# ------------------------------------------------------------------------------
# 7. Significant DE genes for GO analysis
# ------------------------------------------------------------------------------

sig <- res_df[
  !is.na(res_df$padj) &
    res_df$padj < DE_PADJ_CUTOFF &
    !is.na(res_df$log2FoldChange),
  ,
  drop = FALSE
]

up_genes <- sig$gene[
  sig$log2FoldChange > 0
]

down_genes <- sig$gene[
  sig$log2FoldChange < 0
]

cat(
  "\nSignificant DE genes: ",
  nrow(sig),
  "\n",
  "  NF1Mut up: ",
  length(up_genes),
  "\n",
  "  NF1Mut down: ",
  length(down_genes),
  "\n",
  sep = ""
)

# ------------------------------------------------------------------------------
# 8. GO Biological Process enrichment
#
# Historical CD8 logic:
#   - split significant genes into up/down
#   - map each set to Entrez IDs
#   - enrichGO separately
# ------------------------------------------------------------------------------

up_entrez <- convert_to_entrez(up_genes)
down_entrez <- convert_to_entrez(down_genes)

cat(
  "Mapped NF1Mut-up genes: ",
  length(up_entrez),
  "\n",
  "Mapped NF1Mut-down genes: ",
  length(down_entrez),
  "\n",
  sep = ""
)

go_up <- clusterProfiler::enrichGO(
  gene = up_entrez,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2,
  readable = TRUE
)

go_down <- clusterProfiler::enrichGO(
  gene = down_entrez,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2,
  readable = TRUE
)

go_up_df <- as.data.frame(go_up)
go_down_df <- as.data.frame(go_down)

go_up_sig <- go_up_df[
  !is.na(go_up_df$p.adjust) &
    go_up_df$p.adjust < 0.05,
  ,
  drop = FALSE
]

go_down_sig <- go_down_df[
  !is.na(go_down_df$p.adjust) &
    go_down_df$p.adjust < 0.05,
  ,
  drop = FALSE
]

cat(
  "\nSignificant GO BP terms:\n",
  "  NF1Mut up: ",
  nrow(go_up_sig),
  "\n",
  "  NF1Mut down: ",
  nrow(go_down_sig),
  "\n",
  sep = ""
)

# Save complete significant enrichment results before any plotting selection.
utils::write.csv(
  go_up_sig,
  file.path(
    out_dir,
    "Figure3G_GO_BP_all_significant_NF1Mut_up.csv"
  ),
  row.names = FALSE
)

utils::write.csv(
  go_down_sig,
  file.path(
    out_dir,
    "Figure3G_GO_BP_all_significant_NF1Mut_down.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------------------------
# 9. Figure 3G: curated representative GO panel
#
# Full descriptions remain in the saved table. Short labels are used only
# for the visual panel.
# ------------------------------------------------------------------------------

fig3g_down <- go_down_sig[
  go_down_sig$ID %in% FIG3G_DOWN_IDS,
  c(
    "ID",
    "Description",
    "GeneRatio",
    "BgRatio",
    "p.adjust",
    "Count"
  ),
  drop = FALSE
]

fig3g_up <- go_up_sig[
  go_up_sig$ID %in% FIG3G_UP_IDS,
  c(
    "ID",
    "Description",
    "GeneRatio",
    "BgRatio",
    "p.adjust",
    "Count"
  ),
  drop = FALSE
]

missing_down <- setdiff(
  FIG3G_DOWN_IDS,
  fig3g_down$ID
)

missing_up <- setdiff(
  FIG3G_UP_IDS,
  fig3g_up$ID
)

if (length(missing_down) > 0L) {
  warning(
    "Selected NF1Mut-down GO terms not significant/present: ",
    paste(missing_down, collapse = ", ")
  )
}

if (length(missing_up) > 0L) {
  warning(
    "Selected NF1Mut-up GO terms not significant/present: ",
    paste(missing_up, collapse = ", ")
  )
}

fig3g_down$Direction <- "Downregulated"
fig3g_up$Direction <- "Upregulated"

fig3g_down$GeneRatioNum <- ratio_to_num(fig3g_down$GeneRatio)
fig3g_down$BgRatioNum <- ratio_to_num(fig3g_down$BgRatio)
fig3g_down$FoldEnrichment <-
  fig3g_down$GeneRatioNum /
  fig3g_down$BgRatioNum
fig3g_down$SignedFE <- -fig3g_down$FoldEnrichment

fig3g_up$GeneRatioNum <- ratio_to_num(fig3g_up$GeneRatio)
fig3g_up$BgRatioNum <- ratio_to_num(fig3g_up$BgRatio)
fig3g_up$FoldEnrichment <-
  fig3g_up$GeneRatioNum /
  fig3g_up$BgRatioNum
fig3g_up$SignedFE <- fig3g_up$FoldEnrichment

fig3g <- rbind(
  fig3g_down,
  fig3g_up
)

fig3g$log10p <- -log10(fig3g$p.adjust)
fig3g$Label <- unname(FIG3G_LABELS[fig3g$ID])

# Keep the full original GO term in the output table and use only the
# abbreviated Label for plotting.
utils::write.csv(
  fig3g,
  file.path(
    out_dir,
    "Figure3G_selected_GO_BP_terms.csv"
  ),
  row.names = FALSE
)

# Order from strongest negative enrichment at bottom to positive oxygen response
# at top.
fig3g <- fig3g[
  order(fig3g$SignedFE),
  ,
  drop = FALSE
]

fig3g$Label <- factor(
  fig3g$Label,
  levels = fig3g$Label
)

p_fig3g <- ggplot2::ggplot(
  fig3g,
  ggplot2::aes(
    x = SignedFE,
    y = Label,
    color = log10p,
    size = Count
  )
) +
  ggplot2::geom_point(alpha = 0.95) +
  ggplot2::geom_vline(
    xintercept = 0,
    linetype = "dashed",
    color = "grey50"
  ) +
  ggplot2::scale_color_gradient(
    low = "#4575B4",
    high = "#D73027",
    name = expression(-log[10]("adjusted p"))
  ) +
  ggplot2::scale_size_continuous(
    name = "Gene count",
    range = c(2.5, 8)
  ) +
  ggplot2::scale_x_continuous(
    expand = ggplot2::expansion(
      mult = c(0.05, 0.08)
    )
  ) +
  ggplot2::labs(
    x = "Signed fold enrichment",
    y = NULL
  ) +
  ggplot2::theme_classic(base_size = 12) +
  ggplot2::theme(
    axis.text.y = ggplot2::element_text(size = 11),
    legend.title = ggplot2::element_text(size = 10),
    legend.text = ggplot2::element_text(size = 9)
  )

ggplot2::ggsave(
  file.path(
    out_dir,
    "Figure3G_CD8_selected_GO_BP.pdf"
  ),
  p_fig3g,
  width = 8.5,
  height = 6.5,
  useDingbats = FALSE
)

ggplot2::ggsave(
  file.path(
    out_dir,
    "Figure3G_CD8_selected_GO_BP.png"
  ),
  p_fig3g,
  width = 8.5,
  height = 6.5,
  dpi = 300
)

cat("\nFigure 3G selected GO terms:\n")
print(
  fig3g[
    ,
    c(
      "ID",
      "Description",
      "Label",
      "Direction",
      "SignedFE",
      "p.adjust",
      "Count"
    )
  ],
  row.names = FALSE
)

# ------------------------------------------------------------------------------
# 10. Final QC summary
# ------------------------------------------------------------------------------

qc <- data.frame(
  metric = c(
    "analysis_unit",
    "independent_core_pseudobulks",
    "NF1Mut_core_pseudobulks",
    "NF1WT_core_pseudobulks",
    "CD8_cells",
    "genes_entering_DESeq2",
    "significant_DE_genes_padj_lt_0.05",
    "Mut_up_significant_genes",
    "Mut_down_significant_genes",
    "significant_GO_BP_Mut_up",
    "significant_GO_BP_Mut_down",
    "Figure3G_terms_plotted"
  ),
  value = c(
    "tissue_core",
    ncol(pseudobulk_counts),
    sum(core_meta$NF1 == "Mut"),
    sum(core_meta$NF1 == "WT"),
    ncol(cd8),
    nrow(pseudobulk_counts),
    nrow(sig),
    length(up_genes),
    length(down_genes),
    nrow(go_up_sig),
    nrow(go_down_sig),
    nrow(fig3g)
  ),
  stringsAsFactors = FALSE
)

utils::write.csv(
  qc,
  file.path(
    out_dir,
    "CD8_core_pseudobulk_QC.csv"
  ),
  row.names = FALSE
)

cat("\nCD8 core-pseudobulk QC:\n")
print(qc, row.names = FALSE)

cat(
  "\nDONE: Figure 3F-G CD8 core-level pseudobulk analysis complete.\n"
)
