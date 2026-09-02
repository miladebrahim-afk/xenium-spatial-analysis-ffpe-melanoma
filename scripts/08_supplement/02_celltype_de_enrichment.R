# ==============================================================================
# 02_celltype_de_enrichment.R
#
# Figure S5A-D
# Cell type-specific NF1Mut vs NF1WT transcriptional programs in Xenium data.
#
# Recovered from the historical analysis scripts:
#   - Melanoma_De.R
#   - CAFs_DE_nov_25.R
#   - Myeloid_DE.R
#
# Panels
# ------
# S5A  Melanoma Hallmark GSEA
# S5B  Selected melanoma-gene heatmap (exact historical pheatmap scale='row')
# S5C  CAF GO:BP GSEA
# S5D  Myeloid selected GO:BP enrichment
#
# IMPORTANT
# ---------
# `patient_id` in this Xenium object is used operationally as the tissue-core /
# specimen identifier. Independent cores are not combined.
#
# The public outputs do not expose original tissue-core IDs.
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(Matrix)
  library(pheatmap)
  library(fgsea)
  library(msigdbr)
  library(gprofiler2)
})

# ------------------------------------------------------------------------------
# 0. Paths / constants
# ------------------------------------------------------------------------------

in_file <- file.path(
  "data",
  "processed",
  "merged_obj_annotated_with_meta_niche_published_reference.rds"
)

out_dir <- file.path(
  "results",
  "supplement",
  "figureS5_celltype_DE"
)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(in_file)) {
  stop("Missing input: ", in_file)
}

MELANOMA_TYPES <- c(
  "Melanoma Melanocytic",
  "Melanoma NC",
  "Melanoma Proliferative",
  "Melanoma Mesenchymal",
  "Melanoma Intermediate"
)

CAF_TYPES <- c(
  "CAFS Inflammatory",
  "CAFS Myofibroblast"
)

MYELOID_TYPES <- c(
  "Macrophages M2",
  "Macrophages IFN-γ–activated",
  "Activated Myeloid cells"
)

S5A_PATHWAYS <- c(
  "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
  "HALLMARK_UV_RESPONSE_DN",
  "HALLMARK_HYPOXIA",
  "HALLMARK_TGF_BETA_SIGNALING",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "HALLMARK_IL2_STAT5_SIGNALING",
  "HALLMARK_INFLAMMATORY_RESPONSE",
  "HALLMARK_INTERFERON_ALPHA_RESPONSE",
  "HALLMARK_IL6_JAK_STAT3_SIGNALING",
  "HALLMARK_INTERFERON_GAMMA_RESPONSE"
)

S5B_GENES <- c(
  "EGFR",
  "PDGFRB",
  "VEGFC",
  "TGFB1",
  "TGFB2",
  "NOTCH1",
  "MET",
  "VEGFA",
  "NGFR",
  "ERBB2",
  "TGFBR2",
  "TGFB3",
  "SOX2",
  "IFNL1",
  "IL10",
  "TNFSF18",
  "IFNA17",
  "IL21",
  "IFNG"
)

S5C_GO_TERMS <- c(
  "GOBP_RESPONSE_TO_EPIDERMAL_GROWTH_FACTOR",
  "GOBP_SIGNAL_TRANSDUCTION_IN_RESPONSE_TO_DNA_DAMAGE",
  "GOBP_DNA_REPLICATION",
  "GOBP_CELL_DIVISION",
  "GOBP_REGULATION_OF_CELL_CYCLE",
  "GOBP_ACTIVATION_OF_IMMUNE_RESPONSE",
  "GOBP_ADAPTIVE_IMMUNE_RESPONSE",
  "GOBP_CYTOKINE_MEDIATED_SIGNALING_PATHWAY",
  "GOBP_CYTOKINE_PRODUCTION",
  "GOBP_REGULATION_OF_INFLAMMATORY_RESPONSE",
  "GOBP_REGULATION_OF_T_CELL_ACTIVATION",
  "GOBP_T_CELL_ACTIVATION",
  "GOBP_REGULATION_OF_ANTIGEN_PROCESSING_AND_PRESENTATION"
)

S5D_UP_TERMS <- c(
  "cell surface receptor signaling pathway",
  "response to growth factor",
  "angiogenesis",
  "transforming growth factor beta receptor signaling pathway"
)

S5D_SELECTED_TERMS <- c(
  S5D_UP_TERMS,
  "defense response",
  "regulation of immune system process",
  "immune response",
  "inflammatory response",
  "response to cytokine",
  "programmed cell death",
  "T cell activation",
  "immune effector process"
)

# ------------------------------------------------------------------------------
# Helpers
# ------------------------------------------------------------------------------

normalize_nf1 <- function(x) {
  x <- trimws(as.character(x))

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

subset_celltypes <- function(obj, types) {
  md <- obj@meta.data

  if (!"celltype" %in% colnames(md)) {
    stop("The merged object does not contain metadata column `celltype`.")
  }

  keep <- rownames(md)[as.character(md$celltype) %in% types]

  if (length(keep) == 0L) {
    stop(
      "No cells found for requested cell type(s): ",
      paste(types, collapse = ", ")
    )
  }

  subset(obj, cells = keep)
}

prepare_nf1_object <- function(obj) {
  obj$NF1_clean <- normalize_nf1(obj$NF1)
  obj$NF1_clean <- factor(obj$NF1_clean, levels = c("WT", "Mut"))
  Idents(obj) <- "NF1_clean"
  DefaultAssay(obj) <- "RNA"

  obj <- NormalizeData(
    obj,
    assay = "RNA",
    normalization.method = "LogNormalize",
    scale.factor = 10000,
    verbose = FALSE
  )

  obj
}

flatten_leading_edge <- function(x) {
  vapply(
    x,
    function(z) paste(z, collapse = "/"),
    character(1)
  )
}

clean_hallmark_label <- function(x) {
  x <- sub("^HALLMARK_", "", x)
  gsub("_", " ", x)
}

clean_gobp_label <- function(x) {
  x <- sub("^GOBP_", "", x)
  x <- gsub("_", " ", x)
  tools::toTitleCase(tolower(x))
}

get_hallmark_sets <- function() {
  m <- msigdbr::msigdbr(
    species = "Homo sapiens",
    collection = "H"
  )

  split(m$gene_symbol, m$gs_name)
}

get_gobp_sets <- function() {
  m <- msigdbr::msigdbr(
    species = "Homo sapiens",
    collection = "C5",
    subcollection = "GO:BP"
  )

  split(m$gene_symbol, m$gs_name)
}

run_fgsea_10k <- function(pathways, ranks) {
  ranks <- ranks[
    !is.na(ranks) &
      !is.na(names(ranks)) &
      names(ranks) != ""
  ]

  # Historical analyses used one value per gene.
  rank_df <- tibble::tibble(
    gene = names(ranks),
    rank = as.numeric(ranks)
  ) |>
    dplyr::group_by(gene) |>
    dplyr::summarise(
      rank = mean(rank, na.rm = TRUE),
      .groups = "drop"
    )

  ranks <- rank_df$rank
  names(ranks) <- rank_df$gene
  ranks <- sort(ranks, decreasing = TRUE)

  set.seed(42)

  fgsea::fgsea(
    pathways = pathways,
    stats = ranks,
    minSize = 10,
    maxSize = 500,
    nperm = 10000
  ) |>
    as.data.frame()
}

save_gg <- function(p, prefix, width = 7, height = 6) {
  ggplot2::ggsave(
    file.path(out_dir, paste0(prefix, ".pdf")),
    p,
    width = width,
    height = height,
    useDingbats = FALSE
  )

  ggplot2::ggsave(
    file.path(out_dir, paste0(prefix, ".png")),
    p,
    width = width,
    height = height,
    dpi = 300
  )
}

# ------------------------------------------------------------------------------
# 1. Load publication-reference object
# ------------------------------------------------------------------------------

cat("\nLoading publication-reference merged object...\n")
obj <- readRDS(in_file)

required_meta <- c("celltype", "NF1", "patient_id")
missing_meta <- setdiff(required_meta, colnames(obj@meta.data))

if (length(missing_meta) > 0L) {
  stop(
    "Missing required metadata: ",
    paste(missing_meta, collapse = ", ")
  )
}

cat(
  "Loaded ",
  format(ncol(obj), big.mark = ","),
  " cells.\n",
  sep = ""
)

# ==============================================================================
# FIGURE S5A / S5B — MELANOMA
# ==============================================================================

cat("\n============================================================\n")
cat("Figure S5A/S5B: melanoma analysis\n")
cat("============================================================\n")

mel <- subset_celltypes(
  obj,
  MELANOMA_TYPES
)

mel <- prepare_nf1_object(mel)

cat(
  "Melanoma cells: ",
  format(ncol(mel), big.mark = ","),
  "\n",
  sep = ""
)

print(table(mel$NF1_clean))

# ------------------------------------------------------------------------------
# S5A. Melanoma DE + Hallmark GSEA
#
# Historical settings:
#   RNA LogNormalize(scale.factor = 10000)
#   FindMarkers Mut vs WT
#   Wilcoxon
#   logfc.threshold = 0
#   min.pct = 0
#   only.pos = FALSE
#   rank ALL genes by avg_log2FC
#   Hallmark fgsea, nperm = 10000
# ------------------------------------------------------------------------------

cat("\nRunning melanoma Mut-vs-WT DE...\n")

mel_de <- FindMarkers(
  mel,
  ident.1 = "Mut",
  ident.2 = "WT",
  assay = "RNA",
  test.use = "wilcox",
  logfc.threshold = 0,
  min.pct = 0,
  only.pos = FALSE,
  verbose = FALSE
)

mel_de_df <- mel_de |>
  tibble::rownames_to_column("gene") |>
  dplyr::arrange(p_val_adj)

utils::write.csv(
  mel_de_df,
  file.path(
    out_dir,
    "FigureS5A_melanoma_DE_Mut_vs_WT.csv"
  ),
  row.names = FALSE
)

mel_ranks <- mel_de_df$avg_log2FC
names(mel_ranks) <- mel_de_df$gene

cat("Running melanoma Hallmark GSEA...\n")

hallmark_sets <- get_hallmark_sets()

mel_gsea <- run_fgsea_10k(
  hallmark_sets,
  mel_ranks
)

mel_gsea_csv <- mel_gsea

if ("leadingEdge" %in% colnames(mel_gsea_csv)) {
  mel_gsea_csv$leadingEdge <- flatten_leading_edge(
    mel_gsea_csv$leadingEdge
  )
}

utils::write.csv(
  mel_gsea_csv,
  file.path(
    out_dir,
    "FigureS5A_melanoma_Hallmark_GSEA_all.csv"
  ),
  row.names = FALSE
)

s5a <- mel_gsea |>
  dplyr::filter(
    pathway %in% S5A_PATHWAYS
  ) |>
  dplyr::mutate(
    gene_count = lengths(leadingEdge),
    padj_plot = pmax(
      padj,
      .Machine$double.xmin,
      na.rm = FALSE
    ),
    negLogPadj = -log10(padj_plot),
    pathway_label = clean_hallmark_label(pathway),
    pathway_label = factor(
      pathway_label,
      levels = rev(
        clean_hallmark_label(S5A_PATHWAYS)
      )
    )
  )

missing_s5a <- setdiff(
  S5A_PATHWAYS,
  s5a$pathway
)

if (length(missing_s5a) > 0L) {
  warning(
    "S5A requested pathway(s) not returned by fgsea: ",
    paste(missing_s5a, collapse = ", ")
  )
}

utils::write.csv(
  s5a |>
    dplyr::mutate(
      leadingEdge = flatten_leading_edge(leadingEdge)
    ),
  file.path(
    out_dir,
    "FigureS5A_selected_Hallmark_pathways.csv"
  ),
  row.names = FALSE
)

p_s5a <- ggplot2::ggplot(
  s5a,
  ggplot2::aes(
    x = NES,
    y = pathway_label
  )
) +
  ggplot2::geom_point(
    ggplot2::aes(
      size = gene_count,
      color = negLogPadj
    )
  ) +
  ggplot2::geom_vline(
    xintercept = 0,
    linetype = "dashed"
  ) +
  ggplot2::scale_color_gradient(
    low = "blue",
    high = "red",
    name = expression(-log[10]("FDR"))
  ) +
  ggplot2::labs(
    x = "Normalized Enrichment Score (NES)",
    y = NULL,
    size = "Leading-edge\ngenes",
    title = "Melanoma: NF1Mut vs NF1WT"
  ) +
  ggplot2::theme_classic(base_size = 13) +
  ggplot2::theme(
    axis.text.y = ggplot2::element_text(size = 10),
    legend.position = "right"
  )

print(p_s5a)

save_gg(
  p_s5a,
  "FigureS5A_melanoma_Hallmark_GSEA_FINAL",
  width = 7,
  height = 6
)

# ------------------------------------------------------------------------------
# S5B. Selected melanoma-gene heatmap — EXACT HISTORICAL LOGIC
#
# Recovered from the original Melanoma_De.R:
#   1) selected RNA raw counts
#   2) log1p(counts)
#   3) mean expression per tissue core
#   4) order columns NF1Mut then NF1WT
#   5) pheatmap(scale = "row")
#   6) blue -> white -> red, cluster rows, do not cluster columns
#
# IMPORTANT:
# Do NOT pre-scale the matrix and then use scale = "none".
# The historical pheatmap call itself performs row z-scoring.
# ------------------------------------------------------------------------------

cat("\nBuilding Figure S5B selected-gene heatmap (exact historical logic)...\n")

genes_present <- intersect(
  S5B_GENES,
  rownames(mel[["RNA"]])
)

genes_missing <- setdiff(
  S5B_GENES,
  genes_present
)

if (length(genes_missing) > 0L) {
  warning(
    "S5B genes absent from RNA assay: ",
    paste(genes_missing, collapse = ", ")
  )
}

if (length(genes_present) == 0L) {
  stop("None of the selected S5B genes are present in the RNA assay.")
}

# Historical code:
# expr_counts <- GetAssayData(..., layer = "counts")[genes_present, ]
# expr_log <- log1p(expr_counts)

expr_counts <- SeuratObject::LayerData(
  mel,
  assay = "RNA",
  layer = "counts"
)[
  genes_present,
  ,
  drop = FALSE
]

expr_log <- log1p(expr_counts)

core_id <- as.character(mel$patient_id)
core_nf1 <- as.character(mel$NF1_clean)

core_info <- data.frame(
  patient_id = core_id,
  NF1 = core_nf1,
  stringsAsFactors = FALSE
) |>
  dplyr::distinct()

# Ensure one NF1 label per tissue core.
nf1_conflicts <- core_info |>
  dplyr::count(patient_id) |>
  dplyr::filter(n > 1)

if (nrow(nf1_conflicts) > 0L) {
  stop("At least one tissue core has conflicting NF1 labels.")
}

# Historical final order: Mut first, WT second.
core_info <- core_info |>
  dplyr::mutate(
    NF1 = factor(
      NF1,
      levels = c("Mut", "WT")
    )
  ) |>
  dplyr::arrange(
    NF1,
    patient_id
  )

core_order <- core_info$patient_id

cat("Averaging selected-gene expression within each tissue core...\n")

expr_agg_patient <- vapply(
  core_order,
  function(core) {
    idx <- which(core_id == core)

    Matrix::rowMeans(
      expr_log[, idx, drop = FALSE]
    )
  },
  numeric(length(genes_present))
)

rownames(expr_agg_patient) <- genes_present
colnames(expr_agg_patient) <- core_order

# Historical NF1 annotation.
annotation_col <- data.frame(
  NF1 = factor(
    as.character(core_info$NF1),
    levels = c("Mut", "WT")
  ),
  row.names = core_order,
  check.names = FALSE
)

ann_colors <- list(
  NF1 = c(
    Mut = "red",
    WT = "blue"
  )
)

heatmap_colors <- grDevices::colorRampPalette(
  c("blue", "white", "red")
)(50)

pdf_file_s5b <- file.path(
  out_dir,
  "FigureS5B_selected_melanoma_genes_heatmap_FINAL.pdf"
)

png_file_s5b <- file.path(
  out_dir,
  "FigureS5B_selected_melanoma_genes_heatmap_FINAL.png"
)

cat("Writing Figure S5B PDF...\n")

pheatmap::pheatmap(
  expr_agg_patient,
  scale = "row",
  color = heatmap_colors,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  annotation_col = annotation_col,
  annotation_colors = ann_colors,
  annotation_legend = FALSE,
  show_rownames = TRUE,
  show_colnames = FALSE,
  main = "Selected melanoma genes",
  filename = pdf_file_s5b,
  width = 6,
  height = 6
)

cat("Writing Figure S5B PNG...\n")

pheatmap::pheatmap(
  expr_agg_patient,
  scale = "row",
  color = heatmap_colors,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  annotation_col = annotation_col,
  annotation_colors = ann_colors,
  annotation_legend = FALSE,
  show_rownames = TRUE,
  show_colnames = FALSE,
  main = "Selected melanoma genes",
  filename = png_file_s5b,
  width = 6,
  height = 6
)

if (
  !file.exists(png_file_s5b) ||
  file.info(png_file_s5b)$size == 0
) {
  stop("Figure S5B PNG was not written correctly.")
}

# Save a PUBLIC anonymous row-z-score matrix for reproducibility.
# This is derived only after reproducing the historical pheatmap input.
z_s5b <- t(
  scale(
    t(expr_agg_patient)
  )
)

anon_names <- sprintf(
  "%s_%02d",
  as.character(core_info$NF1),
  ave(
    seq_len(nrow(core_info)),
    as.character(core_info$NF1),
    FUN = seq_along
  )
)

colnames(z_s5b) <- anon_names

utils::write.csv(
  data.frame(
    gene = rownames(z_s5b),
    z_s5b,
    check.names = FALSE
  ),
  file.path(
    out_dir,
    "FigureS5B_selected_gene_heatmap_rowZ_anonymous.csv"
  ),
  row.names = FALSE
)

observed_z_range <- range(
  z_s5b,
  finite = TRUE
)

cat(
  "Figure S5B observed row-z-score range: ",
  round(observed_z_range[1], 3),
  " to ",
  round(observed_z_range[2], 3),
  "\n",
  sep = ""
)

cat(
  "Figure S5B PNG size: ",
  round(file.info(png_file_s5b)$size / 1024^2, 3),
  " MB\n",
  sep = ""
)

rm(
  expr_counts,
  expr_log,
  expr_agg_patient,
  z_s5b,
  mel,
  mel_de,
  mel_de_df,
  mel_ranks,
  mel_gsea,
  mel_gsea_csv,
  s5a
)
gc(verbose = FALSE)

# ==============================================================================
# FIGURE S5C — CAFs
# ==============================================================================

cat("\n============================================================\n")
cat("Figure S5C: CAF GO:BP enrichment\n")
cat("============================================================\n")

caf <- subset_celltypes(
  obj,
  CAF_TYPES
)

caf <- prepare_nf1_object(caf)

cat(
  "CAF cells: ",
  format(ncol(caf), big.mark = ","),
  "\n",
  sep = ""
)

print(table(caf$NF1_clean))

# Historical CAF DE settings:
# logfc.threshold = 0
# min.pct = 0
# Wilcoxon
cat("\nRunning CAF Mut-vs-WT DE...\n")

caf_de <- FindMarkers(
  caf,
  ident.1 = "Mut",
  ident.2 = "WT",
  assay = "RNA",
  test.use = "wilcox",
  logfc.threshold = 0,
  min.pct = 0,
  only.pos = FALSE,
  verbose = FALSE
)

caf_de_df <- caf_de |>
  tibble::rownames_to_column("gene") |>
  dplyr::arrange(p_val_adj)

utils::write.csv(
  caf_de_df,
  file.path(
    out_dir,
    "FigureS5C_CAF_DE_Mut_vs_WT.csv"
  ),
  row.names = FALSE
)

# Historical CAF GSEA block ranked the significant DE genes.
caf_sig <- caf_de_df |>
  dplyr::filter(
    p_val_adj < 0.05
  ) |>
  dplyr::arrange(
    dplyr::desc(avg_log2FC)
  )

if (nrow(caf_sig) == 0L) {
  stop("No significant CAF genes available for historical S5C GSEA.")
}

caf_ranks <- caf_sig$avg_log2FC
names(caf_ranks) <- caf_sig$gene

cat("Running CAF GO:BP GSEA...\n")

gobp_sets <- get_gobp_sets()

caf_gsea <- run_fgsea_10k(
  gobp_sets,
  caf_ranks
)

caf_gsea_csv <- caf_gsea

if ("leadingEdge" %in% colnames(caf_gsea_csv)) {
  caf_gsea_csv$leadingEdge <- flatten_leading_edge(
    caf_gsea_csv$leadingEdge
  )
}

utils::write.csv(
  caf_gsea_csv,
  file.path(
    out_dir,
    "FigureS5C_CAF_GO_BP_GSEA_all.csv"
  ),
  row.names = FALSE
)

s5c <- caf_gsea |>
  dplyr::filter(
    pathway %in% S5C_GO_TERMS
  ) |>
  dplyr::mutate(
    padj_plot = pmax(
      padj,
      .Machine$double.xmin,
      na.rm = FALSE
    ),
    negLogPadj = -log10(padj_plot),
    pathway_label = clean_gobp_label(pathway)
  )

missing_s5c <- setdiff(
  S5C_GO_TERMS,
  s5c$pathway
)

if (length(missing_s5c) > 0L) {
  warning(
    "S5C requested GO term(s) not returned by fgsea: ",
    paste(missing_s5c, collapse = ", ")
  )
}

utils::write.csv(
  s5c |>
    dplyr::mutate(
      leadingEdge = flatten_leading_edge(leadingEdge)
    ),
  file.path(
    out_dir,
    "FigureS5C_selected_CAF_GO_BP.csv"
  ),
  row.names = FALSE
)

p_s5c <- ggplot2::ggplot(
  s5c,
  ggplot2::aes(
    x = NES,
    y = reorder(pathway_label, NES)
  )
) +
  ggplot2::geom_point(
    ggplot2::aes(
      size = size,
      color = negLogPadj
    )
  ) +
  ggplot2::geom_vline(
    xintercept = 0,
    linetype = "dashed"
  ) +
  ggplot2::scale_color_gradient(
    low = "blue",
    high = "red",
    name = expression(-log[10]("adjusted p"))
  ) +
  ggplot2::labs(
    x = "Normalized Enrichment Score (NES)",
    y = NULL,
    size = "Gene-set size",
    title = "CAFs: NF1Mut vs NF1WT"
  ) +
  ggplot2::theme_classic(base_size = 13) +
  ggplot2::theme(
    axis.text.y = ggplot2::element_text(size = 9),
    legend.position = "right"
  )

print(p_s5c)

save_gg(
  p_s5c,
  "FigureS5C_CAF_GO_BP_FINAL",
  width = 8,
  height = 7
)

rm(
  caf,
  caf_de,
  caf_de_df,
  caf_sig,
  caf_ranks,
  caf_gsea,
  caf_gsea_csv,
  s5c
)
gc(verbose = FALSE)

# ==============================================================================
# FIGURE S5D — MYELOID
# ==============================================================================

cat("\n============================================================\n")
cat("Figure S5D: myeloid GO:BP enrichment\n")
cat("============================================================\n")

myeloid <- subset_celltypes(
  obj,
  MYELOID_TYPES
)

myeloid <- prepare_nf1_object(myeloid)

cat(
  "Myeloid cells: ",
  format(ncol(myeloid), big.mark = ","),
  "\n",
  sep = ""
)

print(table(myeloid$NF1_clean))

# Historical executed Myeloid_DE.R settings:
# logfc.threshold = 0.25
# min.pct = 0.1
# Wilcoxon
cat("\nRunning myeloid Mut-vs-WT DE...\n")

myeloid_de <- FindMarkers(
  myeloid,
  ident.1 = "Mut",
  ident.2 = "WT",
  assay = "RNA",
  test.use = "wilcox",
  logfc.threshold = 0.25,
  min.pct = 0.1,
  only.pos = FALSE,
  verbose = FALSE
)

myeloid_de_df <- myeloid_de |>
  tibble::rownames_to_column("gene") |>
  dplyr::arrange(p_val_adj)

utils::write.csv(
  myeloid_de_df,
  file.path(
    out_dir,
    "FigureS5D_myeloid_DE_Mut_vs_WT.csv"
  ),
  row.names = FALSE
)

myeloid_sig <- myeloid_de_df |>
  dplyr::filter(
    p_val_adj < 0.05
  ) |>
  dplyr::arrange(
    dplyr::desc(avg_log2FC)
  )

up_genes <- myeloid_sig |>
  dplyr::filter(
    avg_log2FC > 0
  ) |>
  dplyr::pull(gene)

down_genes <- myeloid_sig |>
  dplyr::filter(
    avg_log2FC < 0
  ) |>
  dplyr::pull(gene)

if (length(up_genes) == 0L || length(down_genes) == 0L) {
  stop(
    "Historical S5D enrichment requires both significant up- and downregulated myeloid genes."
  )
}

# Historical final block converted symbols to ENSG before g:Profiler.
cat(
  "Running g:Profiler GO:BP enrichment for S5D...\n",
  "This step requires network access from the HPC.\n",
  sep = ""
)

conv_up <- gprofiler2::gconvert(
  up_genes,
  organism = "hsapiens",
  target = "ENSG"
)

conv_down <- gprofiler2::gconvert(
  down_genes,
  organism = "hsapiens",
  target = "ENSG"
)

up_mapped <- unique(
  conv_up$target[
    !is.na(conv_up$target) &
      conv_up$target != ""
  ]
)

down_mapped <- unique(
  conv_down$target[
    !is.na(conv_down$target) &
      conv_down$target != ""
  ]
)

if (length(up_mapped) == 0L || length(down_mapped) == 0L) {
  stop("g:Profiler gene conversion returned no mapped genes.")
}

gost_up <- gprofiler2::gost(
  query = up_mapped,
  organism = "hsapiens",
  sources = "GO:BP",
  user_threshold = 0.05,
  correction_method = "fdr",
  evcodes = TRUE
)

gost_down <- gprofiler2::gost(
  query = down_mapped,
  organism = "hsapiens",
  sources = "GO:BP",
  user_threshold = 0.05,
  correction_method = "fdr",
  evcodes = TRUE
)

if (
  is.null(gost_up$result) ||
  is.null(gost_down$result)
) {
  stop("g:Profiler returned no GO:BP enrichment table.")
}

# Preserve historical final selection logic:
# bind UP results first, then DOWN; keep one occurrence per selected term;
# manually assign the four growth-factor programs as Upregulated and the
# remaining selected programs as Downregulated.
s5d <- dplyr::bind_rows(
  gost_up$result,
  gost_down$result
) |>
  dplyr::filter(
    term_name %in% S5D_SELECTED_TERMS
  ) |>
  dplyr::distinct(
    term_name,
    .keep_all = TRUE
  ) |>
  dplyr::mutate(
    Direction = ifelse(
      term_name %in% S5D_UP_TERMS,
      "Upregulated",
      "Downregulated"
    ),
    log10p = -log10(
      pmax(
        p_value,
        .Machine$double.xmin
      )
    ),
    SignedFE = ifelse(
      Direction == "Upregulated",
      intersection_size / query_size,
      -(intersection_size / query_size)
    )
  ) |>
  dplyr::arrange(
    dplyr::desc(abs(SignedFE))
  )

missing_s5d <- setdiff(
  S5D_SELECTED_TERMS,
  s5d$term_name
)

if (length(missing_s5d) > 0L) {
  warning(
    "S5D selected GO term(s) not returned by current g:Profiler: ",
    paste(missing_s5d, collapse = ", ")
  )
}

# Save only summary fields needed for reproducibility; no tissue-core IDs.
s5d_public <- s5d |>
  dplyr::select(
    term_id,
    term_name,
    p_value,
    intersection_size,
    query_size,
    term_size,
    Direction,
    log10p,
    SignedFE
  )

utils::write.csv(
  s5d_public,
  file.path(
    out_dir,
    "FigureS5D_selected_myeloid_GO_BP.csv"
  ),
  row.names = FALSE
)

p_s5d <- ggplot2::ggplot(
  s5d,
  ggplot2::aes(
    x = SignedFE,
    y = reorder(term_name, SignedFE)
  )
) +
  ggplot2::geom_point(
    ggplot2::aes(
      color = log10p,
      size = intersection_size
    )
  ) +
  ggplot2::geom_vline(
    xintercept = 0,
    linetype = "dashed"
  ) +
  ggplot2::scale_color_gradient(
    low = "blue",
    high = "red",
    name = expression(-log[10](p))
  ) +
  ggplot2::scale_size_continuous(
    name = "Gene count"
  ) +
  ggplot2::labs(
    x = "Signed enrichment",
    y = NULL,
    title = "Myeloid: NF1Mut vs NF1WT"
  ) +
  ggplot2::theme_classic(base_size = 13) +
  ggplot2::theme(
    axis.text.y = ggplot2::element_text(size = 9),
    legend.position = "right"
  )

print(p_s5d)

save_gg(
  p_s5d,
  "FigureS5D_myeloid_GO_BP_FINAL",
  width = 8,
  height = 7
)

# ------------------------------------------------------------------------------
# Final QC
# ------------------------------------------------------------------------------

cat("\n============================================================\n")
cat("Figure S5 analysis complete\n")
cat("============================================================\n")
cat("Outputs written to: ", out_dir, "\n", sep = "")
cat("S5A: melanoma Hallmark GSEA\n")
cat("S5B: selected melanoma-gene heatmap\n")
cat("S5C: CAF GO:BP GSEA\n")
cat("S5D: myeloid selected GO:BP enrichment\n")
