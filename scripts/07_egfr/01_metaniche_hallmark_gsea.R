# ==============================================================================
# 01_metaniche_hallmark_gsea.R
#
# Figure 5A: NF1Mut vs NF1WT Hallmark pathway enrichment within each meta-niche.
#
# PROVENANCE
# ----------
# Clean reconstruction of the surviving historical script MN_Dot_Mut_VS_WT.R.
# The historical analysis:
#   1. analyzed each meta-niche separately;
#   2. summed raw RNA counts within each tissue core;
#   3. treated each tissue core as an independent pseudobulk sample;
#   4. filtered genes to those detected in >= 2 pseudobulk samples;
#   5. used edgeR normalization + limma-voom with design ~ NF1;
#   6. tested coef = "NF1Mut";
#   7. ranked genes by limma logFC;
#   8. ran Hallmark fgsea with 1,000 permutations.
#
# IMPORTANT ANALYSIS UNIT
# -----------------------
# `patient_id` is used operationally as the independent Xenium tissue-core ID.
# Paired-looking core IDs (for example, an ID with a terminal "_2") are NOT
# merged. This matches the analysis decision used for the reconstructed paper.
#
# INPUT
# -----
# data/processed/merged_obj_annotated_with_meta_niche_published_reference.rds
#
# OUTPUTS
# -------
# results/egfr/figure5A_metaniche_hallmark/
#   DE/
#   GSEA/
#   Figure5A_sample_QC_by_meta_niche.csv
#   Figure5A_all_limma_DE.csv
#   Figure5A_all_Hallmark_GSEA.csv
#   Figure5A_selected_Hallmark_GSEA.csv
#   Figure5A_selected_Hallmark_dotplot.pdf/png
#
# INTERPRETATION
# --------------
# Positive logFC / NES = enriched in NF1Mut relative to NF1WT.
# Negative logFC / NES = enriched in NF1WT relative to NF1Mut.
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(edgeR)
  library(limma)
  library(fgsea)
  library(msigdbr)
  library(dplyr)
  library(tibble)
  library(ggplot2)
  library(stringr)
})

# ------------------------------------------------------------------------------
# 0. Paths and constants
# ------------------------------------------------------------------------------

in_file <- file.path(
  "data",
  "processed",
  "merged_obj_annotated_with_meta_niche_published_reference.rds"
)

out_dir <- file.path(
  "results",
  "egfr",
  "figure5A_metaniche_hallmark"
)

de_dir <- file.path(out_dir, "DE")
gsea_dir <- file.path(out_dir, "GSEA")

dir.create(de_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(gsea_dir, recursive = TRUE, showWarnings = FALSE)

MN_LEVELS <- paste0("MetaNiche_", 1:12)

# Final manuscript-facing meta-niche labels.
MN_LABELS <- c(
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

# Historical display order:
# tumor-rich -> immune-rich -> stromal-rich -> vascular.
MN_DISPLAY_ORDER <- c(
  "MetaNiche_1",
  "MetaNiche_2",
  "MetaNiche_5",
  "MetaNiche_9",
  "MetaNiche_10",
  "MetaNiche_12",
  "MetaNiche_3",
  "MetaNiche_11",
  "MetaNiche_6",
  "MetaNiche_7",
  "MetaNiche_8",
  "MetaNiche_4"
)

# Exact selected Hallmark set retained from the historical Figure 5A script.
SELECTED_PATHWAYS <- c(
  "HALLMARK_WNT_BETA_CATENIN_SIGNALING",
  "HALLMARK_UV_RESPONSE_UP",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "HALLMARK_TGF_BETA_SIGNALING",
  "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY",
  "HALLMARK_PI3K_AKT_MTOR_SIGNALING",
  "HALLMARK_P53_PATHWAY",
  "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
  "HALLMARK_NOTCH_SIGNALING",
  "HALLMARK_MTORC1_SIGNALING",
  "HALLMARK_KRAS_SIGNALING_UP",
  "HALLMARK_KRAS_SIGNALING_DN",
  "HALLMARK_INTERFERON_ALPHA_RESPONSE",
  "HALLMARK_INFLAMMATORY_RESPONSE",
  "HALLMARK_IL6_JAK_STAT3_SIGNALING",
  "HALLMARK_IL2_STAT5_SIGNALING",
  "HALLMARK_HYPOXIA",
  "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
  "HALLMARK_E2F_TARGETS",
  "HALLMARK_APOPTOSIS",
  "HALLMARK_APICAL_JUNCTION"
)

# Historical analysis constants.
MIN_PSEUDOBULKS_GENE_DETECTED <- 2L
MIN_GENES_AFTER_FILTER <- 10L
FGSEA_NPERM <- 1000L

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
      paste(sort(unique(x[is.na(out)])), collapse = ", ")
    )
  }

  out
}

get_hallmark_sets <- function() {

  # msigdbr changed its preferred argument names across versions.
  # Try the current API first, then fall back to the historical API.
  h <- tryCatch(
    msigdbr::msigdbr(
      species = "Homo sapiens",
      collection = "H"
    ),
    error = function(e) {
      msigdbr::msigdbr(
        species = "Homo sapiens",
        category = "H"
      )
    }
  )

  required <- c("gs_name", "gene_symbol")

  if (!all(required %in% colnames(h))) {
    stop(
      "msigdbr Hallmark table is missing: ",
      paste(setdiff(required, colnames(h)), collapse = ", ")
    )
  }

  split(
    x = as.character(h$gene_symbol),
    f = as.character(h$gs_name)
  )
}

flatten_leading_edge <- function(x) {
  if (is.list(x)) {
    vapply(
      x,
      function(z) paste(z, collapse = ";"),
      character(1)
    )
  } else {
    x
  }
}

# ------------------------------------------------------------------------------
# 2. Load publication-reference object
# ------------------------------------------------------------------------------

if (!file.exists(in_file)) {
  stop("Missing input: ", in_file)
}

cat("Loading publication-reference merged object...\n")
obj <- readRDS(in_file)

required_meta <- c(
  "patient_id",
  "NF1",
  "meta_niche"
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

DefaultAssay(obj) <- "RNA"

rna_layers <- SeuratObject::Layers(obj[["RNA"]])

cat(
  "RNA layers: ",
  paste(rna_layers, collapse = ", "),
  "\n",
  sep = ""
)

if ("counts" %in% rna_layers) {

  raw_counts <- SeuratObject::LayerData(
    obj,
    assay = "RNA",
    layer = "counts"
  )

} else {

  count_layers <- grep(
    "^counts",
    rna_layers,
    value = TRUE
  )

  if (length(count_layers) == 0L) {
    stop("No RNA counts layer found.")
  }

  cat(
    "Joining ",
    length(count_layers),
    " RNA counts layers...\n",
    sep = ""
  )

  obj[["RNA"]] <- SeuratObject::JoinLayers(
    obj[["RNA"]],
    layers = count_layers,
    new = "counts"
  )

  raw_counts <- SeuratObject::LayerData(
    obj,
    assay = "RNA",
    layer = "counts"
  )
}

md <- obj@meta.data |>
  as.data.frame()

md$sample_id <- as.character(md$patient_id)
md$NF1_clean <- normalize_nf1(md$NF1)
md$meta_niche_clean <- as.character(md$meta_niche)

if (!identical(colnames(raw_counts), rownames(md))) {
  stop("RNA count columns do not align with metadata rows.")
}

cat(
  "Merged cells: ",
  format(ncol(obj), big.mark = ","),
  "\n",
  "Cells with published meta-niche: ",
  format(sum(md$meta_niche_clean %in% MN_LEVELS), big.mark = ","),
  "\n",
  "Independent tissue-core IDs: ",
  dplyr::n_distinct(md$sample_id),
  "\n",
  sep = ""
)

rm(obj)
gc(verbose = FALSE)

# ------------------------------------------------------------------------------
# 3. Load Hallmark gene sets
# ------------------------------------------------------------------------------

cat("\nLoading Hallmark gene sets from msigdbr...\n")
hallmark_sets <- get_hallmark_sets()

cat(
  "Hallmark pathways loaded: ",
  length(hallmark_sets),
  "\n",
  sep = ""
)

# ------------------------------------------------------------------------------
# 4. Per-meta-niche core pseudobulk + limma-voom + fgsea
# ------------------------------------------------------------------------------

de_results <- list()
gsea_results <- list()
qc_rows <- list()

for (mn in MN_LEVELS) {

  cat(
    "\n============================================================\n",
    "Processing ",
    mn,
    " | ",
    MN_LABELS[[mn]],
    "\n",
    sep = ""
  )

  cell_idx <- which(
    md$meta_niche_clean == mn
  )

  if (length(cell_idx) == 0L) {
    warning("No cells found for ", mn, "; skipping.")
    next
  }

  sample_ids <- md$sample_id[cell_idx]
  nf1_values <- md$NF1_clean[cell_idx]

  samples_present <- unique(sample_ids)

  # Verify every core has exactly one NF1 label.
  sample_nf1 <- vapply(
    samples_present,
    function(sid) {
      z <- unique(
        nf1_values[sample_ids == sid]
      )

      if (length(z) != 1L) {
        stop(
          "Core ",
          sid,
          " has ",
          length(z),
          " NF1 labels within ",
          mn
        )
      }

      z
    },
    character(1)
  )

  # Preserve independent core IDs exactly; never strip terminal suffixes.
  pseudobulk_list <- lapply(
    samples_present,
    function(sid) {

      cells_sid <- cell_idx[
        sample_ids == sid
      ]

      Matrix::rowSums(
        raw_counts[
          ,
          cells_sid,
          drop = FALSE
        ]
      )
    }
  )

  pseudobulk_counts <- do.call(
    cbind,
    pseudobulk_list
  )

  colnames(pseudobulk_counts) <- samples_present

  coldata <- data.frame(
    sample_id = samples_present,
    NF1 = sample_nf1,
    stringsAsFactors = FALSE
  )

  rownames(coldata) <- coldata$sample_id

  coldata <- coldata[
    colnames(pseudobulk_counts),
    ,
    drop = FALSE
  ]

  coldata$NF1 <- factor(
    coldata$NF1,
    levels = c("WT", "Mut")
  )

  n_mut <- sum(coldata$NF1 == "Mut")
  n_wt <- sum(coldata$NF1 == "WT")

  cat(
    "  cells: ",
    format(length(cell_idx), big.mark = ","),
    "\n",
    "  pseudobulk cores: ",
    ncol(pseudobulk_counts),
    " (Mut=",
    n_mut,
    ", WT=",
    n_wt,
    ")\n",
    sep = ""
  )

  # A two-group NF1 comparison requires both genotypes.
  if (n_mut < 1L || n_wt < 1L) {

    warning(
      mn,
      " lacks one NF1 group and cannot be modeled; skipping."
    )

    qc_rows[[mn]] <- data.frame(
      meta_niche = mn,
      meta_niche_label = MN_LABELS[[mn]],
      n_cells = length(cell_idx),
      n_cores = ncol(pseudobulk_counts),
      n_mut = n_mut,
      n_wt = n_wt,
      genes_before_filter = nrow(pseudobulk_counts),
      genes_after_filter = NA_integer_,
      status = "skipped_missing_NF1_group",
      stringsAsFactors = FALSE
    )

    next
  }

  # Exact historical gene filter:
  # retain genes with non-zero counts in at least 2 pseudobulk samples.
  keep_genes <- Matrix::rowSums(
    pseudobulk_counts > 0
  ) >= MIN_PSEUDOBULKS_GENE_DETECTED

  genes_before <- nrow(pseudobulk_counts)

  pseudobulk_counts <- pseudobulk_counts[
    keep_genes,
    ,
    drop = FALSE
  ]

  genes_after <- nrow(pseudobulk_counts)

  cat(
    "  genes retained: ",
    genes_after,
    "/",
    genes_before,
    "\n",
    sep = ""
  )

  if (genes_after < MIN_GENES_AFTER_FILTER) {

    warning(
      mn,
      " has fewer than ",
      MIN_GENES_AFTER_FILTER,
      " genes after filtering; skipping."
    )

    qc_rows[[mn]] <- data.frame(
      meta_niche = mn,
      meta_niche_label = MN_LABELS[[mn]],
      n_cells = length(cell_idx),
      n_cores = ncol(pseudobulk_counts),
      n_mut = n_mut,
      n_wt = n_wt,
      genes_before_filter = genes_before,
      genes_after_filter = genes_after,
      status = "skipped_too_few_genes",
      stringsAsFactors = FALSE
    )

    next
  }

  # --------------------------------------------------------------------------
  # limma-voom
  # --------------------------------------------------------------------------

  dge <- edgeR::DGEList(
    counts = pseudobulk_counts
  )

  dge <- edgeR::calcNormFactors(
    dge
  )

  design <- stats::model.matrix(
    ~ NF1,
    data = coldata
  )

  if (!"NF1Mut" %in% colnames(design)) {
    stop(
      "Expected design coefficient 'NF1Mut' was not created for ",
      mn,
      ". Design columns: ",
      paste(colnames(design), collapse = ", ")
    )
  }

  v <- limma::voom(
    dge,
    design,
    plot = FALSE
  )

  fit <- limma::lmFit(
    v,
    design
  )

  fit <- limma::eBayes(
    fit
  )

  res <- limma::topTable(
    fit,
    coef = "NF1Mut",
    number = Inf,
    sort.by = "P"
  )

  res <- tibble::rownames_to_column(
    as.data.frame(res),
    "gene"
  )

  res$meta_niche <- mn
  res$meta_niche_label <- MN_LABELS[[mn]]

  de_results[[mn]] <- res

  utils::write.csv(
    res,
    file.path(
      de_dir,
      paste0("DE_", mn, "_NF1Mut_vs_WT_limma_voom.csv")
    ),
    row.names = FALSE
  )

  # --------------------------------------------------------------------------
  # Hallmark fgsea
  # --------------------------------------------------------------------------

  ranks <- res$logFC
  names(ranks) <- res$gene

  keep_rank <- is.finite(ranks) & !is.na(names(ranks))

  ranks <- ranks[
    keep_rank
  ]

  # fgsea requires unique gene names.
  if (anyDuplicated(names(ranks))) {
    ranks <- tapply(
      ranks,
      names(ranks),
      function(z) z[which.max(abs(z))]
    )
  }

  ranks <- sort(
    ranks,
    decreasing = TRUE
  )

  set.seed(42)

  gsea_res <- fgsea::fgsea(
    pathways = hallmark_sets,
    stats = ranks,
    nperm = FGSEA_NPERM
  )

  gsea_df <- as.data.frame(
    gsea_res
  )

  gsea_df$meta_niche <- mn
  gsea_df$meta_niche_label <- MN_LABELS[[mn]]

  gsea_results[[mn]] <- gsea_df

  # Flatten leadingEdge only for CSV output.
  gsea_csv <- gsea_df
  if ("leadingEdge" %in% colnames(gsea_csv)) {
    gsea_csv$leadingEdge <- flatten_leading_edge(
      gsea_csv$leadingEdge
    )
  }

  utils::write.csv(
    gsea_csv,
    file.path(
      gsea_dir,
      paste0("GSEA_", mn, "_Hallmark.csv")
    ),
    row.names = FALSE
  )

  qc_rows[[mn]] <- data.frame(
    meta_niche = mn,
    meta_niche_label = MN_LABELS[[mn]],
    n_cells = length(cell_idx),
    n_cores = ncol(pseudobulk_counts),
    n_mut = n_mut,
    n_wt = n_wt,
    genes_before_filter = genes_before,
    genes_after_filter = genes_after,
    status = "modeled",
    stringsAsFactors = FALSE
  )

  cat(
    "  modeled successfully; Hallmark pathways tested: ",
    nrow(gsea_df),
    "\n",
    sep = ""
  )

  rm(
    pseudobulk_list,
    pseudobulk_counts,
    dge,
    v,
    fit,
    res,
    ranks,
    gsea_res,
    gsea_df
  )
  gc(verbose = FALSE)
}

# ------------------------------------------------------------------------------
# 5. Save combined QC / DE / GSEA tables
# ------------------------------------------------------------------------------

sample_qc <- dplyr::bind_rows(
  qc_rows
)

utils::write.csv(
  sample_qc,
  file.path(
    out_dir,
    "Figure5A_sample_QC_by_meta_niche.csv"
  ),
  row.names = FALSE
)

all_de <- dplyr::bind_rows(
  de_results
)

utils::write.csv(
  all_de,
  file.path(
    out_dir,
    "Figure5A_all_limma_DE.csv"
  ),
  row.names = FALSE
)

all_gsea <- dplyr::bind_rows(
  gsea_results
)

# Make a CSV-safe copy.
all_gsea_csv <- all_gsea

if ("leadingEdge" %in% colnames(all_gsea_csv)) {
  all_gsea_csv$leadingEdge <- flatten_leading_edge(
    all_gsea_csv$leadingEdge
  )
}

utils::write.csv(
  all_gsea_csv,
  file.path(
    out_dir,
    "Figure5A_all_Hallmark_GSEA.csv"
  ),
  row.names = FALSE
)

saveRDS(
  de_results,
  file.path(
    out_dir,
    "Figure5A_limma_voom_DE_per_meta_niche.rds"
  )
)

saveRDS(
  gsea_results,
  file.path(
    out_dir,
    "Figure5A_Hallmark_GSEA_per_meta_niche.rds"
  )
)

# ------------------------------------------------------------------------------
# 6. Selected Figure 5A dot plot
# ------------------------------------------------------------------------------

plot_df <- all_gsea |>
  dplyr::filter(
    pathway %in% SELECTED_PATHWAYS
  ) |>
  dplyr::mutate(
    padj_plot = ifelse(
      is.na(padj),
      1,
      padj
    ),
    significance = -log10(
      pmax(
        padj_plot,
        .Machine$double.xmin
      )
    ),
    significance = pmin(
      significance,
      10
    ),
    meta_niche_label = factor(
      meta_niche_label,
      levels = unname(
        MN_LABELS[
          MN_DISPLAY_ORDER
        ]
      )
    ),
    pathway_label = pathway |>
      stringr::str_remove("^HALLMARK_") |>
      stringr::str_replace_all("_", " "),
    pathway_label = factor(
      pathway_label,
      levels = rev(
        SELECTED_PATHWAYS |>
          stringr::str_remove("^HALLMARK_") |>
          stringr::str_replace_all("_", " ")
      )
    )
  )

missing_selected <- setdiff(
  SELECTED_PATHWAYS,
  unique(
    as.character(
      all_gsea$pathway
    )
  )
)

if (length(missing_selected) > 0L) {
  warning(
    "Selected Hallmark pathway(s) absent from fgsea results: ",
    paste(missing_selected, collapse = ", ")
  )
}

plot_csv <- plot_df

if ("leadingEdge" %in% colnames(plot_csv)) {
  plot_csv$leadingEdge <- flatten_leading_edge(
    plot_csv$leadingEdge
  )
}

utils::write.csv(
  plot_csv,
  file.path(
    out_dir,
    "Figure5A_selected_Hallmark_GSEA.csv"
  ),
  row.names = FALSE
)

p_fig5a <- ggplot2::ggplot(
  plot_df,
  ggplot2::aes(
    x = meta_niche_label,
    y = pathway_label
  )
) +
  ggplot2::geom_point(
    ggplot2::aes(
      size = significance,
      fill = NES
    ),
    shape = 21,
    color = "black",
    stroke = 0.3
  ) +
  ggplot2::scale_fill_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0,
    name = "NES"
  ) +
  ggplot2::scale_size_continuous(
    range = c(1, 7),
    name = expression(
      -log[10]("adjusted p-value")
    )
  ) +
  ggplot2::labs(
    x = "Meta-niche",
    y = "Hallmark pathway"
  ) +
  ggplot2::theme_classic(
    base_size = 12
  ) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(
      angle = 45,
      hjust = 1,
      vjust = 1
    ),
    axis.text.y = ggplot2::element_text(
      size = 9
    ),
    legend.position = "right"
  )

ggplot2::ggsave(
  file.path(
    out_dir,
    "Figure5A_selected_Hallmark_dotplot.pdf"
  ),
  p_fig5a,
  width = 12,
  height = 6,
  useDingbats = FALSE
)

ggplot2::ggsave(
  file.path(
    out_dir,
    "Figure5A_selected_Hallmark_dotplot.png"
  ),
  p_fig5a,
  width = 12,
  height = 6,
  dpi = 300
)

# ------------------------------------------------------------------------------
# 7. Console summary
# ------------------------------------------------------------------------------

cat(
  "\n============================================================\n",
  "Figure 5A meta-niche Hallmark analysis complete\n",
  "============================================================\n",
  sep = ""
)

cat("\nMeta-niche QC:\n")
print(
  sample_qc,
  row.names = FALSE
)

sig_summary <- all_gsea |>
  dplyr::mutate(
    significant = !is.na(padj) & padj < 0.05,
    direction = dplyr::case_when(
      NES > 0 ~ "Mut_up",
      NES < 0 ~ "WT_up",
      TRUE ~ "zero"
    )
  ) |>
  dplyr::filter(
    significant
  ) |>
  dplyr::count(
    meta_niche,
    direction,
    name = "n_significant_Hallmarks"
  ) |>
  dplyr::arrange(
    factor(meta_niche, levels = MN_LEVELS),
    direction
  )

cat("\nSignificant Hallmark pathways (padj < 0.05):\n")
print(
  sig_summary,
  row.names = FALSE
)

cat(
  "\nPositive NES = NF1Mut enriched; negative NES = NF1WT enriched.\n"
)
