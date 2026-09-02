# ==============================================================================
# 01_run_niche_cellchat_models.R
#
# Run CellChat models for the meta-niches used in the final manuscript figures.
#
# Primary manuscript CellChat niches:
#   MN2 = NC+M2+CAF
#   MN3 = IFN-myeloid+T cell
#   MN6 = iCAF+M2
#
# The archived revision script also explored MN4, but MN4 is not required for
# the final Figure 2 CellChat panels. Set INCLUDE_MN4_EXPLORATORY <- TRUE if
# that exploratory model is desired.
#
# IMPORTANT
# ---------
# - Uses the publication-reference merged object.
# - Runs one MN x NF1 model at a time and saves it immediately.
# - Does NOT keep all CellChat models in memory simultaneously.
# - Uses SCT assay and CellChatDB.human, matching the archived workflow.
# - Cell types with <10 cells in a specific MN x NF1 subset are removed before
#   CellChat, matching the archived min.cells = 10 workflow.
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(CellChat)
  library(dplyr)
})

# ------------------------------------------------------------------------------
# 0. Configuration
# ------------------------------------------------------------------------------

merged_file <- file.path(
  "data",
  "processed",
  "merged_obj_annotated_with_meta_niche_published_reference.rds"
)

out_dir <- file.path(
  "results",
  "cellchat",
  "niche_models"
)

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

ASSAY_USE <- "SCT"
GROUP_COL <- "celltype"
MIN_CELLS <- 10L

PRIMARY_MNS <- c(
  "MetaNiche_2",
  "MetaNiche_3",
  "MetaNiche_6"
)

INCLUDE_MN4_EXPLORATORY <- FALSE

META_NICHES_TO_RUN <- if (INCLUDE_MN4_EXPLORATORY) {
  c(
    "MetaNiche_2",
    "MetaNiche_3",
    "MetaNiche_4",
    "MetaNiche_6"
  )
} else {
  PRIMARY_MNS
}

GENOTYPES_TO_RUN <- c(
  "Mut",
  "WT"
)

# ------------------------------------------------------------------------------
# 1. Package/input preflight
# ------------------------------------------------------------------------------

cat("CellChat niche-model workflow\n")
cat("============================\n")
cat("R version:       ", R.version.string, "\n", sep = "")
cat("Seurat version:  ", as.character(packageVersion("Seurat")), "\n", sep = "")
cat("CellChat version:", as.character(packageVersion("CellChat")), "\n", sep = "")

if (!file.exists(merged_file)) {
  stop(
    "Missing publication-reference merged object: ",
    merged_file,
    "\nRun scripts/03_meta_niches/03_apply_published_meta_niche_reference.R first."
  )
}

cat("\nLoading publication-reference merged object...\n")
merged_obj <- readRDS(
  merged_file
)

required_meta <- c(
  "meta_niche",
  "NF1",
  GROUP_COL
)

missing_meta <- setdiff(
  required_meta,
  colnames(merged_obj@meta.data)
)

if (length(missing_meta) > 0L) {
  stop(
    "Merged object is missing required metadata: ",
    paste(missing_meta, collapse = ", ")
  )
}

if (!ASSAY_USE %in% Assays(merged_obj)) {
  stop(
    "Assay '",
    ASSAY_USE,
    "' is missing. Available assays: ",
    paste(Assays(merged_obj), collapse = ", ")
  )
}

cat(
  "Merged cells: ",
  ncol(merged_obj),
  "\n",
  sep = ""
)

# Normalize genotype spelling without changing the saved RDS.
nf1_raw <- as.character(
  merged_obj$NF1
)

merged_obj$NF1_cellchat <- dplyr::case_when(
  nf1_raw %in% c("Mut", "NF1Mut", "NF1_Mut", "NF1 Mut") ~ "Mut",
  nf1_raw %in% c("WT", "NF1WT", "NF1_WT", "NF1 WT") ~ "WT",
  TRUE ~ NA_character_
)

if (anyNA(merged_obj$NF1_cellchat)) {
  stop(
    "Unexpected NF1 metadata value(s): ",
    paste(
      sort(
        unique(
          nf1_raw[
            is.na(merged_obj$NF1_cellchat)
          ]
        )
      ),
      collapse = ", "
    )
  )
}

# ------------------------------------------------------------------------------
# 2. Preflight counts
# ------------------------------------------------------------------------------

cat("\nCell counts before CellChat filtering:\n")

cellchat_cell_counts <- merged_obj@meta.data %>%
  mutate(
    NF1_cellchat = merged_obj$NF1_cellchat
  ) %>%
  filter(
    meta_niche %in% META_NICHES_TO_RUN,
    !is.na(.data[[GROUP_COL]])
  ) %>%
  count(
    meta_niche,
    NF1 = NF1_cellchat,
    celltype = .data[[GROUP_COL]],
    name = "n_cells"
  ) %>%
  arrange(
    factor(
      meta_niche,
      levels = META_NICHES_TO_RUN
    ),
    factor(
      NF1,
      levels = GENOTYPES_TO_RUN
    ),
    desc(n_cells)
  )

write.csv(
  cellchat_cell_counts,
  file.path(
    out_dir,
    "cellchat_niche_cell_counts_before_min10_filter.csv"
  ),
  row.names = FALSE
)

run_preflight <- cellchat_cell_counts %>%
  mutate(
    retained_celltype = n_cells >= MIN_CELLS
  ) %>%
  group_by(
    meta_niche,
    NF1
  ) %>%
  summarise(
    total_cells = sum(n_cells),
    retained_cells = sum(
      n_cells[retained_celltype]
    ),
    n_celltypes_total = n(),
    n_celltypes_retained = sum(
      retained_celltype
    ),
    .groups = "drop"
  ) %>%
  arrange(
    factor(
      meta_niche,
      levels = META_NICHES_TO_RUN
    ),
    factor(
      NF1,
      levels = GENOTYPES_TO_RUN
    )
  )

cat("\nCellChat run preflight:\n")
print(
  run_preflight,
  row.names = FALSE
)

write.csv(
  run_preflight,
  file.path(
    out_dir,
    "cellchat_niche_run_preflight.csv"
  ),
  row.names = FALSE
)

expected_runs <- expand.grid(
  meta_niche = META_NICHES_TO_RUN,
  NF1 = GENOTYPES_TO_RUN,
  stringsAsFactors = FALSE
)

run_check <- expected_runs %>%
  left_join(
    run_preflight,
    by = c(
      "meta_niche",
      "NF1"
    )
  )

if (anyNA(run_check$retained_cells)) {
  bad <- run_check %>%
    filter(
      is.na(retained_cells)
    )

  stop(
    "No cells found for one or more requested MN x NF1 runs:\n",
    paste(
      paste0(
        bad$meta_niche,
        " ",
        bad$NF1
      ),
      collapse = "\n"
    )
  )
}

if (any(run_check$n_celltypes_retained < 2L)) {
  bad <- run_check %>%
    filter(
      n_celltypes_retained < 2L
    )

  stop(
    "Fewer than two cell types with >= ",
    MIN_CELLS,
    " cells in:\n",
    paste(
      paste0(
        bad$meta_niche,
        " ",
        bad$NF1
      ),
      collapse = "\n"
    )
  )
}

# ------------------------------------------------------------------------------
# 3. CellChat helper
# ------------------------------------------------------------------------------

run_cellchat_one <- function(
    seurat_obj,
    mn,
    genotype,
    assay_use = ASSAY_USE,
    group_col = GROUP_COL,
    min_cells = MIN_CELLS
) {

  cat("\n============================================================\n")
  cat("Running CellChat: ", mn, " | NF1 ", genotype, "\n", sep = "")
  cat("============================================================\n")

  md <- seurat_obj@meta.data

  cells_use <- rownames(md)[
    as.character(md$meta_niche) == mn &
      as.character(md$NF1_cellchat) == genotype &
      !is.na(md[[group_col]])
  ]

  cat(
    "Cells before rare-celltype filter: ",
    length(cells_use),
    "\n",
    sep = ""
  )

  if (length(cells_use) < min_cells) {
    stop(
      "Too few cells for ",
      mn,
      " ",
      genotype,
      "."
    )
  }

  celltype_counts <- table(
    droplevels(
      factor(
        md[cells_use, group_col]
      )
    )
  )

  keep_celltypes <- names(
    celltype_counts[
      celltype_counts >= min_cells
    ]
  )

  keep_cells <- cells_use[
    as.character(
      md[cells_use, group_col]
    ) %in% keep_celltypes
  ]

  cat(
    "Cells after rare-celltype filter:  ",
    length(keep_cells),
    "\n",
    sep = ""
  )

  cat(
    "Cell types retained (",
    length(keep_celltypes),
    "): ",
    paste(
      keep_celltypes,
      collapse = ", "
    ),
    "\n",
    sep = ""
  )

  if (length(keep_celltypes) < 2L) {
    stop(
      "Fewer than two cell types remain after min.cells filter."
    )
  }

  # Subset only after deciding the exact cells to retain.
  seurat_sub <- subset(
    seurat_obj,
    cells = keep_cells
  )

  seurat_sub[[group_col]][, 1] <- droplevels(
    factor(
      seurat_sub[[group_col]][, 1]
    )
  )

  # Match the archived workflow:
  # createCellChat(Seurat, assay = "SCT", group.by = "celltype")
  cellchat <- createCellChat(
    object = seurat_sub,
    assay = assay_use,
    group.by = group_col
  )

  cellchat@DB <- CellChatDB.human

  cat("subsetData...\n")
  cellchat <- subsetData(
    cellchat
  )

  cat("identifyOverExpressedGenes...\n")
  cellchat <- identifyOverExpressedGenes(
    cellchat
  )

  cat("identifyOverExpressedInteractions...\n")
  cellchat <- identifyOverExpressedInteractions(
    cellchat
  )

  cat("computeCommunProb...\n")
  cellchat <- computeCommunProb(
    cellchat
  )

  cat("filterCommunication(min.cells = ", min_cells, ")...\n", sep = "")
  cellchat <- filterCommunication(
    cellchat,
    min.cells = min_cells
  )

  cat("computeCommunProbPathway...\n")
  cellchat <- computeCommunProbPathway(
    cellchat
  )

  cat("aggregateNet...\n")
  cellchat <- aggregateNet(
    cellchat
  )

  cat("netAnalysis_computeCentrality...\n")
  cellchat <- netAnalysis_computeCentrality(
    cellchat,
    slot.name = "netP"
  )

  run_metadata <- data.frame(
    meta_niche = mn,
    NF1 = genotype,
    n_cells = length(keep_cells),
    n_celltypes = length(keep_celltypes),
    total_interactions = sum(
      cellchat@net$count
    ),
    total_interaction_weight = sum(
      cellchat@net$weight
    ),
    stringsAsFactors = FALSE
  )

  attr(
    cellchat,
    "github_run_metadata"
  ) <- run_metadata

  rm(
    seurat_sub
  )
  gc(
    verbose = FALSE
  )

  list(
    object = cellchat,
    metadata = run_metadata
  )
}

# ------------------------------------------------------------------------------
# 4. Run one model at a time and save immediately
# ------------------------------------------------------------------------------

cat("\nStarting CellChat model runs...\n")
cat(
  "Requested MNs: ",
  paste(
    META_NICHES_TO_RUN,
    collapse = ", "
  ),
  "\n",
  sep = ""
)

run_summary <- list()
run_index <- 0L

for (mn in META_NICHES_TO_RUN) {

  for (gt in GENOTYPES_TO_RUN) {

    run_index <- run_index + 1L

    result_name <- paste0(
      mn,
      "_",
      gt
    )

    model_file <- file.path(
      out_dir,
      paste0(
        "CellChat_",
        result_name,
        ".rds"
      )
    )

    cat(
      "\nModel ",
      run_index,
      "/",
      length(META_NICHES_TO_RUN) *
        length(GENOTYPES_TO_RUN),
      ": ",
      result_name,
      "\n",
      sep = ""
    )

    result <- run_cellchat_one(
      seurat_obj = merged_obj,
      mn = mn,
      genotype = gt
    )

    saveRDS(
      result$object,
      model_file,
      compress = FALSE
    )

    run_summary[[result_name]] <- result$metadata

    cat(
      "Saved: ",
      model_file,
      "\n",
      sep = ""
    )

    cat(
      "Interactions: ",
      result$metadata$total_interactions,
      " | total weight: ",
      signif(
        result$metadata$total_interaction_weight,
        6
      ),
      "\n",
      sep = ""
    )

    rm(
      result
    )
    gc(
      verbose = FALSE
    )
  }
}

# ------------------------------------------------------------------------------
# 5. Save summary / session information
# ------------------------------------------------------------------------------

run_summary_df <- bind_rows(
  run_summary
)

write.csv(
  run_summary_df,
  file.path(
    out_dir,
    "cellchat_niche_model_summary.csv"
  ),
  row.names = FALSE
)

session_file <- file.path(
  out_dir,
  "cellchat_sessionInfo.txt"
)

capture.output(
  sessionInfo(),
  file = session_file
)

cat("\nCellChat model summary:\n")
print(
  run_summary_df,
  row.names = FALSE
)

cat("\nSaved session information: ", session_file, "\n", sep = "")

cat(
  "\nDONE: niche-specific CellChat model fitting complete.\n"
)

cat(
  "Next step: run the comparison/figure script using the saved CellChat RDS files.\n"
)
