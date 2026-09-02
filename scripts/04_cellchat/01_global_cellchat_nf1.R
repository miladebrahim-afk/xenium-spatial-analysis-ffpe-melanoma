# ==============================================================================
# 01_global_cellchat_nf1.R
#
# Global NF1Mut vs NF1WT CellChat model fitting for Figure S3.
#
# Historical workflow:
#   - all retained Xenium cells, split by NF1 genotype
#   - RNA assay
#   - NormalizeData()
#   - group.by = "celltype"
#   - CellChatDB.human
#   - filterCommunication(min.cells = 10)
#
# This cleaned version fits ONE genotype at a time, saves immediately, and
# removes large intermediate objects before running the second genotype.
#
# Inputs
# ------
# data/processed/merged_obj_annotated_with_meta_niche_published_reference.rds
#
# Outputs
# -------
# results/cellchat/global_models/
#   CellChat_global_Mut.rds
#   CellChat_global_WT.rds
#   global_Mut_interactions.csv
#   global_WT_interactions.csv
#   global_cellchat_model_summary.csv
#   global_cellchat_cell_counts.csv
#   global_cellchat_sessionInfo.txt
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
  "global_models"
)

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

ASSAY_USE <- "RNA"
GROUP_COL <- "celltype"
MIN_CELLS <- 10L

GENOTYPES_TO_RUN <- c(
  "Mut",
  "WT"
)

# FALSE is safer: an existing successfully completed model will not be replaced.
OVERWRITE_EXISTING <- FALSE

# ------------------------------------------------------------------------------
# 1. Input / package preflight
# ------------------------------------------------------------------------------

cat("Global CellChat NF1 workflow\n")
cat("===========================\n")
cat("R version:        ", R.version.string, "\n", sep = "")
cat("Seurat version:   ", as.character(packageVersion("Seurat")), "\n", sep = "")
cat("CellChat version: ", as.character(packageVersion("CellChat")), "\n", sep = "")

if (!file.exists(merged_file)) {
  stop(
    "Missing annotated merged object: ",
    merged_file
  )
}

cat("\nLoading annotated merged object...\n")

merged_obj <- readRDS(
  merged_file
)

required_meta <- c(
  "NF1",
  GROUP_COL
)

missing_meta <- setdiff(
  required_meta,
  colnames(merged_obj@meta.data)
)

if (length(missing_meta) > 0L) {
  stop(
    "Missing required metadata: ",
    paste(missing_meta, collapse = ", ")
  )
}

if (!ASSAY_USE %in% Assays(merged_obj)) {
  stop(
    "RNA assay is missing. Available assays: ",
    paste(Assays(merged_obj), collapse = ", ")
  )
}

cat(
  "Merged cells: ",
  format(ncol(merged_obj), big.mark = ","),
  "\n",
  sep = ""
)

# Normalize NF1 spelling in memory only.
nf1_raw <- as.character(
  merged_obj$NF1
)

merged_obj$NF1_cellchat <- dplyr::case_when(
  nf1_raw %in% c(
    "Mut",
    "NF1Mut",
    "NF1_Mut",
    "NF1 Mut"
  ) ~ "Mut",
  nf1_raw %in% c(
    "WT",
    "NF1WT",
    "NF1_WT",
    "NF1 WT"
  ) ~ "WT",
  TRUE ~ NA_character_
)

if (anyNA(merged_obj$NF1_cellchat)) {
  stop(
    "Unexpected NF1 value(s): ",
    paste(
      sort(
        unique(
          nf1_raw[
            is.na(
              merged_obj$NF1_cellchat
            )
          ]
        )
      ),
      collapse = ", "
    )
  )
}

# ------------------------------------------------------------------------------
# 2. Cell-count QC
# ------------------------------------------------------------------------------

global_cell_counts <- merged_obj@meta.data %>%
  filter(
    !is.na(.data[[GROUP_COL]])
  ) %>%
  mutate(
    NF1_cellchat = merged_obj$NF1_cellchat[
      rownames(.)
    ]
  ) %>%
  count(
    NF1 = NF1_cellchat,
    celltype = .data[[GROUP_COL]],
    name = "n_cells"
  ) %>%
  arrange(
    factor(
      NF1,
      levels = GENOTYPES_TO_RUN
    ),
    desc(n_cells)
  )

# The rowname-based mutate above can be awkward on some dplyr versions.
# Rebuild robustly if necessary.
if (
  nrow(global_cell_counts) == 0L ||
  anyNA(global_cell_counts$NF1)
) {
  tmp_md <- merged_obj@meta.data
  tmp_md$NF1_cellchat <- merged_obj$NF1_cellchat

  global_cell_counts <- tmp_md %>%
    filter(
      !is.na(.data[[GROUP_COL]])
    ) %>%
    count(
      NF1 = NF1_cellchat,
      celltype = .data[[GROUP_COL]],
      name = "n_cells"
    ) %>%
    arrange(
      factor(
        NF1,
        levels = GENOTYPES_TO_RUN
      ),
      desc(n_cells)
    )

  rm(tmp_md)
}

write.csv(
  global_cell_counts,
  file.path(
    out_dir,
    "global_cellchat_cell_counts.csv"
  ),
  row.names = FALSE
)

global_preflight <- global_cell_counts %>%
  group_by(NF1) %>%
  summarise(
    total_cells = sum(n_cells),
    n_celltypes = n(),
    n_celltypes_ge10 = sum(
      n_cells >= MIN_CELLS
    ),
    .groups = "drop"
  )

cat("\nGlobal CellChat preflight:\n")
print(
  global_preflight,
  row.names = FALSE
)

if (
  !all(
    GENOTYPES_TO_RUN %in%
      global_preflight$NF1
  )
) {
  stop(
    "One or both NF1 genotype groups are missing."
  )
}

if (
  any(
    global_preflight$n_celltypes_ge10 < 2L
  )
) {
  stop(
    "Too few cell types with >=10 cells in one genotype."
  )
}

# ------------------------------------------------------------------------------
# 3. Global CellChat helper
# ------------------------------------------------------------------------------

run_global_cellchat <- function(
    seurat_obj,
    genotype,
    assay_use = ASSAY_USE,
    group_col = GROUP_COL,
    min_cells = MIN_CELLS
) {

  cat("\n============================================================\n")
  cat("GLOBAL CELLCHAT | NF1 ", genotype, "\n", sep = "")
  cat("============================================================\n")

  md <- seurat_obj@meta.data

  cells_use <- rownames(md)[
    as.character(
      seurat_obj$NF1_cellchat
    ) == genotype &
      !is.na(md[[group_col]])
  ]

  cat(
    "Cells selected: ",
    format(
      length(cells_use),
      big.mark = ","
    ),
    "\n",
    sep = ""
  )

  if (
    length(cells_use) <
      min_cells
  ) {
    stop(
      "Too few cells for NF1 ",
      genotype
    )
  }

  # Historical global workflow split the merged Seurat object by genotype.
  seu <- subset(
    seurat_obj,
    cells = cells_use
  )

  # Match historical global CellChat workflow:
  # DefaultAssay <- RNA followed by NormalizeData().
  DefaultAssay(seu) <- assay_use

  cat(
    "NormalizeData(RNA)...\n"
  )

  seu <- NormalizeData(
    seu,
    assay = assay_use,
    verbose = FALSE
  )

  # Drop unused identity levels for clean CellChat groups.
  seu[[group_col]][, 1] <- droplevels(
    factor(
      seu[[group_col]][, 1]
    )
  )

  cat(
    "Cell types: ",
    length(
      levels(
        seu[[group_col]][, 1]
      )
    ),
    "\n",
    sep = ""
  )

  cat("createCellChat...\n")

  cellchat <- createCellChat(
    object = seu,
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

  cat(
    "filterCommunication(min.cells = ",
    min_cells,
    ")...\n",
    sep = ""
  )

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
    NF1 = genotype,
    n_cells = ncol(seu),
    n_celltypes = length(
      levels(
        cellchat@idents
      )
    ),
    total_interactions = sum(
      cellchat@net$count
    ),
    total_interaction_weight = sum(
      cellchat@net$weight
    ),
    n_signaling_pathways = length(
      cellchat@netP$pathways
    ),
    stringsAsFactors = FALSE
  )

  attr(
    cellchat,
    "github_run_metadata"
  ) <- run_metadata

  # The Seurat subset is the large object; remove it before returning.
  rm(seu)
  gc(
    verbose = FALSE
  )

  list(
    object = cellchat,
    metadata = run_metadata
  )
}

# ------------------------------------------------------------------------------
# 4. Fit one genotype at a time and save immediately
# ------------------------------------------------------------------------------

run_summary <- list()

for (gt in GENOTYPES_TO_RUN) {

  model_file <- file.path(
    out_dir,
    paste0(
      "CellChat_global_",
      gt,
      ".rds"
    )
  )

  interaction_file <- file.path(
    out_dir,
    paste0(
      "global_",
      gt,
      "_interactions.csv"
    )
  )

  if (
    file.exists(model_file) &&
      !OVERWRITE_EXISTING
  ) {

    cat(
      "\nExisting model found; skipping: ",
      model_file,
      "\n",
      sep = ""
    )

    existing <- readRDS(
      model_file
    )

    existing <- updateCellChat(
      existing
    )

    md_existing <- attr(
      existing,
      "github_run_metadata"
    )

    if (
      is.null(md_existing)
    ) {
      md_existing <- data.frame(
        NF1 = gt,
        n_cells = NA_integer_,
        n_celltypes = length(
          levels(
            existing@idents
          )
        ),
        total_interactions = sum(
          existing@net$count
        ),
        total_interaction_weight = sum(
          existing@net$weight
        ),
        n_signaling_pathways = length(
          existing@netP$pathways
        ),
        stringsAsFactors = FALSE
      )
    }

    run_summary[[gt]] <- md_existing

    rm(existing)
    gc(
      verbose = FALSE
    )

    next
  }

  result <- run_global_cellchat(
    seurat_obj = merged_obj,
    genotype = gt
  )

  saveRDS(
    result$object,
    model_file,
    compress = FALSE
  )

  write.csv(
    subsetCommunication(
      result$object
    ),
    interaction_file,
    row.names = FALSE
  )

  run_summary[[gt]] <- result$metadata

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
    " | pathways: ",
    result$metadata$n_signaling_pathways,
    "\n",
    sep = ""
  )

  rm(result)
  gc(
    verbose = FALSE
  )
}

# ------------------------------------------------------------------------------
# 5. Save run summary / session information
# ------------------------------------------------------------------------------

run_summary_df <- bind_rows(
  run_summary
)

write.csv(
  run_summary_df,
  file.path(
    out_dir,
    "global_cellchat_model_summary.csv"
  ),
  row.names = FALSE
)

capture.output(
  sessionInfo(),
  file = file.path(
    out_dir,
    "global_cellchat_sessionInfo.txt"
  )
)

cat("\nGlobal CellChat model summary:\n")
print(
  run_summary_df,
  row.names = FALSE
)

cat(
  "\nDONE: global NF1 CellChat models complete.\n"
)

cat(
  "Next step: generate Figure S3 differential networks and signaling-pathway panels.\n"
)
