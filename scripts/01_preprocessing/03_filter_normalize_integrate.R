# OPTIONAL UPSTREAM PROVENANCE — not required when starting from merged_obj_filtered.rds
# Filter, normalize, and batch-correct the merged Xenium object
# -------------------------------------------------------------
# Input:  data/processed/merged_obj_unfiltered.rds
# Output: data/processed/merged_obj_integrated.rds
#
# This script reconstructs the manuscript-reported preprocessing sequence using
# parameters that are explicitly reported in the manuscript where available.
# It should be validated once against the original final Seurat object before
# the repository is released publicly.

source(file.path('config', 'project_paths.R'))

suppressPackageStartupMessages({
  library(Seurat)
  library(harmony)
})

if (!file.exists(MERGED_UNFILTERED_FILE)) {
  stop('Run scripts/01_preprocessing/02_merge_core_objects.R first.')
}

merged_obj <- readRDS(MERGED_UNFILTERED_FILE)

# Reported QC threshold in the final Results/manuscript.
if (!'nCount_RNA' %in% colnames(merged_obj@meta.data)) {
  stop('nCount_RNA metadata is required for the <50 transcript filter.')
}
merged_obj <- subset(merged_obj, subset = nCount_RNA >= MIN_TRANSCRIPTS_PER_CELL)

# The final Results report exclusion of cores with <5,000 cells.
core_field <- CORE_ID_FIELD
if (!core_field %in% colnames(merged_obj@meta.data)) {
  stop('Core identifier field not found: ', core_field)
}
core_counts <- table(merged_obj@meta.data[[core_field]])
keep_cores <- names(core_counts[core_counts >= MIN_CELLS_PER_CORE])
merged_obj <- subset(
  merged_obj,
  cells = rownames(merged_obj@meta.data)[merged_obj@meta.data[[core_field]] %in% keep_cores]
)

DefaultAssay(merged_obj) <- 'RNA'
merged_obj <- SCTransform(
  merged_obj,
  assay = 'RNA',
  variable.features.n = SCT_VARIABLE_FEATURES,
  return.only.var.genes = FALSE,
  verbose = FALSE
)

merged_obj <- RunPCA(merged_obj, assay = 'SCT', npcs = N_PCS, verbose = FALSE)

pca_embeddings <- Embeddings(merged_obj, 'pca')[, seq_len(N_PCS), drop = FALSE]
harmony_embeddings <- HarmonyMatrix(
  data_mat = pca_embeddings,
  meta_data = merged_obj@meta.data,
  vars_use = HARMONY_GROUP_FIELD,
  do_pca = FALSE,
  verbose = TRUE
)

merged_obj[['harmony']] <- CreateDimReducObject(
  embeddings = harmony_embeddings,
  key = 'harmony_',
  assay = 'SCT'
)

merged_obj <- FindNeighbors(
  merged_obj,
  reduction = 'harmony',
  dims = seq_len(N_PCS),
  verbose = FALSE
)
merged_obj <- FindClusters(
  merged_obj,
  resolution = CLUSTER_RESOLUTION,
  verbose = FALSE
)
merged_obj <- RunUMAP(
  merged_obj,
  reduction = 'harmony',
  dims = seq_len(N_PCS),
  verbose = FALSE
)

saveRDS(merged_obj, MERGED_INTEGRATED_FILE)
message('Integrated object saved: ', MERGED_INTEGRATED_FILE)
