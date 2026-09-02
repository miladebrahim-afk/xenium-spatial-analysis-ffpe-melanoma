# OPTIONAL: rebuild the final Seurat clustering/UMAP from merged_obj_filtered.rds.
#
# This mirrors the final surviving re-annotation block:
#   SCT assay -> PCA -> neighbors (dims 1:30) -> FindClusters(resolution=0.5) -> UMAP.
#
# The main pipeline does NOT run this automatically because the manuscript's saved
# merged_obj_filtered object already contains the final seurat_clusters and UMAP.
# Preserving those fields is the safest way to redraw the published Figure 1 exactly.
# Run this script only if you intentionally want to rebuild clustering from the
# processed expression object and then use the rebuilt object as the starting file.

source(file.path('config', 'project_paths.R'))
source(file.path('R', 'project_io.R'))

suppressPackageStartupMessages(library(Seurat))

obj <- load_starting_object()

if (!'SCT' %in% Assays(obj)) {
  stop(
    'The surviving final re-annotation block used the SCT assay. ',
    'This starting object does not contain an SCT assay, so the documented ',
    'clustering step cannot be reconstructed faithfully from this file.'
  )
}

DefaultAssay(obj) <- 'SCT'
obj <- RunPCA(obj, verbose = FALSE)
obj <- FindNeighbors(obj, dims = seq_len(N_PCS))
obj <- FindClusters(obj, resolution = FINAL_CLUSTER_RESOLUTION)
obj <- RunUMAP(obj, dims = seq_len(N_PCS), verbose = FALSE)

message('Clusters after rebuild:')
print(table(obj$seurat_clusters))

reclustered_file <- file.path(PROCESSED_DIR, 'merged_obj_filtered_reclustered.rds')
saveRDS(obj, reclustered_file)
message('Optional rebuilt object saved to: ', reclustered_file)
message(
  'To use it for the remainder of the pipeline, either rename/copy it to ',
  'data/processed/merged_obj_filtered.rds or update STARTING_OBJECT_FILE.'
)
