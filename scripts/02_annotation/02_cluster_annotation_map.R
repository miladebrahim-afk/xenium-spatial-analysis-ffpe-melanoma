# Document the manual biological annotation decision.
# The numeric clusters were assigned final labels after review of FindAllMarkers
# results and canonical marker expression. This is intentionally explicit rather
# than automated: the biological annotation was a manual interpretation step.
source(file.path('config','project_paths.R'))
source(file.path('R','constants.R'))

out <- file.path(RESULTS_DIR,'annotation')
dir.create(out, recursive=TRUE, showWarnings=FALSE)
map <- data.frame(seurat_cluster=names(CLUSTER_ANNOTATION), final_celltype=unname(CLUSTER_ANNOTATION), stringsAsFactors=FALSE)
write.csv(map, file.path(out,'cluster_to_final_celltype.csv'), row.names=FALSE)
print(map)
