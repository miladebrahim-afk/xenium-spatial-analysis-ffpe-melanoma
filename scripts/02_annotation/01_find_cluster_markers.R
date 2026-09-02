# Figure 1 annotation provenance: marker discovery BEFORE manual biological labels.
# Source: final block of re_annotation_filtering.R.
source(file.path('config','project_paths.R'))
source(file.path('R','project_io.R'))
source(file.path('R','constants.R'))
source(file.path('R','analysis_helpers.R'))

suppressPackageStartupMessages({library(Seurat);library(dplyr);library(ggplot2)})
obj <- load_starting_object()
require_metadata(obj,c('patient_id','NF1','seurat_clusters'))
Idents(obj) <- 'seurat_clusters'
out <- ensure_dir(file.path(RESULTS_DIR,'annotation'))
figout <- ensure_dir(file.path(FIGURES_DIR,'annotation'))

# Memory-safe execution for the million-cell marker step. This changes only
# parallelization/global-size handling, not the marker test or its parameters.
.old_future_plan <- future::plan()
.old_future_max <- getOption('future.globals.maxSize')
future::plan(future::sequential)
options(future.globals.maxSize = 150 * 1024^3)

# The retained final notebook used the SCT assay and PrepSCTFindMarkers before
# FindAllMarkers. If SCT is absent, fall back transparently to normalized RNA.
if ('SCT' %in% Assays(obj)) {
  DefaultAssay(obj) <- 'SCT'
  obj <- PrepSCTFindMarkers(obj,verbose=FALSE)
} else {
  obj <- normalize_rna_if_needed(obj)
}

cluster_markers <- FindAllMarkers(
  object=obj,
  only.pos=TRUE,
  min.pct=0,
  logfc.threshold=0
)
write.csv(cluster_markers,file.path(out,'cluster_markers_all.csv'),row.names=FALSE)

top20 <- cluster_markers |>
  group_by(cluster) |>
  slice_max(order_by=avg_log2FC,n=20,with_ties=FALSE) |>
  ungroup()
write.csv(top20,file.path(out,'cluster_markers_top20.csv'),row.names=FALSE)

# Marker-evidence DotPlot while identities are still numeric clusters.
obj <- normalize_rna_if_needed(obj)
marker_vec <- intersect(unique(unlist(MARKERS_BY_BIOLOGY,use.names=FALSE)),rownames(obj[['RNA']]))
p <- DotPlot(obj,features=marker_vec,group.by='seurat_clusters',assay='RNA',dot.scale=5) +
  RotatedAxis() + ggtitle('Canonical markers by numeric Seurat cluster')
save_gg(p,file.path(figout,'cluster_canonical_marker_dotplot.pdf'),18,9)

# Restore caller settings.
future::plan(.old_future_plan)
options(future.globals.maxSize = .old_future_max)
