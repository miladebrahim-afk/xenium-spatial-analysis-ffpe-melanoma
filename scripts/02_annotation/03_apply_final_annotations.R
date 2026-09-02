# Apply the final 0-25 cluster-to-cell-type map and save the annotated object.
source(file.path('config','project_paths.R'))
source(file.path('R','project_io.R'))
source(file.path('R','constants.R'))
suppressPackageStartupMessages({library(Seurat); library(dplyr)})

obj <- load_starting_object()
require_metadata(obj, c('patient_id','NF1','seurat_clusters'))
clusters <- as.character(obj$seurat_clusters)
unmapped <- setdiff(unique(clusters), names(CLUSTER_ANNOTATION))
if (length(unmapped)) stop('Unmapped Seurat clusters: ', paste(unmapped, collapse=', '))

obj$celltype <- unname(CLUSTER_ANNOTATION[clusters])
obj$celltype <- factor(obj$celltype, levels=CELLTYPE_ORDER)
obj$major_celltype <- unname(BROAD_CELLTYPE_MAP[as.character(obj$celltype)])
obj$major_celltype <- factor(obj$major_celltype, levels=BROAD_CELLTYPE_ORDER)
Idents(obj) <- 'celltype'

if (anyNA(obj$celltype)) stop('NA celltype labels created during annotation.')
saveRDS(obj, ANNOTATED_OBJECT_FILE)
message('Annotated object saved: ', ANNOTATED_OBJECT_FILE)
