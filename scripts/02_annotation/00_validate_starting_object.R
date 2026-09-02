# Validate the reproducible analysis starting object.
source(file.path('config','project_paths.R'))
source(file.path('R','project_io.R'))
source(file.path('R','constants.R'))
suppressPackageStartupMessages(library(Seurat))

obj <- load_starting_object()
require_metadata(obj, c('patient_id','NF1','seurat_clusters'))
if (!'RNA' %in% Assays(obj)) stop('RNA assay is required.')
clusters <- sort(unique(as.character(obj$seurat_clusters)))
expected <- names(CLUSTER_ANNOTATION)
missing_expected <- setdiff(expected,clusters)
unexpected <- setdiff(clusters,expected)
if(length(missing_expected) || length(unexpected)) {
  stop('The starting object does not contain the final 0-25 cluster set used for manuscript annotation. ',
       'Missing: ',paste(missing_expected,collapse=', '),'; unexpected: ',paste(unexpected,collapse=', '),
       '. If this is an earlier object, run/validate 00_optional_rebuild_clusters.R before continuing.')
}
if (!'umap' %in% Reductions(obj)) warning('No UMAP reduction present. The publication UMAP cannot be redrawn exactly; see README.')
message('Cells: ', ncol(obj))
message('Patients/cores represented: ', length(unique(obj$patient_id)))
message('NF1 groups: ', paste(names(table(obj$NF1)), table(obj$NF1), sep='=', collapse=', '))
message('Validated final Seurat clusters: 0-25')
