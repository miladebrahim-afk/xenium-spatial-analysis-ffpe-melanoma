# Optional validation of the processed objects used as the analysis entry point.
source(file.path('config','project_paths.R')); source(file.path('R','project_io.R'))
suppressPackageStartupMessages(library(Seurat))
obj <- load_starting_object(); require_metadata(obj,c('patient_id','NF1','seurat_clusters'))
xlist <- load_xenium_list_input(); stopifnot(length(xlist)>0)
message('Processed starting inputs found. No re-merging is required for the downstream pipeline.')
