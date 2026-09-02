# Check the two processed inputs and software packages before running the pipeline.
source(file.path('config','project_paths.R'))
source(file.path('R','project_io.R'))

required <- c('Seurat','SeuratObject','dplyr','tidyr','tibble','ggplot2','Matrix','dbscan','FNN',
              'DESeq2','edgeR','limma','fgsea','msigdbr','ggpubr','pheatmap',
              'patchwork','openxlsx','ggrepel','stringr','purrr','scales',
              'clusterProfiler','org.Hs.eg.db','gprofiler2','readxl','future')
optional <- c('CellChat','ComplexHeatmap','circlize','arrow','data.table','harmony','renv')

missing_required <- required[!vapply(required,requireNamespace,logical(1),quietly=TRUE)]
missing_optional <- optional[!vapply(optional,requireNamespace,logical(1),quietly=TRUE)]
if (length(missing_required)) stop('Missing required R packages: ',paste(missing_required,collapse=', '))
if (length(missing_optional)) message('Optional packages not installed: ',paste(missing_optional,collapse=', '))
if (!file.exists(STARTING_OBJECT_FILE)) stop('Missing: ',STARTING_OBJECT_FILE)
if (!file.exists(XENIUM_LIST_INPUT_FILE)) stop('Missing: ',XENIUM_LIST_INPUT_FILE)

obj <- readRDS(STARTING_OBJECT_FILE)
require_metadata(obj,c('patient_id','NF1','seurat_clusters'))
if (!'RNA' %in% Seurat::Assays(obj)) stop('RNA assay is required in merged_obj_filtered.rds.')
x <- readRDS(XENIUM_LIST_INPUT_FILE)
if (!is.list(x) || !length(x)) stop('xenium_list.rds must contain a non-empty list of per-core Seurat objects.')
message('Input check passed: ',ncol(obj),' merged cells; ',length(x),' per-core objects.')
