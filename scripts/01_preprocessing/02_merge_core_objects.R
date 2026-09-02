# OPTIONAL UPSTREAM PROVENANCE — not required when starting from merged_obj_filtered.rds
# Merge individual core Seurat objects
# ------------------------------------
# Input:  data/per_core/*.rds
# Output: data/processed/xenium_list.rds
#         data/processed/merged_obj_unfiltered.rds
#
# The surviving analysis loaded individual core RDS objects, prefixed cell
# barcodes with patient/core information, stripped images/tools from a copy for
# memory-efficient merging, and merged the objects stepwise. The original
# spatial objects are retained here as xenium_list for spatial analyses.

source(file.path('config', 'project_paths.R'))

suppressPackageStartupMessages(library(Seurat))

core_files <- list.files(PER_CORE_DIR, pattern = '\\.rds$', full.names = TRUE)
if (!length(core_files)) {
  stop('No per-core RDS files found in: ', PER_CORE_DIR)
}

xenium_list <- lapply(core_files, readRDS)
names(xenium_list) <- tools::file_path_sans_ext(basename(core_files))

for (i in seq_along(xenium_list)) {
  obj <- xenium_list[[i]]
  if (!'patient_id' %in% colnames(obj@meta.data)) {
    stop('Missing patient_id in ', basename(core_files[i]))
  }
  if (!'NF1' %in% colnames(obj@meta.data)) {
    stop('Missing NF1 in ', basename(core_files[i]))
  }

  pid <- as.character(unique(obj$patient_id))
  if (length(pid) != 1L || is.na(pid)) {
    stop('Expected one patient_id per core in ', basename(core_files[i]))
  }

  # Prefixing protects barcode identity across cores while keeping the original
  # per-core object spatially usable.
  if (!all(startsWith(colnames(obj), paste0('Patient', pid, '_')))) {
    obj <- RenameCells(obj, add.cell.id = paste0('Patient', pid))
  }
  xenium_list[[i]] <- obj
}

saveRDS(xenium_list, XENIUM_LIST_FILE)

# Create a merge copy without spatial images/tools to control memory. The
# spatially intact objects remain available in xenium_list.rds.
merge_list <- lapply(xenium_list, function(obj) {
  obj@images <- list()
  obj@tools <- list()
  obj
})

merged_obj <- merge_list[[1]]
if (length(merge_list) > 1L) {
  for (i in 2:length(merge_list)) {
    merged_obj <- merge(
      x = merged_obj,
      y = merge_list[[i]],
      merge.data = TRUE
    )
    gc(verbose = FALSE)
  }
}

saveRDS(merged_obj, MERGED_UNFILTERED_FILE)
message('Saved ', length(xenium_list), ' per-core objects and merged object.')
