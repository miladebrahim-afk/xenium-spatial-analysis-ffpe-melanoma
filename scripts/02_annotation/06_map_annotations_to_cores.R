# Map regenerated integrated cell-type annotations back to the spatially intact
# per-core Xenium objects using a memory-safe two-stage workflow.
#
# IMPORTANT:
#   - The merged million-cell Seurat object and the 39-core Xenium list are
#     NEVER kept in memory at the same time.
#   - This script ONLY maps annotations and saves the annotated core list.
#     Spatial plotting is intentionally performed in a separate script.
#   - Matching first uses direct cell names when possible. Otherwise it uses
#     patient_id + Xenium barcode and, when needed, a shared core/sample field.
#     Ambiguous duplicate barcode matches cause an explicit error rather than
#     silently assigning the wrong annotation.

source(file.path("config", "project_paths.R"))
source(file.path("R", "project_io.R"))
source(file.path("R", "constants.R"))

suppressPackageStartupMessages({
  library(Seurat)
})

out_dir <- ensure_dir(file.path(RESULTS_DIR, "annotation"))
map_file <- file.path(out_dir, "merged_celltype_annotation_map.rds")
report_file <- file.path(out_dir, "core_annotation_mapping_report.csv")

# -----------------------------------------------------------------------------
# Stage 1: load ONLY the merged annotated object and extract a lightweight map.
# -----------------------------------------------------------------------------
message("Stage 1/3: extracting annotation map from merged_obj_annotated.rds ...")
merged <- load_merged_object()
require_metadata(merged, c("patient_id", "NF1", "celltype", "major_celltype"))

mmd <- merged@meta.data

# Candidate core/sample identifiers. These are carried into the lightweight map
# only when they actually exist in the merged object's metadata.
candidate_id_fields <- intersect(
  c(
    "core_id", "core", "core_name",
    "sample_id", "sample", "sample_name",
    "region_id", "region", "orig.ident"
  ),
  colnames(mmd)
)

mmap <- data.frame(
  merged_cell = rownames(mmd),
  patient_id = as.character(mmd$patient_id),
  NF1 = as.character(mmd$NF1),
  barcode = short_barcode(rownames(mmd)),
  celltype = as.character(mmd$celltype),
  major_celltype = as.character(mmd$major_celltype),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

for (field in candidate_id_fields) {
  mmap[[field]] <- as.character(mmd[[field]])
}

saveRDS(mmap, map_file, compress = FALSE)
message("Lightweight annotation map saved: ", map_file)
message(
  "Merged map: ", format(nrow(mmap), big.mark = ","),
  " cells; candidate core/sample fields: ",
  if (length(candidate_id_fields)) paste(candidate_id_fields, collapse = ", ") else "none"
)

# Critical memory step: remove the million-cell Seurat object BEFORE loading
# xenium_list.rds.
rm(merged, mmd)
invisible(gc())

# -----------------------------------------------------------------------------
# Helpers for safe matching.
# -----------------------------------------------------------------------------
try_barcode_map <- function(core_barcodes, candidate_map, method_label) {
  relevant <- candidate_map[candidate_map$barcode %in% core_barcodes, , drop = FALSE]

  if (!nrow(relevant)) {
    return(list(ok = FALSE, matched = 0L, idx = NULL, method = method_label,
                reason = "no matching barcodes"))
  }

  dup_barcode <- duplicated(relevant$barcode) |
    duplicated(relevant$barcode, fromLast = TRUE)

  ambiguous <- unique(relevant$barcode[dup_barcode])
  ambiguous <- intersect(ambiguous, core_barcodes)

  if (length(ambiguous)) {
    return(list(
      ok = FALSE,
      matched = 0L,
      idx = NULL,
      method = method_label,
      reason = paste0(length(ambiguous), " ambiguous duplicated barcode(s)")
    ))
  }

  idx <- match(core_barcodes, candidate_map$barcode)
  list(
    ok = TRUE,
    matched = sum(!is.na(idx)),
    idx = idx,
    method = method_label,
    reason = ""
  )
}

# -----------------------------------------------------------------------------
# Stage 2: load ONLY the clean per-core list and map regenerated annotations.
# -----------------------------------------------------------------------------
message("Stage 2/3: loading clean xenium_list.rds ...")
xlist <- load_xenium_list_input()

mapping_report <- vector("list", length(xlist))

for (i in seq_along(xlist)) {
  obj <- xlist[[i]]
  require_metadata(obj, c("patient_id", "NF1"))

  pid_values <- unique(na.omit(as.character(obj$patient_id)))
  if (length(pid_values) != 1L) {
    stop(
      "Core ", i, " must contain exactly one patient_id; found: ",
      paste(pid_values, collapse = ", ")
    )
  }
  pid <- pid_values[[1]]
  core_name <- if (!is.null(names(xlist)) && nzchar(names(xlist)[i])) {
    names(xlist)[i]
  } else {
    paste0("core_", i)
  }

  core_cells <- rownames(obj@meta.data)
  core_barcodes <- short_barcode(core_cells)

  # Candidate strategies are scored by the number of cells matched. We prefer
  # direct full-cell-name matching when available.
  strategies <- list()

  direct_idx <- match(core_cells, mmap$merged_cell)
  strategies[[length(strategies) + 1L]] <- list(
    ok = any(!is.na(direct_idx)),
    matched = sum(!is.na(direct_idx)),
    idx = direct_idx,
    method = "direct_cell_name",
    reason = ""
  )

  patient_map <- mmap[mmap$patient_id == pid, , drop = FALSE]
  if (!nrow(patient_map)) {
    stop("No cells for patient_id '", pid, "' were found in the merged annotation map.")
  }

  strategies[[length(strategies) + 1L]] <- try_barcode_map(
    core_barcodes,
    patient_map,
    "patient_id+barcode"
  )

  # If patient+barcode is ambiguous because a patient contributes more than one
  # core, try any shared metadata field that uniquely identifies this core.
  shared_fields <- intersect(candidate_id_fields, colnames(obj@meta.data))

  for (field in shared_fields) {
    core_value <- unique(na.omit(as.character(obj@meta.data[[field]])))
    if (length(core_value) != 1L) next

    restricted <- patient_map[
      !is.na(patient_map[[field]]) & patient_map[[field]] == core_value,
      ,
      drop = FALSE
    ]
    if (!nrow(restricted)) next

    strategies[[length(strategies) + 1L]] <- try_barcode_map(
      core_barcodes,
      restricted,
      paste0("patient_id+", field, "+barcode")
    )
  }

  # Also test the Xenium-list element name against candidate core/sample fields.
  if (!is.null(core_name) && nzchar(core_name)) {
    for (field in candidate_id_fields) {
      restricted <- patient_map[
        !is.na(patient_map[[field]]) & patient_map[[field]] == core_name,
        ,
        drop = FALSE
      ]
      if (!nrow(restricted)) next

      strategies[[length(strategies) + 1L]] <- try_barcode_map(
        core_barcodes,
        restricted,
        paste0("patient_id+", field, "(list_name)+barcode")
      )
    }
  }

  usable <- which(vapply(strategies, function(z) isTRUE(z$ok), logical(1)))
  if (!length(usable)) {
    reasons <- paste(
      vapply(strategies, function(z) paste0(z$method, ": ", z$reason), character(1)),
      collapse = "; "
    )
    stop(
      "Could not map core '", core_name, "' (patient_id=", pid,
      ") without ambiguous barcode assignments. Tried: ", reasons,
      ". A core-specific identifier is required before proceeding."
    )
  }

  matched_counts <- vapply(strategies[usable], function(z) z$matched, numeric(1))
  best <- strategies[[usable[[which.max(matched_counts)]]]]

  if (best$method == "direct_cell_name") {
    idx <- best$idx
  } else {
    # Reconstruct the exact candidate subset corresponding to the selected
    # strategy, then convert local indices back into mmap row indices.
    if (best$method == "patient_id+barcode") {
      candidate_rows <- which(mmap$patient_id == pid)
    } else {
      # Identify the field encoded in the method name.
      matched_field <- shared_fields[
        vapply(
          shared_fields,
          function(f) grepl(paste0("patient_id+", f, "+barcode"), best$method, fixed = TRUE),
          logical(1)
        )
      ]

      if (length(matched_field)) {
        field <- matched_field[[1]]
        core_value <- unique(na.omit(as.character(obj@meta.data[[field]])))[[1]]
        candidate_rows <- which(
          mmap$patient_id == pid &
            !is.na(mmap[[field]]) &
            mmap[[field]] == core_value
        )
      } else {
        list_fields <- candidate_id_fields[
          vapply(
            candidate_id_fields,
            function(f) grepl(paste0("patient_id+", f, "(list_name)+barcode"), best$method, fixed = TRUE),
            logical(1)
          )
        ]
        if (!length(list_fields)) stop("Internal mapping-method parsing error for ", best$method)
        field <- list_fields[[1]]
        candidate_rows <- which(
          mmap$patient_id == pid &
            !is.na(mmap[[field]]) &
            mmap[[field]] == core_name
        )
      }
    }

    candidate_map <- mmap[candidate_rows, , drop = FALSE]
    local_idx <- match(core_barcodes, candidate_map$barcode)
    idx <- rep(NA_integer_, length(local_idx))
    ok <- !is.na(local_idx)
    idx[ok] <- candidate_rows[local_idx[ok]]
  }

  obj$celltype <- mmap$celltype[idx]
  obj$major_celltype <- mmap$major_celltype[idx]
  obj$celltype <- factor(obj$celltype, levels = CELLTYPE_ORDER)
  obj$major_celltype <- factor(obj$major_celltype, levels = BROAD_CELLTYPE_ORDER)

  n_cells <- ncol(obj)
  n_matched <- sum(!is.na(idx))
  pct_matched <- if (n_cells) 100 * n_matched / n_cells else NA_real_

  mapping_report[[i]] <- data.frame(
    xenium_index = i,
    core_name = core_name,
    patient_id = pid,
    n_cells = n_cells,
    n_matched = n_matched,
    percent_matched = pct_matched,
    mapping_method = best$method,
    stringsAsFactors = FALSE
  )

  message(
    "Mapped ", i, "/", length(xlist), " [", core_name, "]: ",
    format(n_matched, big.mark = ","), "/", format(n_cells, big.mark = ","),
    " cells (", sprintf("%.2f", pct_matched), "%) using ", best$method
  )

  xlist[[i]] <- obj
  rm(obj)
}

mapping_report <- do.call(rbind, mapping_report)
write.csv(mapping_report, report_file, row.names = FALSE)

message("Stage 3/3: saving annotated per-core list ...")
saveRDS(xlist, XENIUM_LIST_ANNOTATED_FILE, compress = FALSE)
message("Annotated per-core list saved: ", XENIUM_LIST_ANNOTATED_FILE)
message("Mapping QC report saved: ", report_file)
message(
  "Overall mapping: ",
  format(sum(mapping_report$n_matched), big.mark = ","), "/",
  format(sum(mapping_report$n_cells), big.mark = ","),
  " cells (",
  sprintf("%.2f", 100 * sum(mapping_report$n_matched) / sum(mapping_report$n_cells)),
  "%)."
)
