# Figures 1G-I and downstream analyses: exact reproduction of the final 40-um
# neighborhood / 12 meta-niche workflow.
#
# IMPORTANT PROVENANCE DETAIL:
# The final Niche_revision.R computed frNN() on ALL cells present in each
# per-core Xenium object, including cells whose celltype was NA because they
# were not retained in merged_obj_filtered.rds.
#
# Therefore this script:
#   1) uses all 1,086,037 per-core cells as spatial points / neighborhood centers;
#   2) removes the center cell from each 40-um neighborhood;
#   3) requires >= 3 spatial neighbors, regardless of whether those neighbors
#      have a celltype annotation;
#   4) computes composition from NON-NA neighbor celltypes only (R table()
#      drops NA by default), exactly as the original code did;
#   5) allows an unannotated center cell to receive a meta-niche;
#   6) uses round(..., 6), set.seed(123), k=12, nstart=25;
#   7) maps the regenerated MNs back to all cores, then to retained merged cells.
#
# Starts from:
#   data/processed/xenium_list_annotated.rds
# which contains regenerated celltype annotations for the 1,003,604 retained
# cells and NA celltype for the 82,433 cells that were not in merged_obj_filtered.
#
# It does NOT use any old meta_niche annotation.

source(file.path("config", "project_paths.R"))
source(file.path("R", "project_io.R"))
source(file.path("R", "constants.R"))

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(dbscan)
})

out_dir <- ensure_dir(file.path(RESULTS_DIR, "meta_niches"))

message("Stage 1/5: loading regenerated per-core cell-type annotations ...")
xenium_list <- load_xenium_list()

# Remove any accidental legacy MN metadata.
legacy_mn_fields <- c("meta_niche", "meta_niche_name", "meta_niche_short")
for (i in seq_along(xenium_list)) {
  drop <- intersect(legacy_mn_fields, colnames(xenium_list[[i]]@meta.data))
  if (length(drop)) {
    xenium_list[[i]]@meta.data <- xenium_list[[i]]@meta.data[
      , !colnames(xenium_list[[i]]@meta.data) %in% drop, drop = FALSE
    ]
  }
}

# Validate the regenerated annotation coverage.
map_file <- file.path(RESULTS_DIR, "annotation", "merged_celltype_annotation_map.rds")
expected_annotated <- if (file.exists(map_file)) nrow(readRDS(map_file)) else NA_integer_

n_all_core_cells <- sum(vapply(xenium_list, ncol, integer(1)))
n_annotated <- sum(vapply(
  xenium_list,
  function(o) sum(!is.na(o$celltype) & nzchar(as.character(o$celltype))),
  integer(1)
))
n_unannotated <- n_all_core_cells - n_annotated

message(sprintf(
  "Core list: %s total cells; %s regenerated celltype annotations; %s cells with NA celltype.",
  format(n_all_core_cells, big.mark = ","),
  format(n_annotated, big.mark = ","),
  format(n_unannotated, big.mark = ",")
))

if (!is.na(expected_annotated) && n_annotated != expected_annotated) {
  stop(
    "Annotation coverage mismatch: ", n_annotated,
    " annotated core cells versus ", expected_annotated,
    " cells in the merged annotation map."
  )
}

message(
  "Exact-reproduction mode: ALL per-core cells will participate in spatial ",
  "neighborhood geometry; NA celltypes will not contribute to composition."
)

celltype_levels <- CELLTYPE_ORDER

get_cell_column <- function(coords) {
  if ("cell" %in% colnames(coords)) return(as.character(coords$cell))
  rn <- rownames(coords)
  if (!is.null(rn) && length(rn) == nrow(coords) && all(nzchar(rn))) return(rn)
  stop("GetTissueCoordinates() did not provide a cell column or usable row names.")
}

compute_one_exact <- function(obj, sample_index) {
  require_metadata(obj, c("patient_id", "NF1", "celltype"))

  md <- obj@meta.data
  coords <- GetTissueCoordinates(obj)

  if (!all(c("x", "y") %in% colnames(coords))) {
    stop("Core ", sample_index, " lacks x/y tissue coordinates.")
  }

  coord_cells <- get_cell_column(coords)
  idx <- match(coord_cells, rownames(md))
  if (anyNA(idx)) stop("Unmatched coordinate cells in core ", sample_index, ".")

  md_use <- md[idx, , drop = FALSE]

  if (nrow(coords) != nrow(md_use)) {
    stop("Coordinate/metadata row mismatch in core ", sample_index, ".")
  }

  # All spatial cells are included here, matching the final original workflow.
  xy <- as.matrix(coords[, c("x", "y"), drop = FALSE])
  nn <- dbscan::frNN(xy, eps = NEIGHBORHOOD_RADIUS_UM, sort = FALSE)

  celltypes <- as.character(md_use$celltype)

  # Guard against unexpected non-NA labels.
  unexpected <- setdiff(unique(na.omit(celltypes)), celltype_levels)
  if (length(unexpected)) {
    stop(
      "Unexpected celltype label(s) in core ", sample_index, ": ",
      paste(unexpected, collapse = ", ")
    )
  }

  comp <- matrix(
    0,
    nrow = nrow(coords),
    ncol = length(celltype_levels),
    dimnames = list(coord_cells, celltype_levels)
  )

  valid <- logical(nrow(coords))
  neighbor_n <- integer(nrow(coords))
  annotated_neighbor_n <- integer(nrow(coords))

  for (j in seq_len(nrow(coords))) {
    nbr <- nn$id[[j]]

    # Original code explicitly removed self.
    self_pos <- which(nbr == j)
    if (length(self_pos) > 0L) nbr <- nbr[-self_pos]

    neighbor_n[j] <- length(nbr)

    # The threshold is based on ALL spatial neighbors, including NA-celltype cells.
    if (length(nbr) < MIN_NEIGHBORS) next

    # This intentionally mirrors:
    #   tb <- table(celltypes[nbr_idx])
    # R's table() drops NA values. Thus unannotated neighbors are spatially present
    # but do not contribute to the cell-type fraction denominator.
    tb <- table(celltypes[nbr])
    annotated_neighbor_n[j] <- sum(tb)

    # If every neighbor has NA celltype, the original zero-initialized composition
    # row remains all zeros and is still a valid neighborhood.
    if (length(tb) > 0L && sum(tb) > 0L) {
      comp[j, names(tb)] <- as.numeric(tb / sum(tb))
    }

    valid[j] <- TRUE
  }

  df <- as.data.frame(comp[valid, , drop = FALSE], check.names = FALSE)
  df$cell <- coord_cells[valid]
  df$patient_id <- as.character(md_use$patient_id[valid])
  df$NF1 <- as.character(md_use$NF1[valid])
  df$sample_index <- sample_index
  df$n_neighbors_40um <- neighbor_n[valid]
  df$n_annotated_neighbors_40um <- annotated_neighbor_n[valid]
  df$center_celltype_annotated <- !is.na(celltypes[valid]) & nzchar(celltypes[valid])

  qc <- data.frame(
    sample_index = sample_index,
    patient_id = paste(unique(as.character(md_use$patient_id)), collapse = ";"),
    NF1 = paste(unique(as.character(md_use$NF1)), collapse = ";"),
    n_core_cells = ncol(obj),
    n_celltype_annotated = sum(!is.na(celltypes) & nzchar(celltypes)),
    n_celltype_NA = sum(is.na(celltypes) | !nzchar(celltypes)),
    n_valid_neighborhoods = sum(valid),
    n_excluded_lt3_neighbors = sum(!valid),
    pct_valid = round(100 * mean(valid), 3),
    stringsAsFactors = FALSE
  )

  list(neighborhoods = df, qc = qc)
}

message("Stage 2/5: calculating 40-um neighborhoods on ALL per-core cells ...")
res <- vector("list", length(xenium_list))

for (i in seq_along(xenium_list)) {
  core_label <- names(xenium_list)[i]
  if (is.null(core_label) || is.na(core_label) || !nzchar(core_label)) {
    core_label <- paste0("core_", i)
  }

  message(sprintf(
    "Neighborhoods %d/%d [%s]",
    i, length(xenium_list), core_label
  ))

  res[[i]] <- compute_one_exact(xenium_list[[i]], i)
}

cell_neighborhood_df <- bind_rows(lapply(res, `[[`, "neighborhoods"))
neighbor_qc <- bind_rows(lapply(res, `[[`, "qc"))
rm(res)
gc()

write.csv(
  neighbor_qc,
  file.path(out_dir, "neighborhood_40um_QC_by_core.csv"),
  row.names = FALSE
)

message(sprintf(
  "Valid neighborhoods: %s; excluded for < %d spatial neighbors: %s.",
  format(nrow(cell_neighborhood_df), big.mark = ","),
  MIN_NEIGHBORS,
  format(n_all_core_cells - nrow(cell_neighborhood_df), big.mark = ",")
))

if (nrow(cell_neighborhood_df) == 0L) stop("No valid neighborhoods were generated.")

message("Stage 3/5: collapsing identical vectors and running k-means (k=12) ...")

meta_cols <- c(
  "cell", "patient_id", "NF1", "sample_index",
  "n_neighbors_40um", "n_annotated_neighbors_40um",
  "center_celltype_annotated"
)
celltype_cols <- intersect(celltype_levels, colnames(cell_neighborhood_df))
cell_mat <- as.matrix(cell_neighborhood_df[, celltype_cols, drop = FALSE])

# Exact final logic: round to 6 decimals, collapse identical vectors.
rounded_mat <- round(cell_mat, 6)
rowkey <- apply(rounded_mat, 1, paste, collapse = "_")
unique_rowkeys <- unique(rowkey)
first_idx <- match(unique_rowkeys, rowkey)
unique_mat <- rounded_mat[first_idx, , drop = FALSE]

if (nrow(unique_mat) < N_META_NICHES) {
  stop(
    "Only ", nrow(unique_mat),
    " unique neighborhood vectors; cannot fit k=", N_META_NICHES, "."
  )
}

set.seed(META_NICHE_SEED)
km <- kmeans(
  unique_mat,
  centers = N_META_NICHES,
  nstart = META_NICHE_NSTART
)

unique_meta <- paste0("MetaNiche_", km$cluster)
cell_neighborhood_df$meta_niche <- unique_meta[
  match(rowkey, unique_rowkeys)
]

# Save k-means centroids for reproducibility/auditing.
centers_df <- as.data.frame(km$centers, check.names = FALSE)
centers_df$meta_niche <- paste0("MetaNiche_", seq_len(nrow(centers_df)))
centers_df <- centers_df[, c("meta_niche", celltype_cols), drop = FALSE]
write.csv(
  centers_df,
  file.path(out_dir, "meta_niche_kmeans_centers.csv"),
  row.names = FALSE
)

rm(cell_mat, rounded_mat, rowkey, unique_rowkeys, first_idx, unique_mat, km, unique_meta)
gc()

# Map freshly generated MNs to ALL per-core cells.
for (i in seq_along(xenium_list)) {
  sub <- cell_neighborhood_df[
    cell_neighborhood_df$sample_index == i,
    c("cell", "meta_niche"),
    drop = FALSE
  ]

  mn_map <- setNames(sub$meta_niche, sub$cell)
  cells <- rownames(xenium_list[[i]]@meta.data)
  xenium_list[[i]]$meta_niche <- unname(mn_map[cells])
}

# Patient-level abundance exactly over all valid neighborhood centers.
meta_niche_freq <- cell_neighborhood_df |>
  count(patient_id, NF1, meta_niche, name = "count") |>
  group_by(patient_id, NF1) |>
  mutate(freq = count / sum(count)) |>
  ungroup()

meta_niche_stats <- meta_niche_freq |>
  group_by(meta_niche) |>
  filter(sum(NF1 == "Mut") >= 2, sum(NF1 == "WT") >= 2) |>
  summarise(
    p_val = wilcox.test(freq ~ NF1)$p.value,
    mean_mut = mean(freq[NF1 == "Mut"]),
    mean_wt = mean(freq[NF1 == "WT"]),
    n_mut = sum(NF1 == "Mut"),
    n_wt = sum(NF1 == "WT"),
    .groups = "drop"
  ) |>
  mutate(adj_p = p.adjust(p_val, method = "BH"))

mn_prevalence_qc <- cell_neighborhood_df |>
  distinct(patient_id, meta_niche) |>
  count(meta_niche, name = "n_patient_or_specimen_ids") |>
  arrange(as.integer(sub("MetaNiche_", "", meta_niche)))

write.csv(
  meta_niche_freq,
  file.path(out_dir, "meta_niche_frequency_per_patient.csv"),
  row.names = FALSE
)
write.csv(
  meta_niche_stats,
  file.path(out_dir, "meta_niche_nf1_abundance_stats.csv"),
  row.names = FALSE
)
write.csv(
  mn_prevalence_qc,
  file.path(out_dir, "meta_niche_prevalence_QC.csv"),
  row.names = FALSE
)

# Efficient checkpoint first.
saveRDS(
  cell_neighborhood_df,
  file.path(out_dir, "cell_neighborhoods_40um_k12.rds"),
  compress = FALSE
)

message("Writing neighborhood CSV (this can take time for ~1.1M rows) ...")
write.csv(
  cell_neighborhood_df,
  file.path(out_dir, "cell_neighborhoods_40um_k12.csv"),
  row.names = FALSE
)

message("Stage 4/5: saving per-core meta-niche checkpoint ...")
saveRDS(xenium_list, XENIUM_LIST_MN_FILE, compress = FALSE)

# Build map only for the retained/annotated cells, because those are the cells
# present in merged_obj_annotated.rds.
mn_map_df <- bind_rows(lapply(xenium_list, function(o) {
  md <- o@meta.data
  keep <- !is.na(md$celltype) & nzchar(as.character(md$celltype))

  data.frame(
    key = paste(
      as.character(md$patient_id[keep]),
      short_barcode(rownames(md)[keep]),
      sep = "__"
    ),
    meta_niche = as.character(md$meta_niche[keep]),
    stringsAsFactors = FALSE
  )
}))

if (anyDuplicated(mn_map_df$key)) {
  stop(
    "Duplicate patient+barcode keys among retained cells while mapping ",
    "meta-niches back to the merged object."
  )
}

saveRDS(
  mn_map_df,
  file.path(out_dir, "meta_niche_cell_map_retained.rds"),
  compress = FALSE
)

rm(xenium_list)
gc()

message("Stage 5/5: mapping regenerated MNs back to merged_obj_annotated.rds ...")
merged <- load_merged_object()

mkey <- paste(
  as.character(merged$patient_id),
  short_barcode(rownames(merged@meta.data)),
  sep = "__"
)

if (anyDuplicated(mkey)) {
  stop("Duplicate patient+barcode keys in merged object; core-specific matching is required.")
}

idx <- match(mkey, mn_map_df$key)
if (anyNA(idx)) {
  stop(
    sum(is.na(idx)),
    " merged cells could not be matched back to the retained-cell MN map."
  )
}

merged$meta_niche <- mn_map_df$meta_niche[idx]
saveRDS(merged, MN_OBJECT_FILE, compress = FALSE)

message("Exact meta-niche reconstruction complete.")
message("Saved per-core MN list: ", XENIUM_LIST_MN_FILE)
message("Saved merged MN object: ", MN_OBJECT_FILE)
message("Valid neighborhoods: ", format(nrow(cell_neighborhood_df), big.mark = ","))
message("Expected historical validation target, if using the same final input: 1,085,053 valid neighborhoods.")
message("QC: ", file.path(out_dir, "meta_niche_prevalence_QC.csv"))
