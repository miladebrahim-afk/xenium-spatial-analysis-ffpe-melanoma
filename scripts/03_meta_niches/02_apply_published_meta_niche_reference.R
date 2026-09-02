# Apply the publication-reference meta-niche assignments
#
# The de novo reconstruction remains untouched. This script applies only the
# documented 140 boundary-cell differences needed to reproduce the archived
# final manuscript MN assignments exactly.
#
# Inputs:
#   data/processed/xenium_list_with_meta_niche.rds
#   data/processed/merged_obj_annotated_with_meta_niche.rds
#   data/reference/published_meta_niche_boundary_corrections.csv
#
# Outputs:
#   data/processed/xenium_list_with_meta_niche_published_reference.rds
#   data/processed/merged_obj_annotated_with_meta_niche_published_reference.rds
#   results/meta_niches/figure1H_core_occurrence_published_reference.csv
#   results/meta_niches/published_reference_application_QC.csv

suppressPackageStartupMessages({
  library(dplyr)
})

xlist_in <- file.path("data", "processed", "xenium_list_with_meta_niche.rds")
merged_in <- file.path("data", "processed", "merged_obj_annotated_with_meta_niche.rds")
correction_file <- file.path(
  "data", "reference", "published_meta_niche_boundary_corrections.csv"
)

xlist_out <- file.path(
  "data", "processed", "xenium_list_with_meta_niche_published_reference.rds"
)
merged_out <- file.path(
  "data", "processed", "merged_obj_annotated_with_meta_niche_published_reference.rds"
)

mn_dir <- file.path("results", "meta_niches")
dir.create(mn_dir, recursive = TRUE, showWarnings = FALSE)

expected_occurrence <- c(
  MetaNiche_1  = 37L,
  MetaNiche_2  = 39L,
  MetaNiche_3  = 39L,
  MetaNiche_4  = 39L,
  MetaNiche_5  = 33L,
  MetaNiche_6  = 20L,
  MetaNiche_7  = 39L,
  MetaNiche_8  = 39L,
  MetaNiche_9  = 26L,
  MetaNiche_10 = 30L,
  MetaNiche_11 = 33L,
  MetaNiche_12 = 36L
)
mn_levels <- names(expected_occurrence)

same_label <- function(a, b) {
  (is.na(a) & is.na(b)) | (!is.na(a) & !is.na(b) & a == b)
}

short_barcode <- function(x) sub("^.*_", "", x)

if (!file.exists(xlist_in)) stop("Missing input: ", xlist_in)
if (!file.exists(merged_in)) stop("Missing input: ", merged_in)
if (!file.exists(correction_file)) stop("Missing reference: ", correction_file)

cat("Stage 1: load de novo per-core objects and reference table...\n")
xlist <- readRDS(xlist_in)
corr <- read.csv(correction_file, stringsAsFactors = FALSE)

required_cols <- c(
  "sample_index", "cell", "de_novo_meta_niche", "published_meta_niche"
)
if (!all(required_cols %in% colnames(corr))) {
  stop("Correction file is missing required columns.")
}
if (nrow(corr) != 140L) {
  stop("Expected 140 correction rows; found ", nrow(corr), ".")
}
if (length(xlist) != 39L) {
  stop("Expected 39 per-core objects; found ", length(xlist), ".")
}

corr$sample_index <- as.integer(corr$sample_index)

if (anyDuplicated(paste(corr$sample_index, corr$cell, sep = "__"))) {
  stop("Duplicate sample_index + cell keys in correction table.")
}

cat("Reference rows:", nrow(corr), "\n")

cat("\nStage 2: apply documented boundary corrections...\n")
applied <- 0L

for (i in sort(unique(corr$sample_index))) {
  ii <- which(corr$sample_index == i)
  md <- xlist[[i]]@meta.data
  pos <- match(corr$cell[ii], rownames(md))

  if (anyNA(pos)) {
    stop("Reference cell(s) missing from sample_index ", i)
  }

  current <- as.character(md$meta_niche[pos])
  expected_current <- corr$de_novo_meta_niche[ii]

  if (!all(same_label(current, expected_current))) {
    stop(
      "De novo-label validation failed in sample_index ", i,
      ". The correction table does not match this de novo object."
    )
  }

  md$meta_niche[pos] <- corr$published_meta_niche[ii]
  xlist[[i]]@meta.data <- md
  applied <- applied + length(ii)

  cat(sprintf("  core %02d: applied %d correction(s)\n", i, length(ii)))
}

cat("\nStage 3: verify exact Figure 1H core occurrence...\n")

occurrence <- bind_rows(lapply(xlist, function(obj) {
  md <- obj@meta.data
  data.frame(
    patient_or_specimen_id = as.character(md$patient_id),
    meta_niche = as.character(md$meta_niche),
    stringsAsFactors = FALSE
  )
})) |>
  filter(!is.na(meta_niche)) |>
  distinct(patient_or_specimen_id, meta_niche) |>
  count(meta_niche, name = "n_cores") |>
  arrange(as.integer(sub("MetaNiche_", "", meta_niche)))

print(occurrence, row.names = FALSE)

obs <- setNames(occurrence$n_cores, occurrence$meta_niche)[mn_levels]

if (!identical(as.integer(obs), as.integer(expected_occurrence[mn_levels]))) {
  stop(
    "Figure 1H verification failed. Observed: ",
    paste(obs, collapse = ", ")
  )
}

write.csv(
  occurrence,
  file.path(mn_dir, "figure1H_core_occurrence_published_reference.csv"),
  row.names = FALSE
)

cat("Figure 1H matches the archived final analysis exactly.\n")

cat("\nStage 4: save corrected per-core object...\n")
saveRDS(xlist, xlist_out, compress = FALSE)
cat("Saved: ", xlist_out, "\n", sep = "")

cat("\nStage 5: propagate applicable corrections to merged retained-cell object...\n")

# Resolve the patient/specimen ID for each documented correction while xlist
# is still available. This is done before loading the large merged object.
corr_patient <- vapply(seq_len(nrow(corr)), function(r) {
  i <- corr$sample_index[r]
  md_r <- xlist[[i]]@meta.data
  p <- match(corr$cell[r], rownames(md_r))

  if (is.na(p)) {
    stop(
      "Reference cell missing while resolving patient_id: ",
      corr$cell[r],
      " (sample_index ", i, ")."
    )
  }

  as.character(md_r$patient_id[p])
}, character(1))

corr_short <- short_barcode(corr$cell)
corr_key <- paste(corr_patient, corr_short, sep = "__")

if (anyDuplicated(corr_key)) {
  stop(
    "Correction patient_id + short-barcode mapping is ambiguous. ",
    "At least two correction rows resolve to the same merged-cell key."
  )
}

# The 39 per-core objects are no longer needed. Release them, plus the large
# loop metadata object that remains from Stage 2, BEFORE reading the merged RDS.
rm(xlist)
rm(list = intersect(
  c("md", "pos", "current", "expected_current", "ii", "i"),
  ls()
))
invisible(gc(full = TRUE))

cat("Per-core objects released from memory.\n")
cat("Loading merged retained-cell object...\n")

merged <- readRDS(merged_in)

if (!"patient_id" %in% colnames(merged@meta.data)) {
  stop("Merged object is missing patient_id.")
}
if (!"meta_niche" %in% colnames(merged@meta.data)) {
  stop("Merged object is missing meta_niche.")
}

# Low-memory matching:
# Scan merged cell names in chunks and create keys only for cells whose short
# barcode occurs among the 140 documented corrections. This avoids allocating
# one giant patient_id + barcode character vector for every retained cell.
merged_cells <- rownames(merged@meta.data)
merged_patient <- merged@meta.data$patient_id

n_merged <- length(merged_cells)
chunk_size <- 100000L

merged_pos <- rep(NA_integer_, nrow(corr))
match_count <- integer(nrow(corr))

cat(
  "Matching documented corrections against ",
  format(n_merged, big.mark = ","),
  " retained cells in chunks of ",
  format(chunk_size, big.mark = ","),
  "...\n",
  sep = ""
)

for (start in seq.int(1L, n_merged, by = chunk_size)) {
  end <- min(start + chunk_size - 1L, n_merged)
  idx <- start:end

  chunk_short <- short_barcode(merged_cells[idx])

  # Restrict key construction to barcodes that could correspond to a
  # documented correction.
  candidate_local <- which(chunk_short %in% corr_short)

  if (length(candidate_local) > 0L) {
    candidate_idx <- idx[candidate_local]
    candidate_key <- paste(
      as.character(merged_patient[candidate_idx]),
      chunk_short[candidate_local],
      sep = "__"
    )

    corr_idx <- match(candidate_key, corr_key)
    valid <- which(!is.na(corr_idx))

    if (length(valid) > 0L) {
      for (v in valid) {
        r <- corr_idx[v]
        match_count[r] <- match_count[r] + 1L

        if (match_count[r] == 1L) {
          merged_pos[r] <- candidate_idx[v]
        }
      }
    }
  }

  # Drop chunk temporaries before the next block.
  rm(idx, chunk_short, candidate_local)
  if (exists("candidate_idx")) rm(candidate_idx)
  if (exists("candidate_key")) rm(candidate_key)
  if (exists("corr_idx")) rm(corr_idx)
  if (exists("valid")) rm(valid)
}

if (any(match_count > 1L)) {
  ambiguous <- which(match_count > 1L)
  stop(
    "Merged patient_id + short-barcode mapping is ambiguous for ",
    length(ambiguous),
    " documented correction(s)."
  )
}

matched_corr <- which(match_count == 1L)

if (length(matched_corr) > 0L) {
  current_merged <- as.character(
    merged@meta.data$meta_niche[merged_pos[matched_corr]]
  )
  expected_de_novo <- corr$de_novo_meta_niche[matched_corr]

  if (!all(same_label(current_merged, expected_de_novo))) {
    bad <- matched_corr[
      !same_label(current_merged, expected_de_novo)
    ]
    stop(
      "Merged de novo-label validation failed for ",
      length(bad),
      " correction(s)."
    )
  }

  merged@meta.data$meta_niche[merged_pos[matched_corr]] <-
    corr$published_meta_niche[matched_corr]
}

# Remove full-length helper references before serializing the merged object.
rm(merged_cells, merged_patient, match_count)
invisible(gc(full = TRUE))

cat("Saving corrected merged retained-cell object...\n")
saveRDS(merged, merged_out, compress = FALSE)

cat("Corrections present in retained merged object:", length(matched_corr), "\n")
cat(
  "Corrections only in filtered-out per-core cells:",
  nrow(corr) - length(matched_corr), "\n"
)
cat("Saved: ", merged_out, "\n", sep = "")

qc <- data.frame(
  metric = c(
    "documented_boundary_corrections",
    "corrections_applied_to_per_core_objects",
    "corrections_present_in_retained_merged_object",
    "corrections_only_in_filtered_out_per_core_cells",
    "figure1H_exact_match"
  ),
  value = c(
    nrow(corr),
    applied,
    length(matched_corr),
    nrow(corr) - length(matched_corr),
    TRUE
  ),
  stringsAsFactors = FALSE
)

write.csv(
  qc,
  file.path(mn_dir, "published_reference_application_QC.csv"),
  row.names = FALSE
)

cat("\nApplication QC:\n")
print(qc, row.names = FALSE)

cat("\nDONE.\n")
cat("Use *_published_reference.rds for exact manuscript-figure reproduction.\n")
cat("Keep the unmodified de novo RDS files as the unsupervised reconstruction.\n")
