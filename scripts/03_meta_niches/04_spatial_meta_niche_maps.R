# ==============================================================================
# 04_spatial_meta_niche_maps.R
#
# Spatial maps of ALL 12 publication-reference meta-niches.
#
# Uses:
#   data/processed/xenium_list_with_meta_niche_published_reference.rds
#
# Produces:
#   results/meta_niches/spatial_maps/all_meta_niches/
#   results/meta_niches/spatial_meta_niche_core_summary_published_reference.csv
#   results/meta_niches/spatial_meta_niche_map_QC_published_reference.csv
#
# IMPORTANT
# ---------
# This CLEAN script shows ALL 12 MNs in each analyzed Xenium tissue core.
# It intentionally does NOT generate selected/highlighted MN2/MN3/MN4/MN6 maps.
#
# The manuscript Figure 1I contains representative FOVs. This script generates
# the complete all-MN map set from which representative FOVs can be inspected.
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
})

# ------------------------------------------------------------------------------
# 0. Configuration
# ------------------------------------------------------------------------------

xlist_file <- file.path(
  "data",
  "processed",
  "xenium_list_with_meta_niche_published_reference.rds"
)

out_all <- file.path(
  "results",
  "meta_niches",
  "spatial_maps",
  "all_meta_niches"
)

dir.create(
  out_all,
  recursive = TRUE,
  showWarnings = FALSE
)

# Archived MN colors.
NICHE_COLORS <- c(
  "MetaNiche_1"  = "#b2df8a",
  "MetaNiche_2"  = "#33a02c",
  "MetaNiche_3"  = "#e31a1c",
  "MetaNiche_4"  = "#ff7f00",
  "MetaNiche_5"  = "#fdbf6f",
  "MetaNiche_6"  = "#b15928",
  "MetaNiche_7"  = "#a6cee3",
  "MetaNiche_8"  = "#ffff99",
  "MetaNiche_9"  = "#fb9a99",
  "MetaNiche_10" = "#6a3d9a",
  "MetaNiche_11" = "#cab2d6",
  "MetaNiche_12" = "#1f78b4"
)

MN_LEVELS <- names(NICHE_COLORS)

# ------------------------------------------------------------------------------
# 1. Helpers
# ------------------------------------------------------------------------------

sanitize_filename <- function(x) {
  x <- gsub("[^A-Za-z0-9._-]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  x
}

normalize_nf1 <- function(x) {
  x <- as.character(x)

  out <- dplyr::case_when(
    x %in% c("Mut", "NF1Mut", "NF1_Mut", "NF1 Mut") ~ "Mut",
    x %in% c("WT", "NF1WT", "NF1_WT", "NF1 WT") ~ "WT",
    TRUE ~ NA_character_
  )

  if (anyNA(out)) {
    stop(
      "Unexpected NF1 value(s): ",
      paste(sort(unique(x[is.na(out)])), collapse = ", ")
    )
  }

  out
}

single_value <- function(x, field, sample_index) {
  x <- unique(
    as.character(
      x[!is.na(x)]
    )
  )

  if (length(x) != 1L) {
    stop(
      "Expected one unique ",
      field,
      " in sample_index ",
      sample_index,
      "; found ",
      paste(x, collapse = ", ")
    )
  }

  x
}

# ------------------------------------------------------------------------------
# 2. Load publication-reference per-core objects
# ------------------------------------------------------------------------------

if (!file.exists(xlist_file)) {
  stop(
    "Missing publication-reference Xenium list: ",
    xlist_file,
    "\nRun 02_apply_published_meta_niche_reference.R first."
  )
}

cat("Loading publication-reference Xenium list...\n")
xenium_list <- readRDS(xlist_file)

if (length(xenium_list) != 39L) {
  stop(
    "Expected 39 analyzed Xenium tissue cores; found ",
    length(xenium_list),
    "."
  )
}

cat(
  "Loaded ",
  length(xenium_list),
  " tissue-core objects.\n",
  sep = ""
)

# ------------------------------------------------------------------------------
# 3. Per-core summary of ALL 12 MNs
# ------------------------------------------------------------------------------

cat("\nBuilding per-core MN summary...\n")

core_summary_list <- vector(
  "list",
  length(xenium_list)
)

for (i in seq_along(xenium_list)) {

  obj <- xenium_list[[i]]
  md <- obj@meta.data

  required <- c(
    "patient_id",
    "NF1",
    "meta_niche"
  )

  missing_cols <- setdiff(
    required,
    colnames(md)
  )

  if (length(missing_cols) > 0L) {
    stop(
      "sample_index ",
      i,
      " missing metadata: ",
      paste(missing_cols, collapse = ", ")
    )
  }

  sample_id <- single_value(
    md$patient_id,
    "patient_id",
    i
  )

  nf1 <- single_value(
    normalize_nf1(md$NF1),
    "NF1",
    i
  )

  mn <- as.character(md$meta_niche)
  valid <- !is.na(mn)

  if (!any(valid)) {
    stop(
      "No valid MN assignments in sample_index ",
      i
    )
  }

  tab <- table(
    factor(
      mn[valid],
      levels = MN_LEVELS
    )
  )

  frac <- as.numeric(tab) / sum(tab)

  core_summary_list[[i]] <- data.frame(
    sample_index = i,
    sample_id = sample_id,
    NF1 = nf1,
    total_cells_in_core_object = nrow(md),
    valid_MN_neighborhoods = sum(valid),
    unassigned_MN_cells = sum(!valid),
    meta_niche = MN_LEVELS,
    n_neighborhoods = as.integer(tab),
    fraction = frac,
    stringsAsFactors = FALSE
  )
}

core_summary <- dplyr::bind_rows(
  core_summary_list
)

write.csv(
  core_summary,
  file.path(
    "results",
    "meta_niches",
    "spatial_meta_niche_core_summary_published_reference.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------------------------
# 4. Generate ALL-MN spatial maps
# ------------------------------------------------------------------------------

cat("\nGenerating all-12-MN spatial maps...\n")

plot_qc <- vector(
  "list",
  length(xenium_list)
)

for (i in seq_along(xenium_list)) {

  obj <- xenium_list[[i]]
  md <- obj@meta.data

  sample_id <- single_value(
    md$patient_id,
    "patient_id",
    i
  )

  nf1 <- single_value(
    normalize_nf1(md$NF1),
    "NF1",
    i
  )

  file_id <- sanitize_filename(
    paste0(
      sprintf("core%02d", i),
      "_",
      sample_id,
      "_NF1",
      nf1
    )
  )

  # Fixed 12-MN ordering and colors.
  obj$meta_niche_plot <- factor(
    as.character(obj$meta_niche),
    levels = MN_LEVELS
  )

  all_ok <- TRUE
  all_message <- NA_character_

  tryCatch(
    {
      p_all <- ImageDimPlot(
        obj,
        group.by = "meta_niche_plot",
        cols = NICHE_COLORS,
        alpha = 0.9
      ) +
        ggtitle(
          paste0(
            sample_id,
            " | NF1 ",
            nf1
          )
        ) +
        theme(
          legend.position = "right",
          plot.title = element_text(
            face = "bold",
            hjust = 0.5
          )
        )

      ggsave(
        filename = file.path(
          out_all,
          paste0(
            file_id,
            "_all_meta_niches.png"
          )
        ),
        plot = p_all,
        width = 8,
        height = 7,
        dpi = 300
      )

      rm(p_all)
    },
    error = function(e) {
      all_ok <<- FALSE
      all_message <<- conditionMessage(e)
    }
  )

  plot_qc[[i]] <- data.frame(
    sample_index = i,
    sample_id = sample_id,
    NF1 = nf1,
    all_MN_map_ok = all_ok,
    all_MN_message = all_message,
    stringsAsFactors = FALSE
  )

  cat(
    sprintf(
      "  core %02d | %-10s | NF1 %-3s | all-MN map=%s\n",
      i,
      sample_id,
      nf1,
      ifelse(all_ok, "OK", "FAIL")
    )
  )

  rm(obj)
  gc(verbose = FALSE)
}

plot_qc <- dplyr::bind_rows(
  plot_qc
)

write.csv(
  plot_qc,
  file.path(
    "results",
    "meta_niches",
    "spatial_meta_niche_map_QC_published_reference.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------------------------
# 5. Final QC
# ------------------------------------------------------------------------------

n_all_ok <- sum(
  plot_qc$all_MN_map_ok
)

cat("\nSpatial-map QC:\n")
cat(
  "  all-MN maps successful: ",
  n_all_ok,
  "/",
  nrow(plot_qc),
  "\n",
  sep = ""
)

if (n_all_ok != 39L) {
  warning(
    "One or more all-MN maps failed. ",
    "Check spatial_meta_niche_map_QC_published_reference.csv"
  )
}

cat(
  "\nDONE: all-12-MN spatial maps complete. ",
  "No selected/highlighted MN maps were generated.\n",
  sep = ""
)
