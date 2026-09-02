# ==============================================================================
# 02b_egfr_response_dotplots_from_cached_fig5b.R
#
# Figure 5C and Figure S6A
# EGFR prevalence/intensity by treatment response within meta-niches.
#
# FINAL CLEAN PIPELINE METHOD: uses the already-generated Figure 5B
# raw-count core-level EGFR metrics. No SCT extraction is performed.
#
# FINAL CACHED WORKFLOW
# ---------------------
# This script does NOT reopen the million-cell Seurat object.
# It reuses the already-generated Figure 5B core-level EGFR metrics:
#
#   Figure5B_core_level_EGFR_raw_count_metrics.csv
#
# and joins the historical treatment/response annotation by tissue-core ID.
#
# Figure 5C  = NF1Mut, responder vs non-responder
# Figure S6A = NF1WT,  responder vs non-responder
#
# Dot size  = pooled number of EGFR+ cells
# Dot color = pooled mean EGFR among EGFR+ cells
#
# INPUTS
# ------
# results/egfr/figure5BC_egfr_metaniche/
#   Figure5B_core_level_EGFR_raw_count_metrics.csv
#
# private_provenance/
#   patient_treatment_response_annotation.xlsx
#
# OUTPUTS
# -------
# results/egfr/figure5BC_egfr_metaniche/
#   Figure5C_NF1Mut_EGFR_by_response_meta_niche_FINAL.png/.pdf
#   FigureS6A_NF1WT_EGFR_by_response_meta_niche_FINAL.png/.pdf
#   Figure5C_S6A_EGFR_response_dotplot_values_FINAL.csv
#
# NOTE
# ----
# `patient_id` is used operationally as the tissue-core/specimen ID.
# The private annotation file is not written back out, so specimen IDs are not
# exposed in the public result table.
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readxl)
  library(ggplot2)
})

# ------------------------------------------------------------------------------
# 0. Paths
# ------------------------------------------------------------------------------

fig5b_core_file <- file.path(
  "results",
  "egfr",
  "figure5BC_egfr_metaniche",
  "Figure5B_core_level_EGFR_raw_count_metrics.csv"
)

response_file <- file.path(
  "private_provenance",
  "patient_treatment_response_annotation.xlsx"
)

out_dir <- file.path(
  "results",
  "egfr",
  "figure5BC_egfr_metaniche"
)

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

if (!file.exists(fig5b_core_file)) {
  stop("Missing Figure 5B core-level EGFR file: ", fig5b_core_file)
}

if (!file.exists(response_file)) {
  stop("Missing historical response annotation file: ", response_file)
}

# ------------------------------------------------------------------------------
# 1. Read cached Figure 5B core-level EGFR metrics
# ------------------------------------------------------------------------------

egfr_core <- read.csv(
  fig5b_core_file,
  stringsAsFactors = FALSE
)

required_core <- c(
  "sample_id",
  "NF1_clean",
  "meta_niche_clean",
  "n_EGFR_pos",
  "mean_EGFR_pos"
)

missing_core <- setdiff(
  required_core,
  colnames(egfr_core)
)

if (length(missing_core) > 0L) {
  stop(
    "Figure 5B core-level CSV is missing: ",
    paste(missing_core, collapse = ", ")
  )
}

egfr_core$sample_id <- as.character(
  egfr_core$sample_id
)

# ------------------------------------------------------------------------------
# 2. Read historical treatment/response annotations
# ------------------------------------------------------------------------------

annot <- readxl::read_xlsx(
  response_file
) |>
  as.data.frame()

required_annot <- c(
  "patient_id",
  "NF1",
  "treatment",
  "response"
)

missing_annot <- setdiff(
  required_annot,
  colnames(annot)
)

if (length(missing_annot) > 0L) {
  stop(
    "Response annotation file is missing: ",
    paste(missing_annot, collapse = ", ")
  )
}

normalize_nf1 <- function(x) {
  x <- trimws(as.character(x))

  dplyr::case_when(
    x %in% c("Mut", "NF1Mut", "NF1_Mut", "NF1 Mut") ~ "Mut",
    x %in% c("WT", "NF1WT", "NF1_WT", "NF1 WT") ~ "WT",
    TRUE ~ NA_character_
  )
}

normalize_response <- function(x) {
  x0 <- trimws(tolower(as.character(x)))

  dplyr::case_when(
    x0 %in% c(
      "r",
      "responder",
      "responders",
      "response",
      "responded"
    ) ~ "R",

    x0 %in% c(
      "nr",
      "non-responder",
      "non responder",
      "nonresponder",
      "non-responders",
      "non responders",
      "nonresponders",
      "no response",
      "not responder"
    ) ~ "NR",

    TRUE ~ NA_character_
  )
}

annot_clean <- annot |>
  dplyr::transmute(
    sample_id = as.character(patient_id),
    NF1_annotation = normalize_nf1(NF1),
    treatment = as.character(treatment),
    response_clean = normalize_response(response)
  ) |>
  dplyr::distinct()

if (anyDuplicated(annot_clean$sample_id)) {
  stop(
    "Historical response annotation contains duplicate sample/core IDs."
  )
}

# ------------------------------------------------------------------------------
# 3. Join response annotation to cached EGFR metrics
# ------------------------------------------------------------------------------

egfr_joined <- egfr_core |>
  dplyr::left_join(
    annot_clean,
    by = "sample_id"
  )

# Verify NF1 agreement.
nf1_mismatch <- egfr_joined |>
  dplyr::filter(
    !is.na(NF1_annotation),
    NF1_clean != NF1_annotation
  ) |>
  dplyr::distinct(
    sample_id,
    NF1_clean,
    NF1_annotation
  )

if (nrow(nf1_mismatch) > 0L) {
  stop(
    "NF1 mismatch between cached EGFR data and historical annotation."
  )
}

cat("\nResponse-coded tissue cores by NF1 / treatment / response:\n")

print(
  egfr_joined |>
    dplyr::distinct(
      sample_id,
      NF1_clean,
      treatment,
      response_clean
    ) |>
    dplyr::filter(
      !is.na(response_clean)
    ) |>
    dplyr::count(
      NF1_clean,
      treatment,
      response_clean,
      name = "n_cores"
    )
)

# ------------------------------------------------------------------------------
# 4. Pool EGFR+ cells by NF1 x response x meta-niche
# ------------------------------------------------------------------------------

response_dot <- egfr_joined |>
  dplyr::filter(
    !is.na(response_clean),
    !is.na(meta_niche_clean),
    n_EGFR_pos > 0
  ) |>
  dplyr::group_by(
    NF1_clean,
    response_clean,
    meta_niche_clean
  ) |>
  dplyr::summarise(
    n_EGFRpos = sum(
      n_EGFR_pos,
      na.rm = TRUE
    ),

    mean_EGFR = stats::weighted.mean(
      mean_EGFR_pos,
      w = n_EGFR_pos,
      na.rm = TRUE
    ),

    .groups = "drop"
  )

if (nrow(response_dot) == 0L) {
  stop("No response-coded EGFR+ rows were generated.")
}

# Public output contains no sample/core IDs.
utils::write.csv(
  response_dot,
  file.path(
    out_dir,
    "Figure5C_S6A_EGFR_response_dotplot_values_FINAL.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------------------------
# 5. Meta-niche and response ordering
# ------------------------------------------------------------------------------

MN_PLOT_LEVELS <- c(
  "MetaNiche_1",
  "MetaNiche_10",
  "MetaNiche_11",
  "MetaNiche_12",
  "MetaNiche_2",
  "MetaNiche_3",
  "MetaNiche_4",
  "MetaNiche_5",
  "MetaNiche_6",
  "MetaNiche_7",
  "MetaNiche_8",
  "MetaNiche_9"
)

MN_LABELS <- setNames(
  c(
    "1", "10", "11", "12",
    "2", "3", "4", "5",
    "6", "7", "8", "9"
  ),
  MN_PLOT_LEVELS
)

response_dot <- response_dot |>
  dplyr::mutate(
    response_clean = factor(
      response_clean,
      levels = c("NR", "R")
    ),
    meta_niche_clean = factor(
      meta_niche_clean,
      levels = MN_PLOT_LEVELS
    )
  )

# Use common scales for Figure 5C and S6A so the panels are directly comparable.
size_limits <- range(
  response_dot$n_EGFRpos,
  na.rm = TRUE
)

color_limits <- range(
  response_dot$mean_EGFR,
  na.rm = TRUE
)

# ------------------------------------------------------------------------------
# 6. Plot helper
# ------------------------------------------------------------------------------

make_response_plot <- function(
    nf1_group,
    title_text,
    output_prefix
) {

  z <- response_dot |>
    dplyr::filter(
      NF1_clean == nf1_group
    )

  if (nrow(z) == 0L) {
    stop("No response data for NF1 group: ", nf1_group)
  }

  cat(
    "\n",
    nf1_group,
    " pooled EGFR+ count range: ",
    paste(
      range(z$n_EGFRpos, na.rm = TRUE),
      collapse = " to "
    ),
    "\n",
    nf1_group,
    " pooled mean EGFR range: ",
    paste(
      range(z$mean_EGFR, na.rm = TRUE),
      collapse = " to "
    ),
    "\n",
    sep = ""
  )

  p <- ggplot2::ggplot(
    z,
    ggplot2::aes(
      x = response_clean,
      y = meta_niche_clean
    )
  ) +
    ggplot2::geom_point(
      ggplot2::aes(
        size = n_EGFRpos,
        color = mean_EGFR
      )
    ) +
    ggplot2::scale_y_discrete(
      labels = MN_LABELS,
      drop = FALSE
    ) +
    ggplot2::scale_x_discrete(
      labels = c(
        "NR" = "NR",
        "R" = "R"
      ),
      drop = FALSE
    ) +
    ggplot2::scale_size_continuous(
      range = c(1, 8),
      limits = size_limits,
      name = "EGFR+ cells"
    ) +
    ggplot2::scale_color_gradient(
      low = "grey80",
      high = "purple4",
      limits = color_limits,
      name = "Mean EGFR"
    ) +
    ggplot2::theme_classic(
      base_size = 13
    ) +
    ggplot2::labs(
      title = title_text,
      x = "Response",
      y = "Meta-Niche"
    ) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(
        face = "bold"
      ),
      axis.text.y = ggplot2::element_text(
        size = 11
      ),
      legend.position = "right"
    )

  print(p)

  ggplot2::ggsave(
    file.path(
      out_dir,
      paste0(output_prefix, ".png")
    ),
    p,
    width = 5.3,
    height = 6.3,
    dpi = 300
  )

  ggplot2::ggsave(
    file.path(
      out_dir,
      paste0(output_prefix, ".pdf")
    ),
    p,
    width = 5.3,
    height = 6.3,
    useDingbats = FALSE
  )

  invisible(p)
}

# ------------------------------------------------------------------------------
# 7. Figure 5C and Figure S6A
# ------------------------------------------------------------------------------

p_fig5c <- make_response_plot(
  nf1_group = "Mut",
  title_text = "NF1 Mut",
  output_prefix = "Figure5C_NF1Mut_EGFR_by_response_meta_niche_FINAL"
)

p_figS6a <- make_response_plot(
  nf1_group = "WT",
  title_text = "NF1 WT",
  output_prefix = "FigureS6A_NF1WT_EGFR_by_response_meta_niche_FINAL"
)

cat(
  "\nDONE — Figure 5C and Figure S6A generated from cached Figure 5B EGFR metrics.\n"
)
