# ==============================================================================
# 03_meta_niche_cellchat_figures.R
#
# Generate manuscript-style CellChat circle plots from the saved niche-specific
# CellChat models produced by 02_meta_niche_cellchat_nf1.R.
#
# Expected models (missing non-MN3 models are skipped without crashing):
#   results/cellchat/niche_models/CellChat_MetaNiche_2_Mut.rds
#   results/cellchat/niche_models/CellChat_MetaNiche_2_WT.rds
#   results/cellchat/niche_models/CellChat_MetaNiche_3_Mut.rds
#   results/cellchat/niche_models/CellChat_MetaNiche_3_WT.rds
#   results/cellchat/niche_models/CellChat_MetaNiche_6_Mut.rds
#   results/cellchat/niche_models/CellChat_MetaNiche_6_WT.rds
#
# Outputs:
#   results/cellchat/niche_figures/circle_plots/
#   results/cellchat/niche_figures/mn3_differential/
#
#
# NOTE
# ----
# If an optional model such as MN6 WT is absent (for example because no
# interactions were inferred and no RDS was saved), the script records the
# missing model in QC and continues. MN3 WT and MN3 Mut remain required because
# they are needed for the differential Figure 2G analyses.
#
# Figure correspondence:
#   Fig. 2D: MN2 NF1Mut vs NF1WT
#   Fig. 2E: MN6 NF1Mut vs NF1WT
#   Fig. 2F: MN3 NF1Mut vs NF1WT
#   Fig. 2G: MN3 WT-vs-Mut differential interaction number/strength
# ==============================================================================

suppressPackageStartupMessages({
  library(CellChat)
  library(patchwork)
  library(ggplot2)
})

# ------------------------------------------------------------------------------
# 0. Paths
# ------------------------------------------------------------------------------

model_dir <- file.path(
  "results",
  "cellchat",
  "niche_models"
)

figure_dir <- file.path(
  "results",
  "cellchat",
  "niche_figures"
)

circle_dir <- file.path(
  figure_dir,
  "circle_plots"
)

mn3_diff_dir <- file.path(
  figure_dir,
  "mn3_differential"
)

dir.create(
  circle_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  mn3_diff_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

MNS <- c(
  "MetaNiche_2",
  "MetaNiche_3",
  "MetaNiche_6"
)

GENOTYPES <- c(
  "Mut",
  "WT"
)

# ------------------------------------------------------------------------------
# 1. Helper
# ------------------------------------------------------------------------------

model_path <- function(mn, gt) {
  file.path(
    model_dir,
    paste0(
      "CellChat_",
      mn,
      "_",
      gt,
      ".rds"
    )
  )
}

load_cellchat_model <- function(mn, gt, required = FALSE) {

  f <- model_path(
    mn,
    gt
  )

  if (!file.exists(f)) {
    msg <- paste0(
      "Missing CellChat model: ",
      f
    )

    if (required) {
      stop(msg)
    } else {
      message("SKIP: ", msg)
      return(NULL)
    }
  }

  obj <- readRDS(
    f
  )

  # Allows compatibility with objects created under an older CellChat release.
  obj <- updateCellChat(
    obj
  )

  obj
}

# ------------------------------------------------------------------------------
# 2. Standard circle plots: count + interaction strength
# ------------------------------------------------------------------------------

cat("Generating CellChat circle plots...\n")

qc <- list()
idx <- 0L

for (mn in MNS) {

  for (gt in GENOTYPES) {

    idx <- idx + 1L

    nm <- paste0(
      mn,
      "_",
      gt
    )

    cat(
      "  ",
      nm,
      " ... ",
      sep = ""
    )

    cellchat <- load_cellchat_model(
      mn,
      gt,
      required = FALSE
    )

    if (is.null(cellchat)) {

      qc[[idx]] <- data.frame(
        comparison = nm,
        model_found = FALSE,
        total_interaction_count = NA_real_,
        total_interaction_weight = NA_real_,
        zero_interaction_network = NA,
        count_plot_ok = FALSE,
        weight_plot_ok = FALSE,
        stringsAsFactors = FALSE
      )

      cat("SKIPPED (model not found)\n")
      next
    }

    if (is.null(cellchat@net$count)) {
      stop(
        "CellChat model has no @net$count matrix: ",
        nm
      )
    }

    zero_network <- sum(cellchat@net$count, na.rm = TRUE) == 0

    if (zero_network) {
      cat("node-only / zero-interaction network ... ")
    }

    vertex_weight <- as.numeric(
      table(
        cellchat@idents
      )
    )

    # --------------------------------------------------------------------------
    # Number of interactions
    # --------------------------------------------------------------------------

    count_file <- file.path(
      circle_dir,
      paste0(
        "Circle_count_",
        nm,
        ".pdf"
      )
    )

    pdf(
      count_file,
      width = 8,
      height = 8,
      useDingbats = FALSE
    )

    netVisual_circle(
      cellchat@net$count,
      vertex.weight = vertex_weight,
      weight.scale = TRUE,
      label.edge = FALSE,
      title.name = paste0(
        mn,
        " | NF1 ",
        gt,
        " | Number of interactions"
      )
    )

    dev.off()

    # --------------------------------------------------------------------------
    # Interaction strength
    # --------------------------------------------------------------------------

    weight_file <- file.path(
      circle_dir,
      paste0(
        "Circle_weight_",
        nm,
        ".pdf"
      )
    )

    pdf(
      weight_file,
      width = 8,
      height = 8,
      useDingbats = FALSE
    )

    netVisual_circle(
      cellchat@net$weight,
      vertex.weight = vertex_weight,
      weight.scale = TRUE,
      label.edge = FALSE,
      title.name = paste0(
        mn,
        " | NF1 ",
        gt,
        " | Interaction strength"
      )
    )

    dev.off()

    qc[[idx]] <- data.frame(
      comparison = nm,
      model_found = TRUE,
      total_interaction_count = sum(
        cellchat@net$count
      ),
      total_interaction_weight = sum(
        cellchat@net$weight
      ),
      zero_interaction_network = zero_network,
      count_plot_ok = file.exists(
        count_file
      ),
      weight_plot_ok = file.exists(
        weight_file
      ),
      stringsAsFactors = FALSE
    )

    cat("OK\n")

    rm(
      cellchat,
      vertex_weight
    )
    gc(
      verbose = FALSE
    )
  }
}

qc_df <- do.call(
  rbind,
  qc
)

write.csv(
  qc_df,
  file.path(
    figure_dir,
    "cellchat_circle_plot_QC.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------------------------
# 3. MN3 differential interaction plots: WT vs Mut
#
# Archived workflow merged WT first and Mut second, then used
# netVisual_diffInteraction(). In CellChat's comparison convention, the
# displayed edge differences correspond to group 2 minus group 1 internally;
# therefore titles are kept descriptive rather than manually re-signing edges.
# ------------------------------------------------------------------------------

cat("\nGenerating MN3 differential CellChat plots...\n")

cellchat_mn3_wt <- load_cellchat_model(
  "MetaNiche_3",
  "WT",
  required = TRUE
)

cellchat_mn3_mut <- load_cellchat_model(
  "MetaNiche_3",
  "Mut",
  required = TRUE
)

object_list_mn3 <- list(
  WT = cellchat_mn3_wt,
  Mut = cellchat_mn3_mut
)

cellchat_mn3_wt_mut <- mergeCellChat(
  object_list_mn3,
  add.names = names(
    object_list_mn3
  )
)

# Combined differential count + strength panel.
diff_circle_file <- file.path(
  mn3_diff_dir,
  "MetaNiche_3_CellChat_WT_vs_Mut_diff_circle_count_strength.pdf"
)

pdf(
  diff_circle_file,
  width = 14,
  height = 7,
  useDingbats = FALSE
)

par(
  mfrow = c(1, 2),
  xpd = TRUE
)

netVisual_diffInteraction(
  cellchat_mn3_wt_mut,
  weight.scale = TRUE,
  measure = "count",
  title.name = "MN3 WT vs Mut: differential number of interactions"
)

netVisual_diffInteraction(
  cellchat_mn3_wt_mut,
  weight.scale = TRUE,
  measure = "weight",
  title.name = "MN3 WT vs Mut: differential interaction strength"
)

dev.off()

# ------------------------------------------------------------------------------
# 4. MN3 overall interaction-number / strength comparison
# ------------------------------------------------------------------------------

gg_count <- compareInteractions(
  cellchat_mn3_wt_mut,
  show.legend = FALSE,
  group = c(1, 2),
  measure = "count"
)

gg_weight <- compareInteractions(
  cellchat_mn3_wt_mut,
  show.legend = FALSE,
  group = c(1, 2),
  measure = "weight"
)

p_compare <- gg_count + gg_weight

ggsave(
  filename = file.path(
    mn3_diff_dir,
    "MetaNiche_3_CellChat_WT_vs_Mut_compareInteractions.pdf"
  ),
  plot = p_compare,
  width = 8,
  height = 4,
  useDingbats = FALSE
)

# ------------------------------------------------------------------------------
# 5. MN3 differential heatmaps
# ------------------------------------------------------------------------------

ht_count <- netVisual_heatmap(
  cellchat_mn3_wt_mut,
  measure = "count",
  title.name = "MN3 WT vs Mut: differential interaction number"
)

ht_weight <- netVisual_heatmap(
  cellchat_mn3_wt_mut,
  measure = "weight",
  title.name = "MN3 WT vs Mut: differential interaction strength"
)

pdf(
  file.path(
    mn3_diff_dir,
    "MetaNiche_3_CellChat_WT_vs_Mut_diff_heatmaps.pdf"
  ),
  width = 14,
  height = 7,
  useDingbats = FALSE
)

print(
  ht_count + ht_weight
)

dev.off()

# ------------------------------------------------------------------------------
# 6. Final QC
# ------------------------------------------------------------------------------

cat("\nCircle plot QC:\n")
print(
  qc_df,
  row.names = FALSE
)

cat(
  "\nCircle plots directory:\n  ",
  circle_dir,
  "\n",
  sep = ""
)

cat(
  "\nMN3 differential directory:\n  ",
  mn3_diff_dir,
  "\n",
  sep = ""
)

cat(
  "\nDONE: niche-specific CellChat figures complete.\n"
)
