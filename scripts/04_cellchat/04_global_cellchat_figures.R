# ==============================================================================
# 04_global_cellchat_figures.R
#
# Generate global NF1Mut vs NF1WT CellChat panels for Figure S3.
#
# Inputs
# ------
# results/cellchat/global_models/CellChat_global_Mut.rds
# results/cellchat/global_models/CellChat_global_WT.rds
#
# Outputs
# -------
# results/cellchat/global_figures/
#   FigureS3A_global_diff_interaction_strength.pdf
#   FigureS3B_global_diff_interaction_number.pdf
#   FigureS3C_global_rankNet_information_flow.pdf
#   FigureS3D_global_signaling_patterns_Mut.pdf
#   FigureS3E_global_signaling_patterns_WT.pdf
#   FigureS3DE_global_signaling_patterns_side_by_side.pdf
#   global_pathway_presence.csv
#   global_cellchat_figure_QC.csv
#
# Figure logic follows the archived CellChat workflow:
#   A = differential interaction strength
#   B = differential interaction number
#   C = pathway information-flow comparison
#   D/E = overall signaling-role heatmaps in NF1Mut / NF1WT
#
# For differential heatmaps, models are merged as WT first, Mut second.
# CellChat therefore displays Mut - WT: red = stronger/higher in NF1Mut,
# blue = stronger/higher in NF1WT, matching the final Figure S3 legend.
# ==============================================================================

suppressPackageStartupMessages({
  library(CellChat)
  library(ComplexHeatmap)
  library(circlize)
  library(ggplot2)
  library(patchwork)
  library(grid)
})

model_dir <- file.path("results", "cellchat", "global_models")
out_dir <- file.path("results", "cellchat", "global_figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

mut_file <- file.path(model_dir, "CellChat_global_Mut.rds")
wt_file  <- file.path(model_dir, "CellChat_global_WT.rds")

if (!file.exists(mut_file)) stop("Missing: ", mut_file)
if (!file.exists(wt_file)) stop("Missing: ", wt_file)

NF1_COLORS <- c(
  WT = "#00A6A6",
  Mut = "#E64B35"
)

cat("Loading global CellChat models...\n")

cellchat_mut <- updateCellChat(readRDS(mut_file))
cellchat_wt  <- updateCellChat(readRDS(wt_file))

cat(
  "NF1Mut: ", sum(cellchat_mut@net$count),
  " interactions; total weight ",
  signif(sum(cellchat_mut@net$weight), 6), "\n", sep = ""
)

cat(
  "NF1WT:  ", sum(cellchat_wt@net$count),
  " interactions; total weight ",
  signif(sum(cellchat_wt@net$weight), 6), "\n", sep = ""
)

mut_groups <- levels(cellchat_mut@idents)
wt_groups  <- levels(cellchat_wt@idents)

if (!setequal(mut_groups, wt_groups)) {
  stop("NF1Mut and NF1WT objects do not contain the same cell-type groups.")
}

mut_pathways <- cellchat_mut@netP$pathways
wt_pathways  <- cellchat_wt@netP$pathways
pathway_union <- union(mut_pathways, wt_pathways)

pathway_qc <- data.frame(
  pathway = pathway_union,
  in_NF1Mut = pathway_union %in% mut_pathways,
  in_NF1WT = pathway_union %in% wt_pathways,
  stringsAsFactors = FALSE
)

write.csv(
  pathway_qc,
  file.path(out_dir, "global_pathway_presence.csv"),
  row.names = FALSE
)

cat("\nSignaling-pathway QC:\n")
cat("  NF1Mut pathways: ", length(mut_pathways), "\n", sep = "")
cat("  NF1WT pathways:  ", length(wt_pathways), "\n", sep = "")
cat("  Union pathways:  ", length(pathway_union), "\n", sep = "")

cat(
  "  NF1Mut-only: ",
  paste(sort(setdiff(mut_pathways, wt_pathways)), collapse = ", "),
  "\n", sep = ""
)

cat(
  "  NF1WT-only:  ",
  paste(sort(setdiff(wt_pathways, mut_pathways)), collapse = ", "),
  "\n", sep = ""
)

# WT first, Mut second => differential plots are Mut - WT.
object_list_wt_mut <- list(
  WT = cellchat_wt,
  Mut = cellchat_mut
)

cellchat_wt_mut <- mergeCellChat(
  object_list_wt_mut,
  add.names = names(object_list_wt_mut)
)

# Figure S3A
cat("\nGenerating Figure S3A...\n")

ht_strength <- netVisual_heatmap(
  cellchat_wt_mut,
  measure = "weight",
  title.name = "Differential NF1Mut vs NF1WT interaction strength"
)

pdf(
  file.path(out_dir, "FigureS3A_global_diff_interaction_strength.pdf"),
  width = 9,
  height = 9,
  useDingbats = FALSE
)
ComplexHeatmap::draw(ht_strength)
dev.off()

# Figure S3B
cat("Generating Figure S3B...\n")

ht_count <- netVisual_heatmap(
  cellchat_wt_mut,
  measure = "count",
  title.name = "Differential NF1Mut vs NF1WT interaction number"
)

pdf(
  file.path(out_dir, "FigureS3B_global_diff_interaction_number.pdf"),
  width = 9,
  height = 9,
  useDingbats = FALSE
)
ComplexHeatmap::draw(ht_count)
dev.off()

pdf(
  file.path(out_dir, "FigureS3AB_global_diff_heatmaps_side_by_side.pdf"),
  width = 18,
  height = 9,
  useDingbats = FALSE
)
ComplexHeatmap::draw(
  ht_strength + ht_count,
  ht_gap = grid::unit(0.6, "cm")
)
dev.off()

# Figure S3C
cat("Generating Figure S3C...\n")

p_rank_stacked <- rankNet(
  cellchat_wt_mut,
  mode = "comparison",
  comparison = c(1, 2),
  stacked = TRUE,
  do.stat = TRUE,
  color.use = unname(NF1_COLORS[c("WT", "Mut")])
)

ggsave(
  file.path(out_dir, "FigureS3C_global_rankNet_information_flow.pdf"),
  p_rank_stacked,
  width = 8,
  height = 11,
  useDingbats = FALSE
)

ggsave(
  file.path(out_dir, "FigureS3C_global_rankNet_information_flow.png"),
  p_rank_stacked,
  width = 8,
  height = 11,
  dpi = 300
)

p_rank_unstacked <- rankNet(
  cellchat_wt_mut,
  mode = "comparison",
  comparison = c(1, 2),
  stacked = FALSE,
  do.stat = TRUE,
  color.use = unname(NF1_COLORS[c("WT", "Mut")])
)

ggsave(
  file.path(out_dir, "QC_global_rankNet_unstacked.pdf"),
  p_rank_unstacked,
  width = 9,
  height = 11,
  useDingbats = FALSE
)

# Figure S3D/E
cat("Generating Figure S3D/E...\n")

ht_mut_all <- netAnalysis_signalingRole_heatmap(
  cellchat_mut,
  pattern = "all",
  signaling = pathway_union,
  title = "NF1Mut",
  width = 6,
  height = 12,
  color.heatmap = "OrRd"
)

ht_wt_all <- netAnalysis_signalingRole_heatmap(
  cellchat_wt,
  pattern = "all",
  signaling = pathway_union,
  title = "NF1WT",
  width = 6,
  height = 12,
  color.heatmap = "OrRd"
)

pdf(
  file.path(out_dir, "FigureS3D_global_signaling_patterns_Mut.pdf"),
  width = 8,
  height = 13,
  useDingbats = FALSE
)
ComplexHeatmap::draw(ht_mut_all)
dev.off()

pdf(
  file.path(out_dir, "FigureS3E_global_signaling_patterns_WT.pdf"),
  width = 8,
  height = 13,
  useDingbats = FALSE
)
ComplexHeatmap::draw(ht_wt_all)
dev.off()

pdf(
  file.path(out_dir, "FigureS3DE_global_signaling_patterns_side_by_side.pdf"),
  width = 16,
  height = 13,
  useDingbats = FALSE
)
ComplexHeatmap::draw(
  ht_mut_all + ht_wt_all,
  ht_gap = grid::unit(0.6, "cm")
)
dev.off()

# Transparent QC: outgoing/incoming role views from the archived workflow.
cat("Generating outgoing/incoming QC heatmaps...\n")

for (pattern_use in c("outgoing", "incoming")) {

  color_use <- if (pattern_use == "incoming") "GnBu" else "YlOrRd"

  ht_mut <- netAnalysis_signalingRole_heatmap(
    cellchat_mut,
    pattern = pattern_use,
    signaling = pathway_union,
    title = paste0("NF1Mut - ", pattern_use),
    width = 6,
    height = 12,
    color.heatmap = color_use
  )

  ht_wt <- netAnalysis_signalingRole_heatmap(
    cellchat_wt,
    pattern = pattern_use,
    signaling = pathway_union,
    title = paste0("NF1WT - ", pattern_use),
    width = 6,
    height = 12,
    color.heatmap = color_use
  )

  pdf(
    file.path(
      out_dir,
      paste0("QC_global_signaling_role_", pattern_use, "_side_by_side.pdf")
    ),
    width = 16,
    height = 13,
    useDingbats = FALSE
  )

  ComplexHeatmap::draw(
    ht_mut + ht_wt,
    ht_gap = grid::unit(0.6, "cm")
  )

  dev.off()

  rm(ht_mut, ht_wt)
}

expected_files <- c(
  "FigureS3A_global_diff_interaction_strength.pdf",
  "FigureS3B_global_diff_interaction_number.pdf",
  "FigureS3C_global_rankNet_information_flow.pdf",
  "FigureS3D_global_signaling_patterns_Mut.pdf",
  "FigureS3E_global_signaling_patterns_WT.pdf",
  "FigureS3DE_global_signaling_patterns_side_by_side.pdf"
)

qc <- data.frame(
  output = expected_files,
  exists = file.exists(file.path(out_dir, expected_files)),
  stringsAsFactors = FALSE
)

write.csv(
  qc,
  file.path(out_dir, "global_cellchat_figure_QC.csv"),
  row.names = FALSE
)

cat("\nGlobal CellChat figure QC:\n")
print(qc, row.names = FALSE)

if (!all(qc$exists)) {
  warning("One or more Figure S3 outputs were not generated.")
}

cat("\nDONE: Figure S3 global CellChat panels complete.\n")
