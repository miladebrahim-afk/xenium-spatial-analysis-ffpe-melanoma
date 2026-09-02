# Main manuscript Xenium reproduction flow.
# Run from repository root. CellChat is intentionally separate (run_cellchat.R).

run_script <- function(path) {
  cat("\n============================================================\n")
  cat("RUNNING: ", path, "\n", sep = "")
  cat("============================================================\n")
  source(path, echo = FALSE, chdir = FALSE)
  invisible(gc())
}

run_script("scripts/00_setup/00_check_inputs_and_packages.R")

run_script("scripts/02_annotation/00_validate_starting_object.R")
run_script("scripts/02_annotation/01_find_cluster_markers.R")
run_script("scripts/02_annotation/02_cluster_annotation_map.R")
run_script("scripts/02_annotation/03_apply_final_annotations.R")
run_script("scripts/02_annotation/04_generate_umaps_and_dotplots.R")
run_script("scripts/02_annotation/05_celltype_composition.R")
run_script("scripts/02_annotation/06_map_annotations_to_cores.R")
run_script("scripts/02_annotation/07_spatial_focus_panels.R")

run_script("scripts/03_meta_niches/01_build_meta_niches.R")
run_script("scripts/03_meta_niches/02_apply_published_meta_niche_reference.R")
run_script("scripts/03_meta_niches/03_mn_composition_and_nf1_abundance.R")

# Spatial display of all 12 publication-reference MNs.
# Selected/highlighted MN-only maps are intentionally omitted.
run_script("scripts/03_meta_niches/04_spatial_meta_niche_maps.R")

run_script("scripts/05_spatial_distance/01_distance_analysis.R")
run_script("scripts/05_spatial_distance/02_distance_qc.R")
run_script("scripts/05_spatial_distance/03_cd8_abundance_and_spatial.R")
run_script("scripts/05_spatial_distance/04_cd8_core_pseudobulk_de.R")

run_script("scripts/06_antigen_presentation/01_time_signature_scores.R")
run_script("scripts/06_antigen_presentation/02_mhc_gene_expression.R")

run_script("scripts/07_egfr/01_metaniche_hallmark_gsea.R")

# Figure 5B: direct raw-RNA EGFR extraction. No SCT.
run_script("scripts/07_egfr/02_egfr_metaniche_dotplots.R")

# Figure 5C / S6A: validated paper-matching cached Figure 5B workflow.
response_candidates <- c(
  file.path(
    "private_provenance",
    "patient_treatment_response_annotation.xlsx"
  ),
  "patient_treatment_response_annotation.xlsx"
)

if (any(file.exists(response_candidates))) {
  run_script(
    "scripts/07_egfr/02b_egfr_response_dotplots_from_cached_fig5b.R"
  )
} else {
  cat(
    "\nSKIPPING Figure 5C/S6A: private response annotation is not present.\n"
  )
}

run_script("scripts/07_egfr/03_egfr_pathway_scores.R")
run_script("scripts/07_egfr/04_egfr_gene_correlations.R")

run_script("scripts/08_supplement/02_celltype_de_enrichment.R")
run_script("scripts/09_reproducibility/01_capture_session_info.R")

cat("\nMain Xenium pipeline complete. All generated outputs are under results/.\n")
