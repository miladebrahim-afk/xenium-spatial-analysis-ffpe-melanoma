# Figure-to-code map

| Manuscript panel | Production script |
|---|---|
| Fig. 1B/C/E and annotation QC | `scripts/02_annotation/01_find_cluster_markers.R`, `04_generate_umaps_and_dotplots.R` |
| Fig. 1D/F spatial cell-type maps | `scripts/02_annotation/06_map_annotations_to_cores.R`, `07_spatial_focus_panels.R` |
| Fig. 1G-I / Fig. 2A-C meta-niches | `scripts/03_meta_niches/01_build_meta_niches.R`, `02_apply_published_meta_niche_reference.R`, `03_mn_composition_and_nf1_abundance.R`, `04_spatial_meta_niche_maps.R` |
| Fig. 2D-G / S2 CellChat | `scripts/04_cellchat/02_meta_niche_cellchat_nf1.R`, `03_meta_niche_cellchat_figures.R` |
| Fig. S3 global CellChat | `scripts/04_cellchat/01_global_cellchat_nf1.R`, `04_global_cellchat_figures.R` |
| Fig. 3A | `scripts/05_spatial_distance/01_distance_analysis.R` |
| Fig. S4A | `scripts/05_spatial_distance/02_distance_qc.R` |
| Fig. 3B-E | `scripts/05_spatial_distance/03_cd8_abundance_and_spatial.R` |
| Fig. 3F-G | `scripts/05_spatial_distance/04_cd8_core_pseudobulk_de.R` |
| Fig. 4A | `scripts/06_antigen_presentation/01_time_signature_scores.R` |
| Fig. 4E | `scripts/06_antigen_presentation/02_mhc_gene_expression.R` |
| Fig. 5A | `scripts/07_egfr/01_metaniche_hallmark_gsea.R` |
| Fig. 5B/C and S6A | `scripts/07_egfr/02_egfr_metaniche_dotplots.R` |
| Fig. 5D | `scripts/07_egfr/03_egfr_pathway_scores.R` |
| Fig. 5F | `scripts/07_egfr/04_egfr_gene_correlations.R` |
| Fig. S5A-D | `scripts/08_supplement/02_celltype_de_enrichment.R` |

The cached `02b_egfr_response_dotplots_from_cached_fig5b_OPTIONAL_FAST.R` is retained only as a fast local replot option; it is not part of the primary ordered runner.
