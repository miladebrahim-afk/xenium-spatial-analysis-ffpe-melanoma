# GitHub release checklist

## Commit these

- `README.md`
- `.gitignore`
- `CITATION.cff`
- `R/`
- `config/`
- `scripts/`
- `docs/`
- `run_xenium_pipeline.R`
- `run_cellchat.R`
- `data/README.md`
- `data/reference/published_meta_niche_boundary_corrections.csv`
- `data/reference/figure5D_historical_pathway_modules.rds`
- small templates/README files under `data/`

## Do not commit these to public GitHub

- `data/processed/merged_obj_filtered.rds`
- `data/processed/xenium_list.rds`
- any other large raw/processed Xenium RDS/parquet files unless deposited through an approved data archive
- `private_provenance/patient_treatment_response_annotation.xlsx`
- old annotation backups containing legacy labels or original identifiers
- old/diagnostic/rescue scripts
- generated `results/` unless specific final outputs are intentionally selected for release

## Local-only files needed for a full rerun

Place the large starting objects in `data/processed/` and the private response spreadsheet in `private_provenance/` before the complete local manuscript rerun.
