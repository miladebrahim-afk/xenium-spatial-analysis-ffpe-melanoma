# Public reference files required for exact manuscript reproduction

Place these two small publication-reference files here and track them in GitHub:

1. `published_meta_niche_boundary_corrections.csv`
   - 140 documented boundary-cell corrections that map the independent de novo k-means reconstruction to the archived final manuscript meta-niche assignment.
   - Contains no patient/specimen IDs in the public version.

2. `figure5D_historical_pathway_modules.rds`
   - Frozen 50 Hallmark modules recovered from the historical Figure 5D workspace.
   - Required because regenerating leading-edge modules can change pathway scores.

These files are intentionally treated as reference/provenance inputs, not as generated `results/` files.
