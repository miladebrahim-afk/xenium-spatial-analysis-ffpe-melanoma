# Reproducibility notes

- Final manuscript figures/results are the reference target when historical exploratory scripts contain conflicting comments.
- The public pipeline preserves a de novo meta-niche reconstruction and a transparent publication-reference correction step; the de novo result is not overwritten.
- The 39 analyzed Xenium tissue cores/specimens are treated as independent samples in the reconstructed core-level analyses.
- `patient_id` functions operationally as the core/specimen identifier in the saved Xenium objects and should not be described as 39 biological patients.
- CellChat is run separately because of memory and package-version sensitivity.
- Large RDS input objects and private response metadata are not intended for public Git tracking.
- All generated outputs are under `results/`.
