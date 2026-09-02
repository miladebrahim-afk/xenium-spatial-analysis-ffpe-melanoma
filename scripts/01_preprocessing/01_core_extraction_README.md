# Manual TMA core definition in Xenium Explorer

Core isolation was a **manual preprocessing step**, not a scripted analysis step.

The workflow used for the study was:

1. Open each Xenium region in **Xenium Explorer**.
2. Use the Xenium Explorer annotation tools to manually outline each TMA core.
3. Export the coordinates of each annotated core as a CSV file.
4. Use those exported coordinate boundaries to define the individual cores for downstream Seurat analysis.
5. Save each resulting core as an individual Seurat RDS object for merging and downstream processing.

The manuscript additionally describes reconstruction of the core-specific spatial objects using Seurat spatial classes/functions (`CreateSegmentation()`, `Overlay()`, `CreateCentroids()`, `CreateFOV()`, and `CreateMolecules()`). The important provenance point is that the **core boundaries themselves were selected manually in Xenium Explorer** and exported as coordinate CSV files; they were not inferred computationally.

The public executable pipeline resumes from the resulting individual per-core Seurat objects placed in:

```text
data/per_core/
```

Each per-core object should contain, at minimum:

- an RNA count assay;
- `patient_id` and `NF1` metadata;
- a stable core identifier (`orig.ident` or equivalent);
- the spatial FOV / tissue coordinates required for `GetTissueCoordinates()`.

## Coordinate CSV files

If the de-identified Xenium Explorer core-coordinate CSV files are available and permitted for public sharing, they can be deposited with the public data record or placed under a documented `data/core_coordinates/` directory. They are useful provenance for reproducing the manual core boundaries even though the boundary-drawing step itself cannot be automated retrospectively.

The region-level transcript conversion utility is preserved as `00_convert_transcripts_parquet.R`.
