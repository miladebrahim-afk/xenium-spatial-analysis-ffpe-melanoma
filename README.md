# Xenium Spatial Analysis of Human FFPE Melanoma

This repository contains the reproducible computational workflow for spatial transcriptomic analysis of human formalin-fixed paraffin-embedded (FFPE) melanoma using the **10x Genomics Xenium** platform.

The study integrates cell-type annotation, spatial mapping, **meta-niche analysis of recurrent cellular neighborhoods**, differential gene expression, pathway enrichment, cell–cell communication, spatial distance analyses, and EGFR-focused analyses.

A central component of the workflow is the identification and characterization of **12 recurrent meta-niches (MNs)** derived from **40-µm cellular neighborhoods**. These meta-niches represent recurring spatial states within the melanoma tumor microenvironment and provide a framework for investigating cellular composition, **NF1 loss-associated spatial organization**, transcriptional programs, pathway activity, and EGFR-related biology.

## Overview

The analysis was developed to characterize the cellular, transcriptional, and spatial architecture of the melanoma tumor microenvironment using high-resolution Xenium spatial transcriptomics.

In addition to analyses at the level of individual cell populations, the workflow examines how cells organize into recurrent spatial communities. Local cellular neighborhoods are constructed around individual cells, summarized according to their cellular composition, and grouped into higher-order **meta-niches** representing recurring microenvironmental configurations.

Major components of the analysis include:

* Xenium data preprocessing and quality control
* Cell-type annotation
* Spatial visualization of annotated cell populations
* Construction of **40-µm cellular neighborhoods**
* Identification and characterization of **12 meta-niches**
* Meta-niche cellular composition analysis
* **NF1 loss-associated meta-niche abundance and prevalence analyses**
* Differential gene expression analysis
* Cell-type and meta-niche pathway enrichment analysis
* Gene set enrichment analysis (GSEA)
* Cell–cell communication analysis using CellChat
* Spatial distance and nearest-neighbor analyses
* Immune and tumor-microenvironment analyses
* EGFR expression and spatial analyses
* EGFR-associated meta-niche analyses
* Generation of manuscript and supplementary figures

## Analytical Framework

The overall workflow follows the general structure:

```text
Xenium data
    │
    ▼
Preprocessing and quality control
    │
    ▼
Cell-type annotation
    │
    ▼
Spatial mapping of cell populations
    │
    ▼
40-µm cellular neighborhood construction
    │
    ▼
Meta-niche identification
    │
    ▼
12 recurrent meta-niches
    │
    ├── Cellular composition
    ├── NF1 loss-associated abundance/prevalence
    ├── Differential expression
    ├── Pathway enrichment / GSEA
    └── EGFR-associated analyses
    │
    ▼
Cell–cell communication
Spatial distance / nearest-neighbor analyses
Additional tumor and immune analyses
    │
    ▼
Manuscript figures and supplementary analyses
```

## Meta-Niche Analysis

The meta-niche analysis is a central component of this study.

Local cellular neighborhoods are defined using a **40-µm spatial radius** and represented according to the composition of annotated cell populations within each neighborhood.

These local neighborhoods are subsequently clustered into **12 recurrent meta-niches**, representing distinct spatial configurations of tumor, immune, stromal, and other microenvironmental cell populations.

Downstream meta-niche analyses include:

* Meta-niche cellular composition
* Spatial distribution of meta-niches
* Meta-niche abundance
* Meta-niche prevalence
* **NF1 loss-associated differences in meta-niche representation**
* Meta-niche transcriptional programs
* Meta-niche pseudobulk analyses
* Pathway and gene-set enrichment
* EGFR expression across spatial niches
* Relationships between meta-niches and other tumor-microenvironment features

This framework enables analysis of the tissue not only at the level of individual cells, but also at the level of recurrent multicellular spatial ecosystems.

## Repository Structure

The repository is organized into modular scripts corresponding to the major stages of the analysis.

```text
.
├── README.md
├── LICENSE
├── .gitignore
├── scripts/
│   └── run_xenium_workflow.R
├── R/
├── figures/
├── results/
└── renv.lock
```

The exact directory structure will be finalized following validation of the complete reproducibility workflow.

## Running the Analysis

The complete Xenium analysis workflow can be initiated using:

```r
source("scripts/run_xenium_workflow.R")
```

Individual analysis modules are maintained as separate R scripts to facilitate reproducibility, troubleshooting, and reuse.

The main workflow coordinates the analysis stages in the appropriate order.

## Input Data

The Xenium spatial transcriptomics dataset used in this study is available through the **NCBI Gene Expression Omnibus (GEO)** under accession:

**GSE316387**

The corresponding **response table used in the study is provided as Supplementary Table S1 in the associated manuscript**.

The workflow requires the appropriate Xenium-derived input files and analysis metadata to be available locally before execution.

Large raw data files, intermediate R objects, and processed Xenium datasets are not duplicated within this GitHub repository.

Users reproducing the workflow should obtain the Xenium data from GEO and use the accompanying study metadata as described in the associated manuscript and analysis scripts.

Detailed input-file and directory requirements will be finalized following validation of the complete workflow.

## Software Environment

Analyses were performed in **R** using packages for:

* Single-cell and spatial transcriptomic analysis
* Data manipulation
* Statistical analysis
* Data visualization
* Differential gene expression
* Gene-set and pathway enrichment
* Cell–cell communication
* Spatial neighborhood analysis

The computational environment will be recorded using **renv**.

Exact package versions required to reproduce the analysis will be stored in:

```text
renv.lock
```

After cloning the repository, the R environment can be restored using:

```r
install.packages("renv")
renv::restore()
```

The R session used for the final validated workflow will also be documented using:

```r
sessionInfo()
```

## Reproducibility

To reproduce the analysis:

1. Clone or download this repository.
2. Install the appropriate version of R.
3. Open the project in R or RStudio.
4. Restore the package environment using `renv`.
5. Obtain the Xenium dataset from **GEO accession GSE316387**.
6. Provide the required study metadata, including the response information described in **Supplementary Table S1** of the associated manuscript.
7. Confirm the expected input paths in the workflow configuration.
8. Run:

```r
source("scripts/run_xenium_workflow.R")
```

9. Review the generated results and figures.

## Analysis Outputs

The workflow generates outputs associated with the major analyses in the study, including:

* Cell-type annotation summaries
* Spatial maps of annotated cell populations
* Cellular-neighborhood measurements
* Meta-niche assignments
* Meta-niche composition summaries
* Meta-niche abundance and prevalence statistics
* NF1 loss-associated spatial analyses
* Differential gene expression results
* Pathway enrichment and GSEA results
* CellChat communication results
* Spatial distance and nearest-neighbor statistics
* EGFR-related analyses
* Manuscript-quality figures
* Supplementary analyses and figures

Generated outputs may be written to the `results/` and `figures/` directories depending on the individual analysis module.

## Figures

The repository contains analysis code used to generate the major and supplementary analyses associated with the study.

These include analyses related to:

* Spatial organization of melanoma and tumor-microenvironment cell populations
* Cellular neighborhoods and meta-niches
* **NF1 loss-associated spatial differences**
* Meta-niche composition and prevalence
* Cell–cell communication
* Differential gene expression
* Pathway enrichment and GSEA
* Immune-cell spatial organization
* EGFR expression and spatial context
* EGFR-associated meta-niche analyses
* Nearest-neighbor and distance-based analyses

A detailed mapping between individual scripts and manuscript figure panels will be included following final workflow validation.

## Data Availability

The Xenium spatial transcriptomics data associated with this study are available through the **NCBI Gene Expression Omnibus (GEO)** under accession **GSE316387**.

The study response table is provided as **Supplementary Table S1** in the associated manuscript.

Raw Xenium data are not duplicated within this GitHub repository. This repository contains the analysis code and supporting computational workflow required to reproduce the reported analyses using the corresponding public dataset and study metadata.

Additional data availability information should be referenced from the associated manuscript and GEO record.

## Code Availability

The analysis scripts required for the computational workflow are provided in this repository.

Large datasets, intermediate R objects, locally installed packages, and temporary computational files are intentionally excluded from version control.

## Citation

If you use this code, analysis framework, or meta-niche workflow, please cite the associated publication.

Full citation information will be added following publication.

## License

The analysis code in this repository is released under the **MIT License**.

See the `LICENSE` file for details.

## Contact

Questions regarding the analysis workflow or reproducibility of the computational analyses should be directed to the corresponding authors of the associated study.
