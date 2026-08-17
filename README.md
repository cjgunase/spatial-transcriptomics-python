# Spatial transcriptomics in Python

A learning-first portfolio project that answers a simple question:

> Do gene-expression patterns and spatial relationships recover known anatomical
> regions in a mouse brain Visium section?

The project uses a public, preprocessed 10x Genomics Visium dataset distributed
through Squidpy. It demonstrates how to inspect an `AnnData` object, interpret
spot-level quality metrics, visualize expression in tissue coordinates, build a
spatial-neighbor graph, quantify neighborhood enrichment, and rank spatially
autocorrelated genes with Moran's I.

## Start here

Begin with
[`notebooks/00_spatial_transcriptomics_pipeline.ipynb`](notebooks/00_spatial_transcriptomics_pipeline.ipynb)
for a no-code, diagram-based explanation of the complete path from tissue to a
biological conclusion. Then work through
[`notebooks/01_visium_spatial_analysis.ipynb`](notebooks/01_visium_spatial_analysis.ipynb)
for the first Python analysis. Both lessons include interpretation prompts,
exercises, and short knowledge checks.

## What this demonstrates

- Understanding of spots, genes, tissue coordinates, and histology images
- Awareness that a Visium spot can contain multiple cells
- Expression QC and the distinction between QC and biological signal
- Spatial graphs and neighborhood-enrichment analysis
- Moran's I as a measure of spatial autocorrelation
- Reproducible Python analysis with Scanpy, AnnData, and Squidpy
- Cautious interpretation: clusters are not automatically cell types

## Setup

The first run downloads approximately 314 MB of example data.

```bash
conda env create -f environment.yml
conda activate spatial-transcriptomics-demo
jupyter lab
```

Open the notebook and choose **Kernel → Restart Kernel and Run All Cells**.

## Project layout

```text
.
├── notebooks/       # Narrative, executable lessons
├── src/             # Small reusable validation helpers
├── figures/         # Optional exported figures (generated locally)
├── results/         # Optional result tables (generated locally)
└── environment.yml  # Reproducible Conda environment
```

## Learning path

1. **Notebook 00:** the entire experimental and computational pipeline.
2. **Notebook 01:** one processed Visium section—data structure, QC, spatial graphs, Moran's I.
3. **Next:** start from raw counts and justify filtering and normalization.
4. **Then:** annotate regions using marker genes and differential expression.
5. **Later:** compare samples, integrate a single-cell reference, and turn stable
   steps into a workflow.

This intentionally starts as a notebook, not a production pipeline. A production
workflow becomes worthwhile after the analysis choices are understood and need to
be repeated across samples.

## Moving from Visium to Xenium

Moving to Xenium is not simply an increase in spatial resolution. It changes the
measurement technology, the primary data model, the main quality risks, and the
biological claims that can be supported.

### Central difference

| Visium | Xenium |
|---|---|
| Sequencing-based capture | Imaging-based in situ detection |
| Measures capture spots | Detects individual RNA molecules |
| A spot can contain several cells | Molecules are assigned to segmented cells |
| Broad or whole-transcriptome coverage | Predefined targeted gene panel |
| Primary matrix is spots × genes | Primary data are transcript coordinates; cells × genes is derived |
| H&E image plus spot coordinates | Morphology, transcript points, and cell/nucleus boundaries |
| Main concern is mixed-cell spots | Main concerns are segmentation and transcript assignment |
| Well suited to discovery | Well suited to cellular localization and validation |

In short:

```text
Visium: many genes, lower biological resolution
Xenium: fewer selected genes, much higher spatial resolution
```

Xenium Prime panels can be large, but panel composition still determines which
genes, cell types, and states can be detected.

### What transfers from Visium

The following skills remain useful:

- Expression matrices and AnnData tables
- Quality control and normalization
- PCA, clustering, and marker interpretation
- Spatial-neighbor graphs
- Neighborhood enrichment and Moran's I
- Distinguishing spatial association from biological mechanism

The biological unit changes:

```text
Visium observation       = spot
Xenium observation       = segmented cell
Xenium primary detection = decoded transcript molecule
```

### Xenium output files

Important files include:

| File | Contents |
|---|---|
| `transcripts.parquet` | One record per decoded transcript, with gene, coordinates, quality, and cell assignment |
| `cells.parquet` | Cell-level metadata and summaries |
| `cell_feature_matrix.h5` | Derived cell-by-gene count matrix |
| `cell_boundaries.parquet` | Simplified cell-boundary polygons |
| `nucleus_boundaries.parquet` | Simplified nuclear-boundary polygons |
| `cells.zarr.zip` | Segmentation masks used for transcript assignment |
| `morphology_focus/` | High-resolution morphology images |
| `gene_panel.json` | Targeted genes, controls, and panel-design information |

A transcript record can contain its gene identity, x/y/z position, quality score,
assigned cell, nuclear overlap, and distance from a nucleus. The cell-by-gene
matrix is therefore downstream of transcript decoding, segmentation, and molecule
assignment.

### The new central problem: segmentation

Visium spot geometry is predefined. Xenium must determine where nuclei and cells
are and which transcripts belong to each cell. Errors can produce:

- Two cells merged into one
- One cell divided into several cells
- Transcripts assigned to the wrong neighboring cell
- Cytoplasmic transcripts left unassigned
- Implausibly large or small cells
- False co-expression of markers from different cell types

Do not begin by trusting the cell-feature matrix alone. First overlay and inspect:

```text
morphology image
+ nucleus boundaries
+ cell boundaries
+ transcript locations
```

Review multiple tissue regions, boundaries, densities, and cell morphologies.

### How quality control changes

Typical Visium QC includes UMIs per spot, detected genes per spot, mitochondrial
fraction, tissue overlap, and spatial QC patterns.

Typical Xenium QC adds:

- Transcripts and detected genes per cell
- Cell and nucleus area
- Cell-to-nucleus area ratio
- Fraction of transcripts assigned to cells
- Nuclear versus cytoplasmic transcript distribution
- Negative-control and blank-codeword signals
- Transcript Q-scores
- Segmentation quality
- Spatial concentration of low-quality cells
- Field-of-view or image-boundary artifacts

An unusually high-count Xenium cell may be biologically active, or it may be two
cells incorrectly merged by segmentation.

### Normalization, clustering, and annotation

A Xenium cell-by-gene matrix can still be normalized, reduced with PCA, and
clustered. However:

- The panel is intentionally selected, so automatic highly-variable-gene selection
  can discard useful markers.
- Clusters can only be resolved if the panel includes discriminating genes.
- A missing gene may mean it was not targeted, not biologically absent.
- A broad cell-type panel may not distinguish subtle cellular states.
- Cell-level clusters can be interpreted more directly than Visium spots, but they
  still require marker, morphology, and segmentation validation.

For an initial Xenium workflow, retaining most biological panel genes after
excluding controls can be more appropriate than automatically restricting the
analysis to highly variable genes.

### Biological questions Xenium supports well

Xenium is particularly useful for:

- Locating rare cells
- Resolving cell types within mixed Visium spots
- Measuring distances between cell populations
- Studying tumor–immune and anatomical boundaries
- Testing whether cell types form spatial neighborhoods
- Examining intracellular transcript localization
- Validating spatial programs discovered with Visium

For example:

```text
Visium discovery:
A region has high immune-response expression.

Xenium follow-up:
Which immune cells are present?
Are they inside or outside the tumor?
Which cells express the response genes?
How far are they from tumor cells?
```

### Recommended Python data model

Xenium contains points, polygons, masks, images, and an expression table. A
`SpatialData` object preserves these layers and their coordinate transformations:

```python
from spatialdata_io import xenium

sdata = xenium("path/to/xenium/output")
sdata.write("xenium.zarr")
```

The resulting object can retain transcript point clouds, cell polygons,
segmentation masks, morphology images, cell-level AnnData tables, and coordinate
systems. This is richer than treating Xenium as only a cell-by-gene matrix.

### Suggested Xenium notebook progression

```text
02_load_and_inspect_xenium_outputs.ipynb
03_xenium_cell_and_transcript_qc.ipynb
04_xenium_segmentation_assessment.ipynb
05_xenium_clustering_and_annotation.ipynb
06_xenium_cell_neighborhoods.ipynb
07_visium_xenium_comparison.ipynb
```

The most important new question is:

> Are these cell-level counts biologically real, or were they created by incorrect
> segmentation and transcript assignment?

## Xenium reading path: simple to complex

1. **Understand the output bundle:** [Understanding Xenium Outputs — 10x Genomics](https://www.10xgenomics.com/support/software/xenium-onboard-analysis/latest/analysis/xoa-output-understanding-outputs). Focus on transcript records, segmentation, Q-scores, controls, and primary versus derived files.

2. **Follow a practical Python example:** [Analyze Xenium Data — Squidpy](https://squidpy.readthedocs.io/en/stable/notebooks/tutorials/tutorial_xenium.html). This introduces SpatialData and shows how images, shapes, points, and expression tables fit together.

3. **Study a biological Visium–Xenium bridge:** [Janesick et al. 2023, *High-resolution mapping of the tumor microenvironment*](https://www.nature.com/articles/s41467-023-43458-x). This combines single-cell RNA sequencing, Visium CytAssist, Xenium, panel design, segmentation, and tumor-boundary analysis.

4. **Learn best-practice QC:** [*Optimizing Xenium In Situ data utility by quality assessment and best-practice analysis workflows*](https://www.nature.com/articles/s41592-025-02617-2). Focus on filtering, segmentation assessment, transcript assignment, normalization, and annotation.

5. **Learn cross-platform inference:** [*Joint imputation and deconvolution of gene expression across spatial transcriptomics platforms*](https://pmc.ncbi.nlm.nih.gov/articles/PMC11870578/). This uses paired Visium–Xenium data and distinguishes direct measurements from predicted out-of-panel expression.

6. **Compare platform trade-offs:** [*A technical comparison of spatial transcriptomics platforms across six cancer types*](https://pmc.ncbi.nlm.nih.gov/articles/PMC12888464/). This compares Visium, Visium HD, Xenium, and CosMx using matched tumor tissues.

7. **Study advanced benchmarking:** [*Systematic benchmarking of high-throughput subcellular spatial transcriptomics platforms across human tumors*](https://www.nature.com/articles/s41467-025-64292-3). This covers sensitivity, specificity, segmentation, cell-type recovery, and cross-platform concordance.

## Data and references

- [Squidpy Visium H&E tutorial](https://squidpy.readthedocs.io/en/latest/notebooks/tutorials/tutorial_visium_hne.html)
- [Squidpy documentation](https://squidpy.readthedocs.io/)
- [Scanpy documentation](https://scanpy.readthedocs.io/)
- [SpatialData documentation](https://spatialdata.scverse.org/)

The example data originates from the 10x Genomics public dataset portal and is
downloaded by `squidpy.datasets.visium_hne_adata()`.

## License

Code in this repository is available under the MIT License. The downloaded example
dataset remains subject to its original provider's terms.
