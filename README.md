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
