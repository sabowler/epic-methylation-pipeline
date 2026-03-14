# epic-methylation-pipeline

**Unsupervised analysis pipeline for Illumina EPIC v2 methylation array data.**

Built and maintained by [Scott A. Bowler](https://github.com/sabowler) — Ndhlovu Lab, Weill Cornell Medicine.

Processes raw IDAT files from 850K+ CpG EPIC v2 arrays through QC, normalization, and a full suite of unsupervised analyses — PCA, hierarchical clustering, heatmaps, MDS, and UMAP — with all results exported as publication-ready PDFs and CSV files.

---

## Quickstart

```bash
# 1. Place your sample directories in a working directory:
#    /your/data/
#      SAMPLE001/
#        <basename>_Grn.idat
#        <basename>_Red.idat
#      SAMPLE002002/
#        ...

# 2. Edit the setwd() line in methyl.R to point to your data directory

# 3. Run:
Rscript methyl.R
```

---

## Requirements

- R >= 4.2
- Bioconductor packages: `minfi`, `IlluminaHumanMethylationEPICv2manifest`, `IlluminaHumanMethylationEPICv2anno.20a1.hg38`
- CRAN packages: `ggplot2`, `pheatmap`, `RColorBrewer`, `gridExtra`, `ggrepel`, `umap`

Install on first use by uncommenting the installation block at the top of `methyl.R`.

---

## Input

| Requirement | Detail |
|---|---|
| Array type | Illumina EPIC v2 (935K probes) |
| File format | Raw IDAT pairs (`_Grn.idat` + `_Red.idat`) |
| Directory convention | One folder per sample named `SAMPLE001`, `SAMPLE002`, etc. |

---

## Pipeline Steps

| Step | Method | Output |
|---|---|---|
| 1. IDAT loading | minfi `read.metharray()` | RGChannelSet |
| 2. Quality control | Detection p-values | `qc_report.pdf`, `qc_metrics.csv` |
| 3. Normalization | NOOB background correction | Cleaned beta + M values, `processed_data.RData` |
| 4. PCA | Top 10,000 variable probes | `pca_plots.pdf`, `pca_coordinates.csv` |
| 5. Hierarchical clustering | Complete linkage, k=2–5 | `hierarchical_clustering.pdf`, `cluster_assignments.csv` |
| 6. Heatmap | Top 1,000 variable probes | `heatmap_top_variable_probes.pdf` |
| 7. Correlation matrix | Pearson correlation | `sample_correlation_heatmap.pdf` |
| 8. MDS | Classical multidimensional scaling | `mds_plots.pdf` |
| 9. UMAP | UMAP projection (optional) | `umap_plot.pdf`, `umap_coordinates.csv` |
| 10. Summary | Run metadata | `analysis_summary.txt` |

All outputs are written to `methylation_analysis_results/` in the working directory.

---

## QC Thresholds

| Metric | Threshold | Action |
|---|---|---|
| Mean detection p-value | > 0.01 | Flag sample |
| Failed probes per sample | > 5% | Flag sample |
| Failed probes across samples | > 50% of samples | Remove probe |

SNP-associated and cross-reactive probes are removed using `rmSNPandCH()` where the EPIC v2 annotation package is available.

---

## Publications

Methylation analysis methods developed in support of:

- Bowler SA et al. *A machine learning approach utilizing DNA methylation as an accurate classifier of COVID-19 disease severity.* Scientific Reports (2022)

Full publication list: [Google Scholar](https://scholar.google.com/citations?user=qM-DZhMAAAAJ)

---

## License

MIT
