# HOLOFISHCOLOUR Transcriptomics Report

This repository contains the Quarto/GitHub Pages report for the turbot
(*Scophthalmus maximus*) malpigmentation RNA-seq project.

The statistical analysis is not maintained here. The authoritative analysis
project is:

```text
/mnt/lustre/scratch/nlsas/home/otras/pia/dci/holofish_rnaseq
```

## Main Rule

Use this repository as a publication layer. It consumes curated tables, figures
and MultiQC reports exported from the analysis project.

## Build

From this repository:

```bash
scripts/build_report.sh
```

Useful variants:

```bash
scripts/build_report.sh --export-only
scripts/build_report.sh --no-render
scripts/build_report.sh --skip-export
```

The build script performs:

1. export curated assets from `holofish_rnaseq/analysis`;
2. prepare expression-differential chapter tables;
3. prepare functional chapter assets;
4. render the Quarto book into `docs/`.

## Directory Map

```text
chapters/       Quarto chapter sources
appendices/     Quarto appendix sources
config/         report-side metadata exported from analysis config
figures/        curated report figures
tables/         curated report tables
multiqc/        selected MultiQC reports and data
manifests/      export manifests and build logs
docs/           rendered site published by GitHub Pages
scripts/        report preparation/build scripts
```

## Traceability

- `manifests/quarto_export_manifest.tsv` records what was exported from the
  analysis project.
- `config/contrast_metadata.tsv` is generated from
  `analysis/config/analysis_config.yml` and is the report-side source for
  official contrast labels.
- `references.bib` remains the report bibliography. The Zotero transcriptomics
  export lives in `holofish_rnaseq/analysis/references/`; only missing report
  entries should be appended here to avoid replacing the broader bibliography.
- `appendices/reproducibility.qmd` documents render and analytical decisions.

## GitHub Pages

The workflow in `.github/workflows/pages.yml` deploys the existing `docs/`
folder. It does not re-render Quarto on GitHub.
