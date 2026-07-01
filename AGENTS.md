# Codex Instructions For Transcriptomics Report

This repository is the Quarto publication layer for the HOLOFISHCOLOUR turbot
RNA-seq analysis.

## Start Here

Read first:

```text
README.md
appendices/reproducibility.qmd
manifests/quarto_export_manifest.tsv
```

The authoritative analysis project is:

```text
/mnt/lustre/scratch/nlsas/home/otras/pia/dci/holofish_rnaseq
```

Before changing analysis claims, check the analysis-side handoff files:

```text
/mnt/lustre/scratch/nlsas/home/otras/pia/dci/holofish_rnaseq/analysis/context_handoff/
```

## Rules

- Do not recalculate DESeq2 or enrichment in this repository.
- Export analysis assets with `scripts/build_report.sh --export-only` or with
  `holofish_rnaseq/analysis/scripts/export_quarto_assets.R`.
- Render only when explicitly requested.
- Keep `ENSSMAG` as the primary gene identifier.
- Keep `docs/` as rendered output for GitHub Pages.
- Keep generated manifests and build logs in `manifests/`.

## Main Commands

```bash
scripts/build_report.sh --no-render
scripts/build_report.sh
```

Use `conda run -n r45 Rscript ...` for R scripts.
