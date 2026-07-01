#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'USAGE'
Usage: scripts/build_report.sh [--export-only] [--skip-export] [--no-render]

Build the HOLOFISHCOLOUR Quarto report from curated analysis outputs.

Options:
  --export-only   Export analysis assets and stop.
  --skip-export   Use already exported assets.
  --no-render     Prepare assets but do not run quarto render.
USAGE
}

export_only=0
skip_export=0
render=1

while [ "$#" -gt 0 ]; do
  case "$1" in
    --export-only) export_only=1 ;;
    --skip-export) skip_export=1 ;;
    --no-render) render=0 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown option: $1" >&2; usage; exit 2 ;;
  esac
  shift
done

QUARTO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)"
ANALYSIS_DIR="${HOLOFISH_ANALYSIS_DIR:-/mnt/lustre/scratch/nlsas/home/otras/pia/dci/holofish_rnaseq/analysis}"
R_ENV="${HOLOFISH_R_ENV:-r45}"

cd "$QUARTO_DIR"
mkdir -p manifests/logs
LOG_FILE="manifests/logs/build_report_$(date +%Y%m%d_%H%M%S).log"

{
  echo "[build] quarto_dir=$QUARTO_DIR"
  echo "[build] analysis_dir=$ANALYSIS_DIR"
  echo "[build] r_env=$R_ENV"
  echo "[build] started=$(date -Is)"

  if [ "$skip_export" -eq 0 ]; then
    echo "[build] exporting curated assets"
    cd "$ANALYSIS_DIR"
    conda run -n "$R_ENV" Rscript scripts/export_quarto_assets.R --quarto-dir="$QUARTO_DIR"
  else
    echo "[build] skipping export"
  fi

  if [ "$export_only" -eq 1 ]; then
    echo "[build] export-only requested"
    echo "[build] finished=$(date -Is)"
    exit 0
  fi

  cd "$QUARTO_DIR"
  echo "[build] preparing DE chapter tables"
  conda run -n "$R_ENV" Rscript scripts/prepare_de_chapter_tables.R

  echo "[build] preparing functional chapter assets"
  FUNCTIONAL_SOURCE_DIR="$ANALYSIS_DIR/results/05_functional_enrichment" \
    conda run -n "$R_ENV" Rscript scripts/prepare_functional_chapter_assets.R

  if [ "$render" -eq 1 ]; then
    echo "[build] rendering Quarto"
    quarto render
  else
    echo "[build] render skipped"
  fi

  echo "[build] finished=$(date -Is)"
} 2>&1 | tee "$LOG_FILE"
