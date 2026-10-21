#!/usr/bin/env bash
# ==============================================================================
# submit_multiome_batch.sh
# ==============================================================================
# SLURM submission script for running process_and_merge_atac.R in batch mode
# across all samples defined in samples.csv.
#
# USAGE:
#   sbatch submit_multiome_batch.sh samples.csv
# ==============================================================================
#SBATCH --job-name=multiome_batch
#SBATCH --output=logs/multiome_batch_%j.out
#SBATCH --error=logs/multiome_batch_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=12:00:00

set -euo pipefail

SAMPLE_SHEET="${1:-samples.csv}"

if [[ ! -f "$SAMPLE_SHEET" ]]; then
  echo "ERROR: Sample sheet '$SAMPLE_SHEET' not found!"
  exit 1
fi

mkdir -p logs

echo "=============================================================================="
echo "Starting Multiome Batch Processing on Cluster Node"
echo "Job ID:        ${SLURM_JOB_ID:-interactive}"
echo "Sample Sheet:  $SAMPLE_SHEET"
echo "Date:          $(date)"
echo "=============================================================================="

Rscript process_and_merge_atac.R --sample_sheet "$SAMPLE_SHEET"

echo "Batch processing finished successfully at $(date)."
