#!/usr/bin/env bash
# ==============================================================================
# submit_multiome_array.sh
# ==============================================================================
# SLURM submission script for parallelizing process_and_merge_atac.R across samples
# using a SLURM job array.
#
# USAGE:
#   sbatch --array=1-$(sed 1d samples.csv | wc -l) submit_multiome_array.sh samples.csv
# ==============================================================================
#SBATCH --job-name=multiome_atac
#SBATCH --output=logs/multiome_%A_%a.out
#SBATCH --error=logs/multiome_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=06:00:00

set -euo pipefail

SAMPLE_SHEET="${1:-samples.csv}"

if [[ ! -f "$SAMPLE_SHEET" ]]; then
  echo "ERROR: Sample sheet '$SAMPLE_SHEET' not found!"
  exit 1
fi

mkdir -p logs

# Extract the line corresponding to this SLURM array task ID (skipping header)
TASK_ID="${SLURM_ARRAY_TASK_ID:-1}"
ROW=$(sed 1d "$SAMPLE_SHEET" | sed -n "${TASK_ID}p")

if [[ -z "$ROW" ]]; then
  echo "ERROR: No entry found for Task ID $TASK_ID in $SAMPLE_SHEET"
  exit 1
fi

# Parse CSV row columns: sample_id,rna_rds_path,peaks_path,fragments_path,metadata_path,output_dir
SAMPLE_ID=$(echo "$ROW" | awk -F',' '{print $1}')
RNA_RDS=$(echo "$ROW" | awk -F',' '{print $2}')
PEAKS=$(echo "$ROW" | awk -F',' '{print $3}')
FRAGMENTS=$(echo "$ROW" | awk -F',' '{print $4}')
METADATA=$(echo "$ROW" | awk -F',' '{print $5}')
OUT_DIR=$(echo "$ROW" | awk -F',' '{print $6}')

echo "=============================================================================="
echo "Running Multiome ATAC Processing on Cluster Node"
echo "Job ID:        ${SLURM_ARRAY_JOB_ID:-NA}_${SLURM_ARRAY_TASK_ID:-NA}"
echo "Sample ID:     $SAMPLE_ID"
echo "RNA RDS:       $RNA_RDS"
echo "Peaks:         $PEAKS"
echo "Fragments:     $FRAGMENTS"
echo "Metadata:      $METADATA"
echo "Output Dir:    $OUT_DIR"
echo "=============================================================================="

CMD="Rscript process_and_merge_atac.R \
  --sample_id \"$SAMPLE_ID\" \
  --rna_rds \"$RNA_RDS\" \
  --peaks \"$PEAKS\" \
  --fragments \"$FRAGMENTS\" \
  --output_dir \"$OUT_DIR\""

if [[ -n "$METADATA" && "$METADATA" != "NA" && -f "$METADATA" ]]; then
  CMD="$CMD --metadata \"$METADATA\""
fi

eval "$CMD"
echo "Job completed successfully for sample $SAMPLE_ID."
