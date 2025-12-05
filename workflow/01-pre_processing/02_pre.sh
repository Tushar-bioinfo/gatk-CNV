#!/bin/bash
#SBATCH --job-name=mos
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=12G
#SBATCH --ntasks=1
#SBATCH --partition=rra
#SBATCH --qos=rra
#SBATCH --output=./logs/mosdepth_%A_%a.out
#SBATCH --array=0

set -euo pipefail

# ---- user edits ----
MANIFEST="./sample_pairs.txt"
BAM_DIR="./bams"
BED="autosome_exons_filtered.bed"
outdir="./mosdepth"

source /home/s/st25/miniconda3/etc/profile.d/conda.sh
conda activate sra

mkdir -p "$outdir"

echo "[ $(date) ] Starting .........."

idx="${SLURM_ARRAY_TASK_ID}"
line_num=$((idx + 1))

caseid=$( sed -n "${line_num}p" "$MANIFEST" | cut -f1 -d',' )
filename=$( sed -n "${line_num}p" "$MANIFEST" | cut -f2 -d',' )

BAM_PATH="${BAM_DIR}/${filename}"

echo "[ $(date) ] Running mosdepth for sample: $caseid"
echo "[ $(date) ] BAM: $BAM_PATH"
echo "[ $(date) ] BED: $BED"

# --- Ensure BAM index exists and has correct name ---
if [[ ! -f "${BAM_PATH}.bai" ]]; then
    echo "[ $(date) ] No .bam.bai index — checking for matching .bai..."
    alt_bai="${BAM_PATH%.bam}.bai"
    if [[ -f "$alt_bai" ]]; then
        echo "[ $(date) ] Found $alt_bai — creating symlink ${BAM_PATH}.bai"
        ln -sf "$alt_bai" "${BAM_PATH}.bai"
    else
        echo "[ $(date) ] ERROR: No BAM index found for $filename"
        exit 1
    fi
fi

pushd "$outdir" >/dev/null
mosdepth \
  --by "$BED" \
  --no-per-base \
  --threads 8 \
  "$caseid" \
  "$BAM_PATH"
popd >/dev/null

echo "[ $(date) ] Done .........."
