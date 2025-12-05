#!/bin/bash
#SBATCH --job-name=Denoise
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --partition=rra
#SBATCH --qos=rra
#SBATCH --mem=8G
#SBATCH --output=./logs/Denoise_%A_%a.log
#SBATCH --array=0

### ENV ###
eval "$(pixi shell-hook --manifest-path '/work/pi_gblanck/st25/env/gatk')"

set -euo pipefail
### Constants
base_dir="/work/pi_gblanck/st25/my_work/gatk/gbm"
manifest="${base_dir}/sample_pairs.txt"
counts_dir="${base_dir}/counts"
pon_rc="${base_dir}/pon/readcount_pon.hdf5"
outdir="${base_dir}/denoise" 


echo "[ $(date) ] Starting .........."
mkdir -p "$outdir"

### Read manifest line for this task ###
idx=${SLURM_ARRAY_TASK_ID}
line_num=$((idx + 1))

line=$(sed -n "${line_num}p" "$manifest")

caseid=$(echo "$line" | cut -d',' -f1)
tumor=$(echo "$line"  | cut -d',' -f3)

echo "[ $(date) ] Case ID -> $caseid  Processing ${tumor} ........"

gatk DenoiseReadCounts \
  -I "${counts_dir}/${tumor%.bam}.counts.hdf5" \
  --count-panel-of-normals "$pon_rc" \
  --standardized-copy-ratios "${outdir}/${tumor%.bam}.standardizedCR.tsv" \
  --denoised-copy-ratios "${outdir}/${tumor%.bam}.denoisedCR.tsv"

echo "[ $(date) ] Done .........."