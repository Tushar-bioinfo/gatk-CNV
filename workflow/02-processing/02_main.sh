#!/bin/bash
#SBATCH --job-name=Collect_rc
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --partition=rra
#SBATCH --qos=rra
#SBATCH --mem=8G
#SBATCH --output=./logs/Collect_rc_%A_%a.log
#SBATCH --array=0

set -euo pipefail

### ENV ###
eval "$(pixi shell-hook --manifest-path '/work/pi_gblanck/st25/env/gatk')"

### Paths ###
base_dir="/work/pi_gblanck/st25/my_work/gatk/gbm"
manifest="${base_dir}/sample_pairs.txt"
outdir="${base_dir}/counts"
ref_intervals="${base_dir}/intervals/filtered.interval_list"
bam_dir="${base_dir}/bams"

mkdir -p "$outdir"

echo "[ $(date) ] Starting ......"

### Read manifest line for this task ###
idx=${SLURM_ARRAY_TASK_ID}
line_num=$((idx + 1))

line=$(sed -n "${line_num}p" "$manifest")

caseid=$(echo "$line" | cut -d',' -f1)
normal=$(echo "$line" | cut -d',' -f2)
tumor=$(echo "$line"  | cut -d',' -f3)

echo "[INFO] Case ID     : $caseid"
echo "[INFO] Normal BAM  : $normal"
echo "[INFO] Tumor  BAM  : $tumor"

### Function: Run CollectReadCounts ###
collect_rc() {
    local bam_file="$1"
    local label="$2"

    echo "[ $(date) ] Running CollectReadCounts → ${label}"

    gatk CollectReadCounts \
        -I "${bam_dir}/${bam_file}" \
        -L "${ref_intervals}" \
        --interval-merging-rule OVERLAPPING_ONLY \
        -O "${outdir}/${bam_file%.bam}.counts.hdf5"
}

### Run for tumor ###
collect_rc "$tumor" "tumor"

### Run for normal ###
collect_rc "$normal" "normal"

echo "[ $(date) ] Done ...."
