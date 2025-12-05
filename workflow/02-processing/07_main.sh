#!/bin/bash
#SBATCH --job-name=mseg
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --partition=rra
#SBATCH --qos=rra
#SBATCH --mem=8G
#SBATCH --output=./logs/mseg_%A_%a.log
#SBATCH --array=0

### ENV ###
eval "$(pixi shell-hook --manifest-path '/work/pi_gblanck/st25/env/gatk')"
set -euo pipefail

### Constants ###
base_dir="/work/pi_gblanck/st25/my_work/gatk/gbm"
manifest="${base_dir}/sample_pairs.txt"
denoise_dir="${base_dir}/denoise"
allelic_dir="${base_dir}/allelic"
outdir="${base_dir}/seg"

echo "[ $(date) ] Starting .........."
mkdir -p "$outdir"

### Functions ###

model_segments () {
    local id="$1"
    local tumor="$2"
    local normal="$3"
    echo "[ $(date) ] Running ModelSegments for $id"

    gatk ModelSegments \
      --denoised-copy-ratios "${denoise_dir}/${tumor%.bam}.denoisedCR.tsv" \
      --allelic-counts "${allelic_dir}/${tumor%.bam}.allelic.tsv" \
      --normal-allelic-counts "${allelic_dir}/${normal%.bam}.allelic.tsv" \
      --output-prefix tumor \
      -O "${outdir}/${id}"
}


call_segments (){
    local id="$1"
    echo "[ $(date) ] Running ModelSegments for $id"
    gatk CallCopyRatioSegments \
    -I "${outdir}/${id}/tumor.cr.seg" \
    -O "${outdir}/${id}/tumor.called.seg"
}

### Main ###

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

# Run
model_segments "$id" "$tumor" "$normal"

echo "Calling ..."

call_segments "$id"


echo "[ $(date) ] Done .........."
