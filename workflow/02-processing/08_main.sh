#!/bin/bash
#SBATCH --job-name=plots
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=12G
#SBATCH --partition=rra
#SBATCH --qos=rra
#SBATCH --output=./logs/plots_%A_%a.log
#SBATCH --array=0


module purge
module load apps/apptainer/1.3.5        
set -euo pipefail


# Constants
base_dir="/work/pi_gblanck/st25/my_work/gatk/gbm"
manifest="${base_dir}/sample_pairs.txt"
dict="/work/pi_gblanck/st25/my_work/genomes/gdc/GRCh38.d1.vd1.dict"
denoise_dir="${base_dir}/denoise"
seg_dir="${base_dir}/seg"
outdir="${base_dir}/plots"



echo "[$(date)] Starting ......."
mkdir -p "$outdir"

# use the official image; bind PWD so container sees files
#IMG="docker://broadinstitute/gatk:4.6.2.0"
IMG="/work/pi_gblanck/st25/env/containers/gatk_4.6.2.0.sif"
WD="$PWD"

plotting () {
    local id="$1"
    local tumor="$2"

    echo "[$(date)] Running PlotDenoisedCopyRatios and PlotModeledSegments in one container ..."

    apptainer exec -B "$WD:$WD" --pwd "$WD" "$IMG" bash -c "
        set -euo pipefail

        echo '[INFO] PlotDenoisedCopyRatios ...'
        gatk PlotDenoisedCopyRatios \
          --standardized-copy-ratios '${denoise_dir}/${tumor%.bam}.standardizedCR.tsv' \
          --denoised-copy-ratios     '${denoise_dir}/${tumor%.bam}.denoisedCR.tsv' \
          --sequence-dictionary      '$dict' \
          --output-prefix            tumor \
          -O                         '${outdir}/${id}'

        echo '[INFO] PlotModeledSegments ...'
        gatk PlotModeledSegments \
          --denoised-copy-ratios '${denoise_dir}/${tumor%.bam}.denoisedCR.tsv' \
          --allelic-counts       '${seg_dir}/${id}/tumor.hets.tsv' \
          --segments             '${seg_dir}/${id}/tumor.modelFinal.seg' \
          --sequence-dictionary  '$dict' \
          --output-prefix        tumor \
          -O                     '${outdir}/${id}'
    "
}


### Main ###
### Read manifest line for this task ###
idx=${SLURM_ARRAY_TASK_ID}
line_num=$((idx + 1))

line=$(sed -n "${line_num}p" "$manifest")

caseid=$(echo "$line" | cut -d',' -f1)
# normal=$(echo "$line" | cut -d',' -f2)
tumor=$(echo "$line"  | cut -d',' -f3)

echo "[INFO] Case ID     : $caseid"
echo "[INFO] Tumor  BAM  : $tumor"

plotting "$caseid" "$tumor"

echo "[$(date)] Done .........."