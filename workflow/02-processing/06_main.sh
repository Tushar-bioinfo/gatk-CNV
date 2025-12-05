#!/bin/bash
#SBATCH --job-name=Ale-Counts
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --partition=rra
#SBATCH --qos=rra
#SBATCH --mem=16G
#SBATCH --output=./logs/Ale-Counts_%A_%a.log
#SBATCH --array=0

### ENV ###
eval "$(pixi shell-hook --manifest-path '/work/pi_gblanck/st25/env/gatk')"

set -euo pipefail

### Constants ###
base_dir="/work/pi_gblanck/st25/my_work/gatk/gbm"
ref_fasta="/work/pi_gblanck/st25/my_work/genomes/gdc/GRCh38.d1.vd1.fa"
manifest="${base_dir}/sample_pairs.txt"
ref_dir="${base_dir}/refs"
outdir="${base_dir}/allelic"
bam_dir="${base_dir}/bams"

echo "[ $(date) ] Starting .........."
mkdir -p "$outdir" 

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



compute_ac() {
    local bam_file="$1"
    local label="$2"

    echo "[ $(date) ] Running CollectAllelicCounts for ${label}: ${bam_file}" >&2

    gatk CollectAllelicCounts \
        -I "${bam_dir}/${bam_file}" \
        -L "${ref_dir}/broad_hc_snp_intersect.vcf.gz" \
        -R "$ref_fasta" \
        -O "${outdir}/${bam_file%.bam}.allelic.tsv"
}

### Main ###
echo "[ $(date) Case ID:  $caseid"
### Run for tumor ###
echo "[ $(date) ] Processing Tumor:  ${tumor}"
compute_ac "$tumor" "tumor"

### Run for normal ###
echo "[ $(date) ] Processing Normal:  ${normal}"
compute_ac "$normal" "normal"

echo "[ $(date) ] Done .........."
