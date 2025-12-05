#!/bin/bash
#SBATCH --job-name=PoN
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --partition=rra
#SBATCH --qos=rra
#SBATCH --mem=25G
#SBATCH --output=./logs/PoN%j.log


### ENV ###
eval "$(pixi shell-hook --manifest-path '/work/pi_gblanck/st25/env/gatk')"


set -euo pipefail
### Constants
base_dir="/work/pi_gblanck/st25/my_work/gatk/gbm"
manifest="${base_dir}/sample_pairs.txt"
allelic_dir="${base_dir}/allelic"
counts_dir="${base_dir}/counts"
outdir="${base_dir}/pon" 

echo "[ $(date) ] Starting .........."
# mkdir -p "$allelic_dir"
mkdir -p "$outdir"

echo "[ $(date) ] Creating PoN .........."

gatk CreateReadCountPanelOfNormals \
  $(awk -F ',' -v counts_dir="$counts_dir" '{sub(/\.bam$/,"",$2); printf "-I %s/%s.counts.hdf5 ", counts_dir, $2}' "$manifest") \
  --annotated-intervals "${base_dir}/intervals/intervals.annotated.tsv" \
  -O "${outdir}/readcount_pon.hdf5"


# echo "[ $(date) ] Creating Allelic PoN .........."

# gatk CreateAllelicPanelOfNormals \
#   $(awk -v allelic_dir="$allelic_dir" '{sub(/\.bam$/,"",$0); printf "-I %s/%s.allelic.tsv ", allelic_dir, $0}' "$manifest") \
#   -O "${outdir}/wes_allelic_pon.hdf5"

echo "[ $(date) ] Done .........."