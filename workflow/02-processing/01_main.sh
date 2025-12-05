#!/bin/bash
#SBATCH --job-name=Intervals
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=5G
#SBATCH --partition=rra
#SBATCH --qos=rra
#SBATCH --output=./logs/Intervals_%A_%a.log

### ENV activate 
eval "$(pixi shell-hook --manifest-path '/work/pi_gblanck/st25/env/gatk')"

set -euo pipefail

### Constants
base_dir="/work/pi_gblanck/st25/my_work/gatk/gbm"
ref_dir="/work/pi_gblanck/st25/my_work/genomes/gdc/"
interval_dir="${base_dir}/intervals"
bed_file="${base_dir}/autosome_exons_filtered.bed"

echo "[ $(date) ] Starting .........."

mkdir -p "$interval_dir"


### Main
echo "[ $(date) ] Preprocessing Intervals .........."
gatk PreprocessIntervals \
  -R "${ref_dir}/GRCh38.d1.vd1.fa" \
  -L "${bed_file}" \
  --bin-length 0 \
  --padding 250 \
  --interval-merging-rule OVERLAPPING_ONLY \
  -O "${interval_dir}/preprocessed.interval_list"

if [ -s "${interval_dir}/preprocessed.interval_list" ]; then
    echo "[Log] File Exists Continue ...."
else 
    echo "[Error] preprocessed.interval_list does not Exist"
    exit 1
fi


echo "[ $(date) ] Annotating Intervals .........."
gatk AnnotateIntervals \
  -R "${ref_dir}/GRCh38.d1.vd1.fa" \
  -L "${interval_dir}/preprocessed.interval_list" \
  --interval-merging-rule OVERLAPPING_ONLY \
  -O "${interval_dir}/intervals.annotated.tsv"

if [ -s "${interval_dir}/intervals.annotated.tsv" ]; then
    echo "[Log] File Exists Continue ...."
else 
    echo "[Error] intervals.annotated.tsv does not Exist"
    exit 1
fi

echo "[ $(date) ] Filtering Intervals .........."
  gatk FilterIntervals \
  -L "${interval_dir}/preprocessed.interval_list" \
  --annotated-intervals "${interval_dir}/intervals.annotated.tsv" \
  --interval-merging-rule OVERLAPPING_ONLY \
  -O "${interval_dir}/filtered.interval_list"

if [ -s "${interval_dir}/filtered.interval_list" ]; then
    echo "[Log] File Exists Continue ...."
else 
    echo "[Error] filtered.interval_list does not Exist"
    exit 1
fi

echo "[ $(date) ] Done .........."

