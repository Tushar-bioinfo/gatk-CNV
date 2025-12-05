#!/bin/bash
#SBATCH --job-name=gtf2bed
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=5G
#SBATCH --partition=rra
#SBATCH --qos=rra
#SBATCH --output=./logs/gtf2bed_%A_%a.log

eval "$(pixi shell-hook --manifest-path '/work/pi_gblanck/st25/env/gatk')"

gtf="./gencode.v36.annotation.gtf"
exons="exons.bed"

echo "[ $(date) ] Starting .........."

awk 'BEGIN{OFS="\t"}
  $3=="exon" {
    n=0
    if ($1 ~ /^chr[0-9]+$/)      n = substr($1,4)+0     # chr10 → 10
    else if ($1 ~ /^[0-9]+$/)    n = $1+0               # 10    → 10
    # keep autosomes 1–22 only
    if (n>=1 && n<=22) print "chr" n, $4-1, $5
  }' "$gtf" \
| m -i - \
| bedtools merge -i - \
> "$exons"

echo "[ $(date) ] Done .........."