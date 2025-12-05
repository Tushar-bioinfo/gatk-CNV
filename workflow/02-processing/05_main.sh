#!/bin/bash
#SBATCH --job-name=Select_variants
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --partition=rra
#SBATCH --qos=rra
#SBATCH --mem=50G
#SBATCH --output=./logs/Sel-vars_%j.log

### ENV ###
eval "$(pixi shell-hook --manifest-path '/work/pi_gblanck/st25/env/gatk')"

set -euo pipefail
### Constants
base_dir="/work/pi_gblanck/st25/my_work/gatk/gbm"
ref_fasta="/work/pi_gblanck/st25/my_work/genomes/gdc/GRCh38.d1.vd1.fa"
ref_dir="${base_dir}/refs"


echo "[ $(date) ] Starting .........."


gatk SelectVariants \
  -R "$ref_fasta" \
  -V "${ref_dir}/broad_hc_snp.vcf.gz" \
  --select-type-to-include SNP \
  --restrict-alleles-to BIALLELIC \
  -O "${ref_dir}/broad_hc_snp_biallelic.vcf.gz"


echo "[ $(date) ] Intersecting .........."

bedtools intersect -header \
  -a "${ref_dir}/broad_hc_snp_biallelic.vcf.gz" \
  -b "${base_dir}/autosome_exons_filtered.bed" \
  > "${ref_dir}/broad_hc_snp_intersect.vcf"

bgzip -f "${ref_dir}/broad_hc_snp_intersect.vcf"
tabix -f -p vcf "${ref_dir}/broad_hc_snp_intersect.vcf.gz"


echo "[ $(date) ] Done .........."