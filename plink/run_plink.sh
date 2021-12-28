#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --time=24:00:00
#SBATCH --mem=30g
#SBATCH --mail-type=ALL
#SBATCH --mail-user=haasx092@umn.edu
#SBATCH -p amdsmall
#SBATCH --account=jkimball
#SBATCH -o run_plink.out
#SBATCH -e run_plink.err

cd /scratch.global/haasx092/claudia_analysis/211227_snp_calling_results

module load plink

plink --vcf merged_vcf_files.vcf --mind 0.99 --double-id --allow-extra-chr --recode --out claudia_analysis

# PCA calculation
plink --pca --file claudia_analysis --allow-extra-chr -out claudia_analysis_pca
