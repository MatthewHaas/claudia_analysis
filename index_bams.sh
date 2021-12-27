#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --time=24:00:00
#SBATCH --mem=30g
#SBATCH --mail-type=ALL
#SBATCH --mail-user=haasx092@umn.edu
#SBATCH -p amdsmall
#SBATCH --account=jkimball
#SBATCH -o index_bams.out
#SBATCH -e index_bams.err

cd /scratch.global/haasx092/claudia_analysis

module load samtools

for i in $(cat claudia_sorted_bam_files.txt);
do
samtools index -c $i
done
