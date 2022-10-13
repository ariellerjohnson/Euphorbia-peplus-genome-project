#!/bin/bash

#SBATCH --job-name=EDTA.sh
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=arj66@cornell.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=50gb
#SBATCH --time=999:00:00
#SBATCH --output=joblogs/EDTA_non_chrom_scaffolds.joblog

#make directory for this step
mkdir 38_EDTA_non_chrom_scaffolds
cd 38_EDTA_non_chrom_scaffolds

#files used
genome_file=../06_3D_DNA_post_adjustment/E_peplus_hifiasm.bp.p_ctg.FINAL.fasta

#seqtk path
seqtk_path=/programs/seqtk/seqtk

#genome name
genome_name=E_peplus

#copy genome and select scaffolds only
cp $genome_file ${genome_name}.fa

scaf_numbers=$(seq 9 1242)

for number in ${scaf_numbers[*]}
do
  HiC="HiC_scaffold_"${number}
  echo $HiC >> my_scaffolds.txt
done

$seqtk_path subseq ${genome_name}.fa my_scaffolds.txt > non_chrom_scaffolds_only.fa

#https://github.com/oushujun/EDTA#installation
SINGULARITY_CACHEDIR=./
export SINGULARITY_CACHEDIR
singularity pull EDTA.sif docker://oushujun/edta:2.0.0

#run EDTA
export LC_ALL=C
singularity exec --no-home ./EDTA.sif EDTA.pl --genome non_chrom_scaffolds_only.fa --threads 20 --anno 1
