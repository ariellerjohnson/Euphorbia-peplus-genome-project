#!/bin/bash

#SBATCH --job-name=aligning_HiC_to_chromosomes.sh
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=arj66@cornell.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=75gb
#SBATCH --time=999:00:00
#SBATCH --output=joblogs/aligning_HiC_to_chromosomes.joblog

date

#make directory for this step
mkdir 37_aligning_HiC_to_chromosomes
cd 37_aligning_HiC_to_chromosomes

#files used
HiC_forward=../raw_data/HiC_data/12221_11073_130590_HKT53AFX2_E_peplus_Hi_C_ATGTGGAC_TTGGTGAG_R1.fastq.gz
HiC_reverse=../raw_data/HiC_data/12221_11073_130590_HKT53AFX2_E_peplus_Hi_C_ATGTGGAC_TTGGTGAG_R2.fastq.gz
genome=../22_finalizing_genome_for_upload/E_peplus_chromosomal_assembly/Euphorbia_peplus.fa
genome_all_scaffolds=../22_finalizing_genome_for_upload/E_peplus_assembly_including_all_scaffolds/E_peplus_all_scaffolds.fa

#path to executables
export PATH=/programs/STAR-2.7.5a/bin/Linux_x86_64_static:$PATH
export PATH=/programs/samtools-1.14/bin:$PATH
seqtk_path=/programs/seqtk/seqtk

#copy HiC data
cp $HiC_forward HiC_forward.fastq.gz
cp $HiC_reverse HiC_reverse.fastq.gz
gunzip HiC_*.fastq.gz

#get non-chromosomes only
scaf_numbers=$(seq 9 1242)

for number in ${scaf_numbers[*]}
do
  HiC="HiC_scaffold_"${number}
  echo $HiC >> my_scaffolds.txt
done

$seqtk_path subseq $genome_all_scaffolds my_scaffolds.txt > non_chrom_scaffolds_only.fa

#make genomeDir directories for STAR
mkdir STARgenome_chromosomes
mkdir STARgenome_other_scaffolds

#make STAR genome index
STAR --runMode genomeGenerate --genomeSAindexNbases 12 --runThreadN 20 \
--genomeDir STARgenome_chromosomes --genomeFastaFiles $genome

STAR --runMode genomeGenerate --genomeSAindexNbases 12 --runThreadN 20 \
--genomeDir STARgenome_other_scaffolds --genomeFastaFiles non_chrom_scaffolds_only.fa

#run STAR on each pair of RNAseq files separately
STAR --twopassMode Basic --genomeDir STARgenome_chromosomes \
--runThreadN 20 --readFilesIn HiC_forward.fastq,HiC_reverse.fastq \
--outFileNamePrefix chromosomes --outSAMtype BAM Unsorted \
--limitOutSJcollapsed 5000000 --limitSjdbInsertNsj 2000000

STAR --twopassMode Basic --genomeDir STARgenome_other_scaffolds \
--runThreadN 20 --readFilesIn HiC_forward.fastq,HiC_reverse.fastq \
--outFileNamePrefix other_scaffolds --outSAMtype BAM Unsorted \
--limitOutSJcollapsed 5000000 --limitSjdbInsertNsj 2000000

date
