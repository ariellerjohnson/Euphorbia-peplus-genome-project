#!/bin/bash

#SBATCH --job-name=trimming_RNAseq.sh
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=arj66@cornell.edu
#SBATCH --ntasks=10
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=4gb
#SBATCH --time=999:00:00
#SBATCH --output=joblogs/trimming_RNAseq.joblog

date

#make directory for this step
mkdir 09_trimming_RNAseq
cd 09_trimming_RNAseq

#files used
RNAseq_raw_dir=../raw_data/RNAseq_data
RNAseq_files=( $(ls $RNAseq_raw_dir))

#path to executable
trimmomatic_path=/programs/trimmomatic/trimmomatic-0.39.jar

#other variables
outdir=$(pwd)

for file in ${RNAseq_files[*]}
do
  eval "long_prefix=$(echo ${file} | cut -d_ -f 1-6)"
  eval "short_prefix=$(echo ${file} | cut -d_ -f 5-6)"
  java -jar $trimmomatic_path PE -threads 2 -trimlog ${outdir}/${long_prefix}.trimlog ${RNAseq_raw_dir}/${long_prefix}_R1.fastq.gz ${RNAseq_raw_dir}/${long_prefix}_R2.fastq.gz -baseout ${outdir}/${short_prefix} ILLUMINACLIP:/programs/trimmomatic/adapters/TruSeq3-PE-2.fa:2:30:10:2:keepBothReads SLIDINGWINDOW:5:20 MINLEN:90 &
done
wait

date
