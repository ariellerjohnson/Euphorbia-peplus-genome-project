#!/bin/bash

#SBATCH --job-name=aligning_RNAseq_to_genome.sh
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=arj66@cornell.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH --mem=100gb
#SBATCH --time=999:00:00
#SBATCH --output=joblogs/aligning_RNAseq_to_genome.joblog

date

#make directory for this step
mkdir 11_aligning_RNAseq_to_genome
cd 11_aligning_RNAseq_to_genome

#files used
RNAseq_trimmed_dir=../09_trimming_RNAseq
post_review_assembly=../06_3D_DNA_post_adjustment/E_peplus_hifiasm.bp.p_ctg.FINAL.fasta

#path to executables
export PATH=/programs/STAR-2.7.5a/bin/Linux_x86_64_static:$PATH
export PATH=/programs/samtools-1.14/bin:$PATH

#make genomeDir directory for STAR
mkdir STARgenome

#make STAR genome index
STAR --runMode genomeGenerate --genomeSAindexNbases 12 --runThreadN 40 \
--genomeDir STARgenome --genomeFastaFiles $post_review_assembly

#get the forward read file names
f_list=$RNAseq_trimmed_dir/*1P

#run STAR on each pair of RNAseq files separately
for file in ${f_list[*]}
do
  r_name=$(echo ${file} | awk '{ print substr( $0, 1, length($0)-2 ) }' | awk '{print $0"2P"}')
  out_name=$(echo ${file##*/} | awk '{ print substr( $0, 1, length($0)-2 ) }')
  STAR --twopassMode Basic --genomeDir STARgenome \
  --runThreadN 40 --readFilesIn $file $r_name \
  --outFileNamePrefix $out_name --outSAMtype BAM Unsorted
done

#the aligned files
aligned_list=*_Aligned.out.bam

#sort bam files
for aligned_file in ${aligned_list[*]}
do
  sorted_name=$(echo ${aligned_file} | awk '{ print substr( $0, 1, length($0)-3 ) }' | awk '{print $0"sorted.bam"}')
  samtools sort -@ 40 -O BAM -o $sorted_name $aligned_file
  samtools index -@ 40 $sorted_name
done

#combine bam files
samtools merge -@ 40 STAR_Aligned.out.sorted.bam *sorted.bam

#index combined bam file
samtools index -@ 40 STAR_Aligned.out.sorted.bam

date
