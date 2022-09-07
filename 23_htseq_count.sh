#!/bin/bash

#SBATCH --job-name=htseq_count.sh
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=arj66@cornell.edu
#SBATCH --ntasks=30
#SBATCH --cpus-per-task=1
#SBATCH --mem=50gb
#SBATCH --time=999:00:00
#SBATCH --output=joblogs/htseq_count.joblog

date

#make directory for this step
mkdir 23_htseq_count
cd 23_htseq_count

#path to htseq_count
export PYTHONPATH=/programs/HTSeq-0.11.2/lib64/python3.6/site-packages/
export PATH=/programs/HTSeq-0.11.2/bin:$PATH

#path to gtftools.py
gtftools_path=/programs/GTFtools_0.6.5/gtftools.py

#files used
BAMDIR=../11_aligning_RNAseq_to_genome
BAMs=$BAMDIR/*.sorted.bam
gtf_file=../15_TSEBRA/Euphorbia_peplus.gtf

#Genome name
genome_name=Euphorbia_peplus

#copy gtf
cp $gtf_file ${genome_name}.gtf

#removing the gene and transcript lines
python3 ../23_1_removing_gene_and_transcript_gtf_lines.py ${genome_name}.gtf

#get rid of zeroes in chromosome number since bam files don't have them
awk -v OFS='\t' '{gsub ("_0*", "_", $1)};1' < ${genome_name}_no_gene_or_transcript.gtf > ${genome_name}_no_zero.gtf

#do htseq-count for every file
for bam_file in ${BAMs[*]}
do
  outfile_name=$(echo ${bam_file} | awk '{ print substr( $0, 1, length($0)-4 ) }' | awk '{print $0".counts"}')
  htseq-count -s no -r pos -f bam $bam_file ${genome_name}_no_zero.gtf > $outfile_name &
done
wait

mv $BAMDIR/*.counts ./

#convert gff to ensembl-ish format
awk '{OFS="\t"} $1=1' ${genome_name}.gtf > ensembl.gtf

#use gtftools to get genelength
$gtftools_path -l ${genome_name}.genelength ensembl.gtf

date
