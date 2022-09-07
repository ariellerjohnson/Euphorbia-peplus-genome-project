#!/bin/bash

#SBATCH --job-name=BLAST_for_BLAST2GO.sh
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=arj66@cornell.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH --mem=50gb
#SBATCH --time=999:00:00
#SBATCH --output=joblogs/BLAST_for_BLAST2GO.joblog

date

#make directory for this step
mkdir 19_BLAST_for_BLAST2GO
cd 19_BLAST_for_BLAST2GO

#get UniRef90 database (Cornell BioHPC has it downloaded)
cp /shared_data/genome_db/uniref90.dmnd ./

#Path to diamond
diamond_path=/programs/diamond/diamond

#files used
protein_fasta=../15_TSEBRA/Euphorbia_peplus.aa

#Genome name
genome_name=Euphorbia_peplus

#copy proteins
cp $protein_fasta ${genome_name}.aa

#run diamond blastp on uniref90 database
$diamond_path blastp --db uniref90 --query ${genome_name}.aa \
--outfmt 5 --max-target-seqs 100 --max-hsps 1 --evalue 1e-10  -t ./ --block-size 10 \
--index-chunks 1 -o uniref_blastresults.xml -p 40 > uniref_diamond_log.out

#rm uniref90.dmnd

date
