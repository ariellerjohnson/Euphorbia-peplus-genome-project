#!/bin/bash

#SBATCH --job-name=BRAKER.sh
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=arj66@cornell.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH --mem=100gb
#SBATCH --time=999:00:00
#SBATCH --output=joblogs/BRAKER.joblog

date

#files used
softmasked_genome=../13_RepeatMasker/Euphorbia_peplus.fa.masked
bam_file=../11_aligning_RNAseq_to_genome/STAR_Aligned.out.sorted.bam

#BRAKER path
BRAKER_path=/programs/braker2-2.1.6

#Genome name
genome_name=Euphorbia_peplus

#make directory for this step
# mkdir 14_BRAKER
cd 14_BRAKER

#also downloaded GeneMark-ES/ET/EP software and license key
#from http://exon.gatech.edu/GeneMark/license_download.cgi

#move GeneMark-ES/ET/EP software and license key
# mv ../gm_key_64.gz ../gmes_linux_64.tar.gz ./

#decompress key and GeneMark files and get BRAKER
# zcat gm_key_64.gz > $HOME/.gm_key
# tar xvfz gmes_linux_64.tar.gz
# cp -r $BRAKER_path/* ./

#did this on the advice of Cornell BioHPC: https://biohpc.cornell.edu//Lab/userguide.aspx?a=software&i=684#c
#apparently this step fixes ProtHint
# rm -fr gmes_linux_64/ProtHint
# mv ProtHint gmes_linux_64/


#protein only run
# mkdir protein_only_BRAKER
# cd protein_only_BRAKER
#Downloading OrthoDB dataset
# wget https://v100.orthodb.org/download/odb10_plants_fasta.tar.gz --no-check-certificate
# tar xvf odb10_plants_fasta.tar.gz
# cat plants/Rawdata/* > OrthoDB_proteins.fasta
# cp -r ../config/ ./
# cp -r ../gmes_linux_64/ ./
# cp ../$softmasked_genome $genome_name.fa
# export LC_ALL=C
# ../braker2 --genome=$genome_name.fa \
# --prot_seq=OrthoDB_proteins.fasta --softmasking --species=Euphorbia_peplus --cores 40
# cd ../

#RNA only run
mkdir RNA_only_BRAKER
cd RNA_only_BRAKER
cp -r ../config/ ./
cp -r ../gmes_linux_64/ ./
cp ../$softmasked_genome $genome_name.fa
cp ../$bam_file $genome_name.bam
../braker2 --genome=$genome_name.fa \
--bam=$genome_name.bam --softmasking --species=Euphorbia_peplus --cores 40
cd ../

date
