#!/bin/bash

#SBATCH --job-name=Merqury.sh
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=arj66@cornell.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=50gb
#SBATCH --time=999:00:00
#SBATCH --output=joblogs/Merqury.joblog

date

mkdir 36_Merqury
cd 36_Merqury

#range of possible genome sizes
low_size_estimate=267200000
high_size_estimate=330500000

#get PacBio HiFi reads
cp ../raw_data/HiFi_CCS_reads/PBmixSequel991_2_B01_PCTM_30hours_19kbV2PD_70pM_Euphorbiapeplus.fastq.gz PacBio.fastq.gz

#get genome
cp ../22_finalizing_genome_for_upload/E_peplus_chromosomal_assembly/Euphorbia_peplus.fa chromosome_only_genome.fa

#install meryl
wget https://github.com/marbl/meryl/releases/download/v1.3/meryl-1.3.Linux-amd64.tar.xz
tar -xJf meryl-1.3.*.tar.xz
cd meryl-1.3/bin
export PATH=$PWD:$PATH
cd ../..

#install Merqury
wget https://github.com/marbl/merqury/archive/v1.3.tar.gz
tar -zxvf v1.3.tar.gz
cd merqury-1.3
export MERQURY=$PWD
cd ..
ln -s $MERQURY/merqury.sh

#confirm that best kmer size is 19
$MERQURY/best_k.sh ${low_size_estimate}
$MERQURY/best_k.sh ${high_size_estimate}

#run meryl to get kmer count with k=19
meryl count k=19 PacBio.fastq.gz output pacbio.meryl

#run merqury
./merqury.sh pacbio.meryl chromosome_only_genome.fa merqury_results

date
