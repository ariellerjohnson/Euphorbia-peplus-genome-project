#!/bin/bash

#SBATCH --job-name=setting_up_for_genespace.sh
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=arj66@cornell.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=50gb
#SBATCH --time=999:00:00
#SBATCH --output=joblogs/setting_up_for_genespace.joblog

date

mkdir 27_GENESPACE
cd 27_GENESPACE

#paths to peplus files
peplus_longest_isoform=../24_OrthoFinder/Euphorbia_peplus_longest_isoform_only.fa
peplus_gff3=../22_finalizing_genome_for_upload/E_peplus_chromosomal_assembly/Euphorbia_peplus.gff
Orthofinder_protein_dir=../raw_data/OrthoFinder_protein_files/2_Euphorbiaceae
gff3_dir=../raw_data/GENESPACE_gff_files

mkdir results

mkdir -p genomes/EuphorbiaPeplus/v/annotation/
mkdir -p genomes/EuphorbiaLathyris/v/annotation/
mkdir -p genomes/HeveaBrasiliensis/v/annotation/
mkdir -p genomes/RicinusCommunis/v/annotation/
mkdir -p genomes/ManihotEsculenta/v/annotation/

cp $peplus_longest_isoform genomes/EuphorbiaPeplus/v/annotation/protein.fasta
cp $peplus_gff3 genomes/EuphorbiaPeplus/v/annotation/genome.gff3

cp $Orthofinder_protein_dir/E_lathyris.fa genomes/EuphorbiaLathyris/v/annotation/protein.fasta
cp $Orthofinder_protein_dir/H_brasiliensis.fa genomes/HeveaBrasiliensis/v/annotation/protein.fasta
cp $Orthofinder_protein_dir/M_esculenta.fa genomes/ManihotEsculenta/v/annotation/protein.fasta
cp $Orthofinder_protein_dir/R_communis.fa genomes/RicinusCommunis/v/annotation/protein.fasta

cp $gff3_dir/genome.gff3 genomes/EuphorbiaLathyris/v/annotation/genome.gff3
cp $gff3_dir/GCA_010458925.1_ASM1045892v1_genomic.gff3 genomes/HeveaBrasiliensis/v/annotation/genome.gff3
cp $gff3_dir/Mesculenta_671_v8.1.gene.gff3 genomes/ManihotEsculenta/v/annotation/genome.gff3
cp $gff3_dir/Wild_castor_gene.gff3 genomes/RicinusCommunis/v/annotation/genome.gff3

date
