#!/bin/bash

#SBATCH --job-name=OrthoFinder.sh
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=arj66@cornell.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH --mem=50gb
#SBATCH --time=999:00:00
#SBATCH --output=joblogs/OrthoFinder.joblog

date

mkdir 24_OrthoFinder
cd 24_OrthoFinder

#files used
protein_fasta=../22_finalizing_genome_for_upload/E_peplus_chromosomal_assembly/Euphorbia_peplus.aa
other_plants_protein_fasta_dir=../raw_data/OrthoFinder_protein_files

#Genome name
genome_name=Euphorbia_peplus

#copy proteins
cp $protein_fasta ${genome_name}.fasta

#Getting OrthoFinder
wget https://github.com/davidemms/OrthoFinder/releases/latest/download/OrthoFinder.tar.gz
tar xvfz OrthoFinder.tar.gz

#Get longest isoform only
#https://bioinformatics.stackexchange.com/questions/595/how-can-longest-isoforms-per-gene-be-extracted-from-a-fasta-file/603#603
python3 ../24_1_Ryan_filter.py ${genome_name}.fasta > ${genome_name}_longest_isoform_only.fa

#see how many protein sequences this removed from the file
grep ">" ${genome_name}.aa | wc -l
grep ">" ${genome_name}_longest_isoform_only.fa | wc -l

#checking the 1st 50 lines visually
grep "^>" ${genome_name}.aa | head -n50
grep "^>" ${genome_name}_longest_isoform_only.fa | head -n50

#list of transcript names only
grep "^>" ${genome_name}_longest_isoform_only.fa | sed 's/>//g' > ${genome_name}_longest_isoform_names.out

mkdir 1_all_species
cp $other_plants_protein_fasta_dir/1_all_species/*.fa 1_all_species
cp ${genome_name}_longest_isoform_only.fa 1_all_species/E_peplus.fa
mkdir 2_Euphorbiaceae
cp $other_plants_protein_fasta_dir/2_Euphorbiaceae/*.fa 2_Euphorbiaceae
cp ${genome_name}_longest_isoform_only.fa 2_Euphorbiaceae/E_peplus.fa
mkdir 3_Euphorbia_only
cp $other_plants_protein_fasta_dir/3_Euphorbia_only/*.fa 3_Euphorbia_only
cp ${genome_name}_longest_isoform_only.fa 3_Euphorbia_only/E_peplus.fa

#Activate OrthoFinder conda environment on BioHPC at Cornell
source /programs/miniconda3/bin/activate orthofinder

orthofinder -S diamond -t 40 -a 8 -f 1_all_species

orthofinder -S diamond -t 40 -a 8 -f 2_Euphorbiaceae

orthofinder -S diamond -t 40 -a 8 -f 3_Euphorbia_only

date
