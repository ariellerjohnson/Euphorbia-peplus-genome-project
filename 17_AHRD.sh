#!/bin/bash

#SBATCH --job-name=AHRD.sh
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=arj66@cornell.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH --mem=75gb
#SBATCH --time=999:00:00
#SBATCH --output=joblogs/AHRD.joblog

date

mkdir 17_AHRD
cd 17_AHRD

#get AHRD
git clone https://github.com/groupschoof/AHRD.git
cd AHRD
git checkout tags/v3.3.3
/programs/apache-ant-1.9.9/bin/ant dist
cd ..

#files used
protein_fasta=../15_TSEBRA/Euphorbia_peplus.aa

#Genome name
genome_name=Euphorbia_peplus

#copy proteins
cp $protein_fasta ${genome_name}.aa

#SwissProt
wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
gunzip uniprot_sprot.fasta.gz
# creating a diamond-formatted database file
/programs/diamond/diamond makedb --in uniprot_sprot.fasta -d SwissProt
# running a search in blastp mode
/programs/diamond/diamond blastp -d SwissProt -q ${genome_name}.aa -o query_vs_SwissProt.tsv -p 40 --evalue 1e-5 --outfmt 6 > swissprot_blast.log

#TAIR
wget https://www.arabidopsis.org/download_files/Proteins/TAIR10_protein_lists/TAIR10_pep_20110103_representative_gene_model --no-check-certificate
# creating a diamond-formatted database file
/programs/diamond/diamond makedb --in TAIR10_pep_20110103_representative_gene_model -d TAIR
# running a search in blastp mode
/programs/diamond/diamond blastp -d TAIR -q ${genome_name}.aa -o query_vs_TAIR.tsv -p 40 --evalue 1e-5 --outfmt 6 > tair_blast.log

#run AHRD
java -Xmx2g -jar ./AHRD/dist/ahrd.jar ../17_1_my_ahrd_config.yml > ahrd_log.out

date
