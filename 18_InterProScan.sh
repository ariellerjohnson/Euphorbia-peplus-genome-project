#!/bin/bash

#SBATCH --job-name=InterProScan.sh
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=arj66@cornell.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=50gb
#SBATCH --time=999:00:00
#SBATCH --output=joblogs/InterProScan.joblog

date

#make directory for this step
mkdir 18_InterProScan
cd 18_InterProScan

#get InterProScan
wget http://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.55-88.0/interproscan-5.55-88.0-64-bit.tar.gz
tar xvfz interproscan-5.55-88.0-64-bit.tar.gz

#files used
protein_fasta=../15_TSEBRA/Euphorbia_peplus.aa

#Genome name
genome_name=Euphorbia_peplus

#copy proteins
cp $protein_fasta ${genome_name}.aa

#remove asterisks and blanks in proteins file
sed 's/*//g' ${genome_name}.aa > proteins_no_asterisk.aa
sed '/^$/d' proteins_no_asterisk.aa > proteins_no_asterisk_no_blanks.aa

#run InterProScan
interproscan-5.55-88.0/interproscan.sh -b ips_output -f XML -i proteins_no_asterisk_no_blanks.aa --goterms --pathways --iprlookup -t p -T ./

date
