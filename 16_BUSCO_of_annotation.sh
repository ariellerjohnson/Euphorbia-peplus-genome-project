#!/bin/bash

#SBATCH --job-name=BUSCO_of_annotation.sh
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=arj66@cornell.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=30
#SBATCH --mem=50gb
#SBATCH --time=999:00:00
#SBATCH --output=joblogs/BUSCO_of_annotation.joblog

date

#make directory for this step
mkdir 16_BUSCO_of_annotation
cd 16_BUSCO_of_annotation

#files used
protein_fasta=../15_TSEBRA/Euphorbia_peplus.aa

#other variables
busco_lineage=embryophyta_odb10

#activate busco conda environment (already installed on BioHPC at Cornell)
source /programs/miniconda3/bin/activate busco-5.1.2

busco -i $protein_fasta -l $busco_lineage -o protein_BUSCO_${busco_lineage} -m protein --cpu 30

date
