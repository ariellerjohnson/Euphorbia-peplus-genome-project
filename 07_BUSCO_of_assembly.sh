#!/bin/bash

#SBATCH --job-name=BUSCO_of_assembly.sh
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=arj66@cornell.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=30
#SBATCH --mem=50gb
#SBATCH --time=999:00:00
#SBATCH --output=joblogs/BUSCO_of_assembly.joblog

date

#make directory for this step
mkdir 07_BUSCO_of_assembly
cd 07_BUSCO_of_assembly

#files used
genome_fasta=../06_3D_DNA_post_adjustment/E_peplus_hifiasm.bp.p_ctg.FINAL.fasta

#other variables
busco_lineage=embryophyta_odb10

#activate busco conda environment (already installed on BioHPC at Cornell)
source /programs/miniconda3/bin/activate busco-5.1.2

busco -i $genome_fasta -l $busco_lineage -o genome_BUSCO_${busco_lineage} -m genome --cpu 30

date
