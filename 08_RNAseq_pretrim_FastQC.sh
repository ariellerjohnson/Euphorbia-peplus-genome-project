#!/bin/bash

#SBATCH --job-name=RNAseq_pretrim_FastQC.sh
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=arj66@cornell.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=30
#SBATCH --mem=50gb
#SBATCH --time=999:00:00
#SBATCH --output=joblogs/RNAseq_pretrim_FastQC.joblog

date

#make directory for this step
mkdir 08_RNAseq_pretrim_FastQC
cd 08_RNAseq_pretrim_FastQC

#files used
RNAseq_raw_dir=../raw_data/RNAseq_data

#export executable to path
export PATH=/programs/FastQC-0.11.8:$PATH

#other variables
outdir=$(pwd)

fastqc --threads 30 $RNAseq_raw_dir/*.fastq.gz --outdir $outdir

date
