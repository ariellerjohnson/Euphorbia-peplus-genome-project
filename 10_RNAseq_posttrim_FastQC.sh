#!/bin/bash

#SBATCH --job-name=RNAseq_posttrim_FastQC.sh
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=arj66@cornell.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=30
#SBATCH --mem=50gb
#SBATCH --time=999:00:00
#SBATCH --output=joblogs/RNAseq_posttrim_FastQC.joblog

date

#make directory for this step
mkdir 10_RNAseq_posttrim_FastQC
cd 10_RNAseq_posttrim_FastQC

#files used
RNAseq_trim_dir=../09_trimming_RNAseq

#export executable to path
export PATH=/programs/FastQC-0.11.8:$PATH

#other variables
outdir=$(pwd)

fastqc --threads 30 $RNAseq_trim_dir/*P --outdir $outdir
fastqc --threads 30 $RNAseq_trim_dir/*U --outdir $outdir

date
