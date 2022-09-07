#!/bin/bash

#SBATCH --job-name=RepeatModeler.sh
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=arj66@cornell.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH --mem=100gb
#SBATCH --time=999:00:00
#SBATCH --output=joblogs/RepeatModeler.joblog

date

mkdir 12_RepeatModeler
cd 12_RepeatModeler

#get TETools
singularity pull tetools.sif docker://dfam/tetools:latest

#files used
post_review_assembly=../06_3D_DNA_post_adjustment/E_peplus_hifiasm.bp.p_ctg.FINAL.fasta

#Genome name
genome_name=Euphorbia_peplus

#copy genome
cp $post_review_assembly $genome_name.fa

#build database
export LC_ALL=C
./tetools.sif BuildDatabase -name $genome_name $genome_name.fa

#model repeats
./tetools.sif RepeatModeler -database $genome_name -pa 40 -LTRStruct

date
