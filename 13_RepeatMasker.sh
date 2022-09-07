#!/bin/bash

#SBATCH --job-name=RepeatMasker.sh
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=arj66@cornell.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH --mem=100gb
#SBATCH --time=999:00:00
#SBATCH --output=joblogs/RepeatMasker.joblog

date

#make directory for this step
mkdir 13_RepeatMasker
cd 13_RepeatMasker

#files used
post_review_assembly=../06_3D_DNA_post_adjustment/E_peplus_hifiasm.bp.p_ctg.FINAL.fasta
library_file=../12_RepeatModeler/Euphorbia_peplus-families.fa

#Genome name
genome_name=Euphorbia_peplus

#copy genome and library file
cp $post_review_assembly $genome_name.fa
cp $library_file $genome_name-families.fa

#directory containing TETools
TETools_location=../12_RepeatModeler

export LC_ALL=C
$TETools_location/tetools.sif RepeatMasker -lib $genome_name-families.fa \
-pa 40 -gff -xsmall -nolow $genome_name.fa

date
