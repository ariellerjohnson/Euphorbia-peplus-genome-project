#!/bin/bash

#SBATCH --job-name=3D_DNA_post_adjustment.sh
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=arj66@cornell.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=30
#SBATCH --mem=100gb
#SBATCH --time=999:00:00
#SBATCH --output=joblogs/3D_DNA_post_adjustment.joblog

date

#make directory for this step
mkdir 06_3D_DNA_post_adjustment
cd 06_3D_DNA_post_adjustment

#files used
post_review_assembly=../05_manual_adjustment_in_JBAT/E_peplus_hifiasm.final.review.review.assembly
hifiasm_assembly=/home/FrankLab/Arielle_E_peplus_genome_assembly/2_hifiasm_assembly/hifiasm_output.bp.p_ctg.fa
juicer_run=/home/FrankLab/Arielle_E_peplus_genome_assembly/8_juicer/old_attempts/untrimmed_hic_juicer/topDir/aligned/merged_nodups.txt

#export 3D-DNA path
export PATH=/programs/3d-dna:/programs/lastz-distrib-1.04:$PATH

#run post-review script to generate new fasta for genome
run-asm-pipeline-post-review.sh -r $post_review_assembly $hifiasm_assembly $juicer_run

date
