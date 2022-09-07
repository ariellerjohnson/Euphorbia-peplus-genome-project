#!/bin/bash

#SBATCH --job-name=3D_DNA.sh
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=arj66@cornell.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=30
#SBATCH --mem=100gb
#SBATCH --time=999:00:00
#SBATCH --output=joblogs/3D_DNA.joblog

date

#make directory for this step
mkdir 04_3D_DNA
cd 04_3D_DNA

#files used
juicer_run=../03_juicer
hifiasm_assembly=../01_hifiasm/E_peplus_hifiasm.fa

#export path to executable
export PATH=/programs/3d-dna:/programs/lastz-distrib-1.04:$PATH

#run 3D-DNA
run-asm-pipeline.sh --editor-saturation-centile 10 --editor-coarse-resolution 250000 --editor-coarse-region 400000 --editor-repeat-coverage 50 $hifiasm_assembly ${juicer_run}/topDir/aligned/merged_nodups.txt

date
