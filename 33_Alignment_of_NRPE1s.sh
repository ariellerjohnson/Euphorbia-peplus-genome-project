#!/bin/bash

#SBATCH --job-name=Alignment_of_NRPE1s.sh
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=arj66@cornell.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=50gb
#SBATCH --time=999:00:00
#SBATCH --output=joblogs/Alignment_of_NRPE1s.joblog

date

mkdir 33_Alignment_of_NRPE1s
cd 33_Alignment_of_NRPE1s

#export MAFFT location
export PATH=/programs/mafft/bin:$PATH

#retrieve genes of interest
cp ../25_OrthoFinder/2_Euphorbiaceae/OrthoFinder/Results_May20/Orthogroup_Sequences/OG0005416.fa ./

#do mafft alignment
params="--maxiterate 1000 --localpair"

sequence_file="OG0005416.fa"

outfile_name=$(echo $sequence_file | cut -d '.' -f1 | awk '{print $0"_aligned.fa"}')

mafft $params  $sequence_file > $outfile_name

date
