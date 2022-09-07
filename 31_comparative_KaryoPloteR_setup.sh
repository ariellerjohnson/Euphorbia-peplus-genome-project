#!/bin/bash

#SBATCH --job-name=comparative_KaryoPloteR_setup.sh
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=arj66@cornell.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH --mem=75gb
#SBATCH --time=999:00:00
#SBATCH --output=joblogs/comparative_KaryoPloteR_setup.joblog

date

mkdir 31_comparative_KaryoPloteR_setup
cd 31_comparative_KaryoPloteR_setup

files_directory="30_Kimura_distance"

for species in "E_peplus" "E_lathyris" "H_brasiliensis" "M_esculenta" "R_communis"
do
  #copy the table and gff files
  cp ../${files_directory}/chromosome_only_assemblies/${species}.fa.table ./
  cp ../${files_directory}/chromosome_only_assemblies/${species}.fa.out.gff ./
  #delete table header lines 1-6
  sed "1,6d" ${species}.fa.table > ${species}.fa.table.noheader
  #delete table footer lines starting with blank lines
  awk '/^$/{exit}1' ${species}.fa.table.noheader > ${species}.fa.table.fixed
  #delete extra stuff in gff file attributes
  awk 'BEGIN { OFS=FS="\t" }{ gsub(/^.*\"Motif:|\".*$/, "", $9) }1' ${species}.fa.out.gff > ${species}.fa.out.gff.fixed
  #make array using table and rename column 3 of gff based on the array
  awk 'BEGIN { OFS=FS="\t" } FNR==NR { a[$2]=$1; next } $9 in a { $3=a[$9] }1' ${species}.fa.table.fixed ${species}.fa.out.gff.fixed > ${species}.KaryoPlotR.gff
done

date
