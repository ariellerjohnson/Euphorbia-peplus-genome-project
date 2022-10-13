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

#paths to files
files_directory="30_Kimura_distance"
peplus_genes=../22_finalizing_genome_for_upload/E_peplus_chromosomal_assembly/Euphorbia_peplus.gff
peplus_centro=../21_putative_centromeric_repeat/top_repeat_blast.out.gff

for species in "E_peplus" "E_lathyris" "H_brasiliensis" "M_esculenta" "R_communis"
do
  #copy the table and gff files
  cp ../${files_directory}/chromosome_only_assemblies/${species}.fa.table ./
  cp ../${files_directory}/chromosome_only_assemblies/${species}.fa.out.gff ./
  cp ../${files_directory}/chromosome_only_assemblies/${species}.fa ./
  #delete table header lines 1-6
  sed "1,6d" ${species}.fa.table > ${species}.fa.table.noheader
  #delete table footer lines starting with blank lines
  awk '/^$/{exit}1' ${species}.fa.table.noheader > ${species}.fa.table.fixed
  #delete extra stuff in gff file attributes
  awk 'BEGIN { OFS=FS="\t" }{ gsub(/^.*\"Motif:|\".*$/, "", $9) }1' ${species}.fa.out.gff > ${species}.fa.out.gff.fixed
  #make array using table and rename column 3 of gff based on the array
  awk 'BEGIN { OFS=FS="\t" } FNR==NR { a[$2]=$1; next } $9 in a { $3=a[$9] }1' ${species}.fa.table.fixed ${species}.fa.out.gff.fixed > ${species}.TEs.gff
  #make separate Ty3 and Copia files
  grep 'LTR/Copia' ${species}.TEs.gff > ${species}.Copia.gff
  grep 'LTR/Gypsy' ${species}.TEs.gff > ${species}.Ty3.gff
  #index the genome to prepare to make sliding windows
  samtools faidx ${species}.fa
  #make a text file in the format bedtools wants for -g (which is the chromosome name and its length)
  awk '{ print $1"\t"$2 }' ${species}.fa.fai > ${species}.lengths.bedtools.txt
  #make a text file for loading the chromosomes into karyoploteR
  awk '{ print $1"\t1\t"$2"\tNA\tNA" }' ${species}.fa.fai  > ${species}.lengths.karyoploter.txt
  #make 500k windows
  bedtools makewindows -g ${species}.lengths.bedtools.txt -w 500000 -s 500000 > ${species}.500k.windows.bed
  #calculate coverage
  bedtools coverage -b ${species}.TEs.gff -a ${species}.500k.windows.bed > ${species}.TEs.500kwindows.cov.txt
  bedtools coverage -b ${species}.Copia.gff -a ${species}.500k.windows.bed > ${species}.Copia.500kwindows.cov.txt
  bedtools coverage -b ${species}.Ty3.gff -a ${species}.500k.windows.bed > ${species}.Ty3.500kwindows.cov.txt

done

#calculate TE lengths to check why density and coverage are so different for annotated TEs
awk 'NR >= 4 { $10 = $5 - $4 + 1 } 1' E_peplus.TEs.gff > Ep_TE_lengths.gff
awk 'NR >= 4 { print $1"\t"$4"\t"$5"\t"$10 }' Ep_TE_lengths.gff > Ep_TE_lengths.bed
bedtools sort -i Ep_TE_lengths.bed > sorted_Ep_TE_lengths.bed
bedtools map -a E_peplus.500k.windows.bed -b sorted_Ep_TE_lengths.bed -c 4 -o mean > Ep_TE_lengths.mean.txt

#Peplus genes and centro repeats
cp $peplus_genes E_peplus_genes.gff
cp $peplus_centro E_peplus_centro.gff
sed -i 's/HiC_scaffold_/Chromosome_/g' E_peplus_centro.gff
bedtools coverage -b E_peplus_genes.gff -a E_peplus.500k.windows.bed > E_peplus.genes.500kwindows.cov.txt
bedtools coverage -b E_peplus_centro.gff -a E_peplus.500k.windows.bed > E_peplus.centro.500kwindows.cov.txt

date
