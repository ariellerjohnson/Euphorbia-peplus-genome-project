#!/bin/bash

#SBATCH --job-name=kimura_distance.sh
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=arj66@cornell.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH --mem=75gb
#SBATCH --time=999:00:00
#SBATCH --output=joblogs/kimura_distance.joblog

#date

mkdir 30_Kimura_distance
cd 30_Kimura_distance

#get TETools
singularity pull tetools.sif docker://dfam/tetools:latest

#seqtk path
seqtk_path=/programs/seqtk/seqtk

#paths to peplus files
peplus_unmasked_assembly=../22_finalizing_genome_for_upload/E_peplus_chromosomal_assembly/Euphorbia_peplus.fa
assembly_dir=../raw_data/assemblies_for_Kimura

mkdir chromosome_only_assemblies

#get Euphorbia peplus chromosome only assembly
cp $peplus_unmasked_assembly chromosome_only_assemblies/E_peplus.fa

#get assemblies from other species
cp $assembly_dir/GCA_010458925.1_ASM1045892v1_genomic.fna chromosome_only_assemblies/H_brasiliensis_all.fa
cp $assembly_dir/Wild_castor_genome.fa chromosome_only_assemblies/R_communis_all.fa
cp $assembly_dir/Mesculenta_671_v8.0.fa chromosome_only_assemblies/M_esculenta_all.fa
cp $assembly_dir/genome.fasta chromosome_only_assemblies/E_lathyris_all.fa

#select chromosomes
cd chromosome_only_assemblies
$seqtk_path subseq E_lathyris_all.fa ../../30_1_E_lathyris_chrom.txt > E_lathyris.fa
$seqtk_path subseq H_brasiliensis_all.fa ../../30_2_H_brasiliensis_chrom.txt > H_brasiliensis.fa
$seqtk_path subseq M_esculenta_all.fa ../../30_3_M_esculenta_chrom.txt > M_esculenta.fa
$seqtk_path subseq R_communis_all.fa ../../30_4_R_communis_chrom.txt > R_communis.fa

rm *_all.fa

#"undo" softmasking of Hevea genome by making everything uppercase
sed -i 's/[a-z]/\U&/g' H_brasiliensis.fa

#build databases
cd ..
export LC_ALL=C
./tetools.sif BuildDatabase -name E_peplus chromosome_only_assemblies/E_peplus.fa
./tetools.sif BuildDatabase -name E_lathyris chromosome_only_assemblies/E_lathyris.fa
./tetools.sif BuildDatabase -name H_brasiliensis chromosome_only_assemblies/H_brasiliensis.fa
./tetools.sif BuildDatabase -name M_esculenta chromosome_only_assemblies/M_esculenta.fa
./tetools.sif BuildDatabase -name R_communis chromosome_only_assemblies/R_communis.fa

#model repeats
./tetools.sif RepeatModeler -database E_peplus -pa 40 -LTRStruct
./tetools.sif RepeatModeler -database E_lathyris -pa 40 -LTRStruct
./tetools.sif RepeatModeler -database H_brasiliensis -pa 40 -LTRStruct
./tetools.sif RepeatModeler -database M_esculenta -pa 40 -LTRStruct
./tetools.sif RepeatModeler -database R_communis -pa 40 -LTRStruct

#mask repeats
./tetools.sif RepeatMasker -a -lib E_peplus-families.fa -pa 40 -gff -xsmall chromosome_only_assemblies/E_peplus.fa
./tetools.sif RepeatMasker -a -lib E_lathyris-families.fa -pa 40 -gff -xsmall chromosome_only_assemblies/E_lathyris.fa
./tetools.sif RepeatMasker -a -lib H_brasiliensis-families.fa -pa 40 -gff -xsmall chromosome_only_assemblies/H_brasiliensis.fa
./tetools.sif RepeatMasker -a -lib M_esculenta-families.fa -pa 40 -gff -xsmall chromosome_only_assemblies/M_esculenta.fa
./tetools.sif RepeatMasker -a -lib R_communis-families.fa -pa 40 -gff -xsmall chromosome_only_assemblies/R_communis.fa

#get Kimura matrix
./tetools.sif calcDivergenceFromAlign.pl -s chromosome_only_assemblies/E_peplus.fa.table -a chromosome_only_assemblies/new_name1.align chromosome_only_assemblies/E_peplus.fa.align
./tetools.sif calcDivergenceFromAlign.pl -s chromosome_only_assemblies/E_lathyris.fa.table -a chromosome_only_assemblies/new_name2.align chromosome_only_assemblies/E_lathyris.fa.align
./tetools.sif calcDivergenceFromAlign.pl -s chromosome_only_assemblies/H_brasiliensis.fa.table -a chromosome_only_assemblies/new_name3.align chromosome_only_assemblies/H_brasiliensis.fa.align
./tetools.sif calcDivergenceFromAlign.pl -s chromosome_only_assemblies/M_esculenta.fa.table -a chromosome_only_assemblies/new_name4.align chromosome_only_assemblies/M_esculenta.fa.align
./tetools.sif calcDivergenceFromAlign.pl -s chromosome_only_assemblies/R_communis.fa.table -a chromosome_only_assemblies/new_name5.align chromosome_only_assemblies/R_communis.fa.align

cd chromosome_only_assemblies
sed -e '1,/Coverage for each repeat class and divergence (Kimura)/ d' E_peplus.fa.table > E_peplus_Kimura_table.txt
sed -e '1,/Coverage for each repeat class and divergence (Kimura)/ d' E_lathyris.fa.table > E_lathyris_Kimura_table.txt
sed -e '1,/Coverage for each repeat class and divergence (Kimura)/ d' H_brasiliensis.fa.table > H_brasiliensis_Kimura_table.txt
sed -e '1,/Coverage for each repeat class and divergence (Kimura)/ d' M_esculenta.fa.table > M_esculenta_Kimura_table.txt
sed -e '1,/Coverage for each repeat class and divergence (Kimura)/ d' R_communis.fa.table > R_communis_Kimura_table.txt

date
