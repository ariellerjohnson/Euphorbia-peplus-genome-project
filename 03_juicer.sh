#!/bin/bash

#SBATCH --job-name=juicer.sh
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=arj66@cornell.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=30
#SBATCH --mem=100gb
#SBATCH --time=999:00:00
#SBATCH --output=joblogs/juicer.joblog

date

#make directory for this step
mkdir 03_juicer
cd 03_juicer

#files used
raw_HiC_R1=../../raw_data/HiC_data/12221_11073_130590_HKT53AFX2_E_peplus_Hi_C_ATGTGGAC_TTGGTGAG_R1.fastq.gz
raw_HiC_R2=../../raw_data/HiC_data/12221_11073_130590_HKT53AFX2_E_peplus_Hi_C_ATGTGGAC_TTGGTGAG_R2.fastq.gz
hifiasm_assembly=../../01_hifiasm/E_peplus_hifiasm.fa

#other variables
enzyme=Sau3AI
outname=E_peplus_hifiasm
site=MboI

#paths to executable (as a variable here)
juicer_location=/programs/juicer-1.6

#export paths to executables
export PATH=/programs/samblaster:$PATH

#setting up juicer directories
#(see https://github.com/aidenlab/juicer/wiki/Installation#quick-start)
mkdir topDir
cd topDir

mkdir fastq
cp $raw_HiC_R1 fastq/
cp $raw_HiC_R2 fastq/

cd ..
mkdir juicerDir
cd juicerDir

mkdir references
cp $hifiasm_assembly references/
cd references
bwa index *.fa #this indexes the genome, you could also copy from the HiC_QC
cd ..

ln -s ${juicer_location}/CPU scripts

mkdir restriction_sites
cd restriction_sites

# running script to generate restriction sites
python ${juicer_location}/misc/generate_site_positions.py $enzyme $outname ../references/*.fa

#getting contig sizes
awk 'BEGIN{OFS="\t"}{print $1, $NF}' ${outname}_${enzyme}.txt > ${outname}.contig.sizes

cd ../..

#working directory
myworkdir=$(pwd)

#running juicer
$myworkdir/juicerDir/scripts/juicer.sh -D $myworkdir/juicerDir -d $myworkdir/topDir -y $myworkdir/juicerDir/restriction_sites/${outname}_${enzyme}.txt -g "${outname}" -s $site -z $myworkdir/juicerDir/references/*.fa -p $myworkdir/juicerDir/restriction_sites/${outname}.contig.sizes

date
