#!/bin/bash

#SBATCH --job-name=hifiasm.sh
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=arj66@cornell.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=30
#SBATCH --mem=100gb
#SBATCH --time=999:00:00
#SBATCH --output=joblogs/hifiasm.joblog

date

#make directory for this step
mkdir 01_hifiasm
cd 01_hifiasm

#files used
ccs_hifi_file=../raw_data/HiFi_CCS_reads/PBmixSequel991_2_B01_PCTM_30hours_19kbV2PD_70pM_Euphorbiapeplus.fastq.gz
outname=E_peplus

#installing hifiasm
git clone https://github.com/chhylp123/hifiasm
cd hifiasm && make
cd ..

#running hifiasm on the hifi CCS data
hifiasm/hifiasm -version
hifiasm/hifiasm -o $outname -t 30 $ccs_hifi_file

#converting to fasta format
awk '/^S/{print ">"$2;print $3}' *.bp.p_ctg.gfa > E_peplus_hifiasm.fa

date
