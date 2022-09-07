#!/bin/bash

#SBATCH --job-name=HiC_QC.sh
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=arj66@cornell.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=30
#SBATCH --mem=100gb
#SBATCH --time=999:00:00
#SBATCH --output=joblogs/HiC_QC.joblog

date

#make directory for this step
mkdir 02_HiC_QC
cd 02_HiC_QC

#files used
raw_HiC_R1=../raw_data/HiC_data/12221_11073_130590_HKT53AFX2_E_peplus_Hi_C_ATGTGGAC_TTGGTGAG_R1.fastq.gz
raw_HiC_R2=../raw_data/HiC_data/12221_11073_130590_HKT53AFX2_E_peplus_Hi_C_ATGTGGAC_TTGGTGAG_R2.fastq.gz
hifiasm_assembly=../01_hifiasm/E_peplus_hifiasm.fa

#export paths to executables
export PATH=/programs/samblaster:$PATH
#bwa is already in path on BioHPC at Cornell

#download necessary programs and scripts
git clone https://github.com/phasegenomics/hic_qc.git #instructions at https://github.com/phasegenomics/hic_qc
cd hic_qc/
source $HOME/miniconda3/bin/activate #this activates conda on BioHPC at Cornell
conda env create -n hic_qc --file env.yml
conda activate hic_qc
python setup.py install --user
cd ..

#aligning Hi-C data with hifiasm assembly
bwa index $hifiasm_assembly
bwa mem -5SP $hifiasm_assembly $raw_HiC_R1 $raw_HiC_R2 -t 30 | samblaster | samtools view -S -h -b -F 2316 > hic_aligned_with_hifiasm.bam

#running QC python script
hic_qc/hic_qc.py -b hic_aligned_with_hifiasm.bam --outfile_prefix hic_qc_aligned_with_hifiasm

date
