#!/bin/bash

#SBATCH --job-name=TSEBRA.sh
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=arj66@cornell.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH --mem=50gb
#SBATCH --time=999:00:00
#SBATCH --output=joblogs/TSEBRA.joblog

date

#files used
post_review_assembly=../06_3D_DNA_post_adjustment/E_peplus_hifiasm.bp.p_ctg.FINAL.fasta
RNA_only_BRAKER_dir=../14_BRAKER/RNA_only_BRAKER/braker
protein_only_BRAKER_dir=../14_BRAKER/protein_only_BRAKER/braker

#Genome name
genome_name=Euphorbia_peplus

#add genometools to path
export PATH=/programs/genometools-1.5.9/bin:$PATH

#make directory for this step
mkdir 15_TSEBRA
cd 15_TSEBRA


git clone https://github.com/Gaius-Augustus/TSEBRA

./TSEBRA/bin/fix_gtf_ids.py --gtf ${RNA_only_BRAKER_dir}/braker.gtf --out braker1_fixed.gtf
./TSEBRA/bin/fix_gtf_ids.py --gtf ${protein_only_BRAKER_dir}/braker.gtf --out braker2_fixed.gtf

./TSEBRA/bin/tsebra.py -g braker1_fixed.gtf,braker2_fixed.gtf -c TSEBRA/config/default.cfg \
-e ${RNA_only_BRAKER_dir}/hintsfile.gff,${protein_only_BRAKER_dir}/hintsfile.gff \
-o braker1+2_combined.gtf

#see https://github.com/Gaius-Augustus/BRAKER/issues/194
# Sort braker gtf file
sort -k1,1 -k4,4n -k5,5n braker1+2_combined.gtf > braker1+2_combined_sorted.gtf

#make the gene and transcript names the same
./TSEBRA/bin/rename_gtf.py --gtf braker1+2_combined_sorted.gtf --prefix Ep \
--translation_tab rename_gtf_translation_tab.txt --out braker1+2_combined_sorted_renamed.gtf

#adding the scaffold name to each gene and transcript name
python3 ../15_1_add_scaf_name_to_gtf.py braker1+2_combined_sorted_renamed.gtf

cp braker1+2_combined_sorted_renamed_withScaf.gtf ${genome_name}.gtf

#get genome
cp $post_review_assembly genome_no_zeroes.fa

python3 ../15_2_add_zeroes_to_genome.py genome_no_zeroes.fa $genome_name.fa

# Get coding sequences and translation
/programs/Augustus-3.3.3/scripts/getAnnoFastaFromJoingenes.py \
-g $genome_name.fa \
-o $genome_name -f ${genome_name}.gtf
#produces .aa and .codingseq files

#convert to gff3
#see https://github.com/Gaius-Augustus/BRAKER/issues/275
gt gtf_to_gff3 <(grep -P "\tCDS\t|\texon\t" Euphorbia_peplus.gtf ) > unsorted.gff
gt gff3 -sort yes unsorted.gff > ${genome_name}.gff

date
