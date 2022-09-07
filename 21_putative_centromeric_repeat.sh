#!/bin/bash

#SBATCH --job-name=putative_centromeric_repeat.sh
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=arj66@cornell.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH --mem=50gb
#SBATCH --time=999:00:00
#SBATCH --output=joblogs/putative_centromeric_repeat.joblog

date

mkdir 21_putative_centromeric_repeat
cd 21_putative_centromeric_repeat

#files used
post_review_assembly=../06_3D_DNA_post_adjustment/E_peplus_hifiasm.bp.p_ctg.FINAL.fasta

#Genome name
genome_name=Euphorbia_peplus

#copy genome
cp $post_review_assembly ${genome_name}.fasta

#Run Tandem Repeats Finder (Cornell BioHPC already has it installed and in the path)
trf ${genome_name}.fasta 2 7 7 80 10 50 2000 -h

#get only the data lines
grep '^[0-9]' ${genome_name}.fasta.2.7.7.80.10.50.2000.dat > data_lines_only.dat
sort -k4 -nr data_lines_only.dat > sorted_data_lines_only.dat

#select only the top repeat
awk 'NR == 1 {print $14}' sorted_data_lines_only.dat > top_repeat.fa

#make a fasta header
sed -i '1 s/^/>top_repeat\n/' top_repeat.fa

#blast for the repeat in genome
blastn -query ${genome_name}.fasta -subject top_repeat.fa -evalue 1e-10 -outfmt 6 > top_repeat_blast.out

#convert blast results to gff
git clone https://github.com/raymondkiu/blastoutput2gff/
chmod +x blastoutput2gff/blastoutput2gff.sh
./blastoutput2gff/blastoutput2gff.sh top_repeat_blast.out

date
