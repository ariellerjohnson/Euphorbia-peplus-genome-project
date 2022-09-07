#!/bin/bash

#SBATCH --job-name=finalizing_genome_for_upload.sh
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=arj66@cornell.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=50gb
#SBATCH --time=999:00:00
#SBATCH --output=joblogs/finalizing_genome_for_upload.joblog

date

#files used
gff_file=../15_TSEBRA/Euphorbia_peplus.gff
protein_file=../18_InterProScan/proteins_no_asterisk_no_blanks.aa
coding_seq_file=../15_TSEBRA/Euphorbia_peplus.codingseq
genome_file=../06_3D_DNA_post_adjustment/E_peplus_hifiasm.bp.p_ctg.FINAL.fasta
genome_with_zeroes=../15_TSEBRA/Euphorbia_peplus.fa

#seqtk path
seqtk_path=/programs/seqtk/seqtk

#make directory for this step
mkdir 22_finalizing_genome_for_upload
cd 22_finalizing_genome_for_upload

#make subdirectories for the two collections of files
mkdir E_peplus_chromosomal_assembly
mkdir E_peplus_assembly_including_all_scaffolds


###GFF file cleanup
cp $gff_file gff_file.gff
cp $gff_file E_peplus_assembly_including_all_scaffolds/tentative.gff

output_name=Euphorbia_peplus

cp gff_file.gff $output_name

chrom_numbers=$(seq 8)

for number in ${chrom_numbers[*]}
do
  HiC="HiC_scaffold_000"${number}
  Chromosome="Chromosome_"${number}
  sc="_sc000"${number}
  chr="_chr"${number}
  sed -i "s/${HiC}/${Chromosome}/g" $output_name
  sed -i "s/${sc}/${chr}/g" $output_name
done

sed -i '/HiC_scaffold/d' $output_name

uniq $output_name > $output_name.gff

rm $output_name

cp $output_name.gff E_peplus_chromosomal_assembly/tentative.gff

rm $output_name.gff
rm gff_file.gff

###Protein cleanup
cp $protein_file protein_file.fasta
cp $protein_file E_peplus_assembly_including_all_scaffolds/tentative.aa

python3 ../22_1_selecting_chromosomes_and_removing_zeroes.py protein_file.fasta
cp protein_file.fasta_chromosomes.fasta E_peplus_chromosomal_assembly/tentative.aa

rm protein_file.fasta_chromosomes.fasta
rm protein_file.fasta

###Coding sequence cleanup
cp $coding_seq_file coding_seq.fasta
cp $coding_seq_file E_peplus_assembly_including_all_scaffolds/tentative.codingseq

python3 ../22_1_selecting_chromosomes_and_removing_zeroes.py coding_seq.fasta
cp coding_seq.fasta_chromosomes.fasta E_peplus_chromosomal_assembly/tentative.codingseq

rm coding_seq.fasta_chromosomes.fasta
rm coding_seq.fasta

###Genome cleanup
cp $genome_file genome_file.fasta

$seqtk_path subseq genome_file.fasta ../22_2_chromosome_names.txt > chromosomes_only.fasta

sed -i 's/HiC_scaffold/Chromosome/g' chromosomes_only.fasta

mv chromosomes_only.fasta E_peplus_chromosomal_assembly/${output_name}.fa

cp $genome_with_zeroes E_peplus_assembly_including_all_scaffolds/E_peplus_all_scaffolds.fa

rm genome_file.fasta

#cleaning out fake gene in E_peplus_chromosomal_assembly
cd E_peplus_chromosomal_assembly
blastp -query tentative.aa -subject ../../22_3_potential_fake_gene.fasta -evalue 1e-10 -outfmt 6 > fake_blast.out
awk '{print $1}' fake_blast.out > fake_genes.txt
$seqtk_path subseq tentative.aa fake_genes.txt > removed_genes.fasta
sed -i 's/.t1//g' fake_genes.txt
grep -Fvf fake_genes.txt tentative.gff > tentative_no_fake.gff
uniq tentative_no_fake.gff ${output_name}.gff
cat tentative.aa | paste - - - - | grep -Fvf fake_genes.txt | tr "\t" "\n" > ${output_name}.aa
cat tentative.codingseq | paste - - - - | grep -Fvf fake_genes.txt | tr "\t" "\n" > ${output_name}.codingseq
rm *tentative*
rm *fake*
cd ..

#cleaning out fake gene in E_peplus_assembly_including_all_scaffolds
cd E_peplus_assembly_including_all_scaffolds
blastp -query tentative.aa -subject ../../22_3_potential_fake_gene.fasta -evalue 1e-10 -outfmt 6 > fake_blast.out
awk '{print $1}' fake_blast.out > fake_genes.txt
$seqtk_path subseq tentative.aa fake_genes.txt > removed_genes.fasta
sed -i 's/.t1//g' fake_genes.txt
grep -Fvf fake_genes.txt tentative.gff > tentative_no_fake.gff
uniq tentative_no_fake.gff E_peplus_all_scaffolds.gff
cat tentative.aa | paste - - - - | grep -Fvf fake_genes.txt | tr "\t" "\n" > E_peplus_all_scaffolds.aa
cat tentative.codingseq | paste - - - - | grep -Fvf fake_genes.txt | tr "\t" "\n" > E_peplus_all_scaffolds.codingseq
rm *tentative*
rm *fake*

#(added later) editing names in E_peplus_assembly_including_all_scaffolds to include "chromosomes" also
chrom_numbers=$(seq 8)
for number in ${chrom_numbers[*]}
do
  HiC="HiC_scaffold_000"${number}
  Chromosome="Chromosome_"${number}
  sc="_sc000"${number}
  chr="_chr"${number}
  sed -i "s/${HiC}/${Chromosome}/g" *.*
  sed -i "s/${sc}/${chr}/g" *.*
done

date
