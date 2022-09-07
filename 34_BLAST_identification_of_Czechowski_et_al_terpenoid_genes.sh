#!/bin/bash

#SBATCH --job-name=BLAST_identification_of_Czechowski_et_al_terpenoid_genes.sh
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=arj66@cornell.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=50gb
#SBATCH --time=999:00:00
#SBATCH --output=joblogs/BLAST_identification_of_Czechowski_et_al_terpenoid_genes.joblog

date

mkdir 34_BLAST_identification_of_Czechowski_et_al_terpenoid_genes
cd 34_BLAST_identification_of_Czechowski_et_al_terpenoid_genes

#copy fasta files for genes of interest
cp -r ../raw_data/genes_of_interest/Czechowski_et_al_terpenoid_genes_of_interest/* ./

directories=$(echo *)

#copy proteins
cp ../22_finalizing_genome_for_upload/E_peplus_chromosomal_assembly/Euphorbia_peplus.aa ./

for directory in ${directories[*]}
do
  cd ${directory}
  outfile_name=$(echo ${directory} |  awk '{print "putative_"$0"_blast.out"}')
  tblastn -query ../Euphorbia_peplus.aa -subject ./*.fasta -qcov_hsp_perc 80 \
  -evalue 1e-10 -outfmt "6 qcovs qseqid sseqid pident length mismatch gapopen evalue" > ${outfile_name}
  sort -k 4nr ${outfile_name} > sorted_${outfile_name}
  cd ..
done

date
