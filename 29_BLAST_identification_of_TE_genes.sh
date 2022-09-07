#!/bin/bash

#SBATCH --job-name=BLAST_identification_of_TE_genes.sh
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=arj66@cornell.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=50gb
#SBATCH --time=999:00:00
#SBATCH --output=joblogs/BLAST_identification_of_TE_genes.joblog

date

mkdir 29_BLAST_identification_of_TE_genes
cd 29_BLAST_identification_of_TE_genes

#copy fasta files for genes of interest
cp -r ../raw_data/genes_of_interest/TE_regulation_genes_of_interest/* ./

#copy proteins
cp ../22_finalizing_genome_for_upload/E_peplus_chromosomal_assembly/Euphorbia_peplus.aa ./

cd AGO1
blastp -query ../Euphorbia_peplus.aa -subject Arabidopsis_AGO1.fasta -qcov_hsp_perc 80 \
-evalue 1e-10 -outfmt "6 qcovs qseqid sseqid pident length mismatch gapopen evalue" > putative_AGO1_blast.out
sort -k 4nr putative_AGO1_blast.out > sorted_putative_AGO1_blast.out
cd ..

cd AGO9
blastp -query ../Euphorbia_peplus.aa -subject Arabidopsis_AGO9.fasta -qcov_hsp_perc 80 \
-evalue 1e-10 -outfmt "6 qcovs qseqid sseqid pident length mismatch gapopen evalue" > putative_AGO9_blast.out
sort -k 4nr putative_AGO9_blast.out > sorted_putative_AGO9_blast.out
cd ..

cd ATXR5
blastp -query ../Euphorbia_peplus.aa -subject Arabidopsis_ATXR5.fasta -qcov_hsp_perc 80 \
-evalue 1e-10 -outfmt "6 qcovs qseqid sseqid pident length mismatch gapopen evalue" > putative_ATXR5_blast.out
sort -k 4nr putative_ATXR5_blast.out > sorted_putative_ATXR5_blast.out
cd ..

cd ATXR6
blastp -query ../Euphorbia_peplus.aa -subject Arabidopsis_ATXR6.fasta -qcov_hsp_perc 80 \
-evalue 1e-10 -outfmt "6 qcovs qseqid sseqid pident length mismatch gapopen evalue" > putative_ATXR6_blast.out
sort -k 4nr putative_ATXR6_blast.out > sorted_putative_ATXR6_blast.out
cd ..

cd CMT2
blastp -query ../Euphorbia_peplus.aa -subject Arabidopsis_CMT2.fasta -qcov_hsp_perc 80 \
-evalue 1e-10 -outfmt "6 qcovs qseqid sseqid pident length mismatch gapopen evalue" > putative_CMT2_blast.out
sort -k 4nr putative_CMT2_blast.out > sorted_putative_CMT2_blast.out
cd ..

cd CMT3
blastp -query ../Euphorbia_peplus.aa -subject Arabidopsis_CMT3.fasta -qcov_hsp_perc 80 \
-evalue 1e-10 -outfmt "6 qcovs qseqid sseqid pident length mismatch gapopen evalue" > putative_CMT3_blast.out
sort -k 4nr putative_CMT3_blast.out > sorted_putative_CMT3_blast.out
cd ..

cd DDM1
blastp -query ../Euphorbia_peplus.aa -subject Arabidopsis_DDM1.fasta -qcov_hsp_perc 80 \
-evalue 1e-10 -outfmt "6 qcovs qseqid sseqid pident length mismatch gapopen evalue" > putative_DDM1_blast.out
sort -k 4nr putative_DDM1_blast.out > sorted_putative_DDM1_blast.out
cd ..

cd DRM1
blastp -query ../Euphorbia_peplus.aa -subject Arabidopsis_DRM1.fasta -qcov_hsp_perc 80 \
-evalue 1e-10 -outfmt "6 qcovs qseqid sseqid pident length mismatch gapopen evalue" > putative_DRM1_blast.out
sort -k 4nr putative_DRM1_blast.out > sorted_putative_DRM1_blast.out
cd ..

cd DRM2
blastp -query ../Euphorbia_peplus.aa -subject Arabidopsis_DRM2.fasta -qcov_hsp_perc 80 \
-evalue 1e-10 -outfmt "6 qcovs qseqid sseqid pident length mismatch gapopen evalue" > putative_DRM2_blast.out
sort -k 4nr putative_DRM2_blast.out > sorted_putative_DRM2_blast.out
cd ..

cd MET1
blastp -query ../Euphorbia_peplus.aa -subject Arabidopsis_MET1.fasta -qcov_hsp_perc 80 \
-evalue 1e-10 -outfmt "6 qcovs qseqid sseqid pident length mismatch gapopen evalue" > putative_MET1_blast.out
sort -k 4nr putative_MET1_blast.out > sorted_putative_MET1_blast.out
cd ..

cd MET2a
blastp -query ../Euphorbia_peplus.aa -subject Arabidopsis_MET2a.fasta -qcov_hsp_perc 80 \
-evalue 1e-10 -outfmt "6 qcovs qseqid sseqid pident length mismatch gapopen evalue" > putative_MET2a_blast.out
sort -k 4nr putative_MET2a_blast.out > sorted_putative_MET2a_blast.out
cd ..

cd NRPD1
blastp -query ../Euphorbia_peplus.aa -subject Arabidopsis_NRPD1.fasta -qcov_hsp_perc 80 \
-evalue 1e-10 -outfmt "6 qcovs qseqid sseqid pident length mismatch gapopen evalue" > putative_NRPD1_blast.out
sort -k 4nr putative_NRPD1_blast.out > sorted_putative_NRPD1_blast.out
cd ..

cd NRPE1
blastp -query ../Euphorbia_peplus.aa -subject Arabidopsis_NRPE1.fasta -qcov_hsp_perc 80 \
-evalue 1e-10 -outfmt "6 qcovs qseqid sseqid pident length mismatch gapopen evalue" > putative_NRPE1_blast.out
sort -k 4nr putative_NRPE1_blast.out > sorted_putative_NRPE1_blast.out
cd ..

cd NRPE4
blastp -query ../Euphorbia_peplus.aa -subject Arabidopsis_NRPE4.fasta -qcov_hsp_perc 80 \
-evalue 1e-10 -outfmt "6 qcovs qseqid sseqid pident length mismatch gapopen evalue" > putative_NRPE4_blast.out
sort -k 4nr putative_NRPE4_blast.out > sorted_putative_NRPE4_blast.out
cd ..

cd NRPE5
blastp -query ../Euphorbia_peplus.aa -subject Arabidopsis_NRPE5.fasta -qcov_hsp_perc 80 \
-evalue 1e-10 -outfmt "6 qcovs qseqid sseqid pident length mismatch gapopen evalue" > putative_NRPE5_blast.out
sort -k 4nr putative_NRPE5_blast.out > sorted_putative_NRPE5_blast.out
cd ..

cd RDR6
blastp -query ../Euphorbia_peplus.aa -subject Arabidopsis_RDR6.fasta -qcov_hsp_perc 80 \
-evalue 1e-10 -outfmt "6 qcovs qseqid sseqid pident length mismatch gapopen evalue" > putative_RDR6_blast.out
sort -k 4nr putative_RDR6_blast.out > sorted_putative_RDR6_blast.out
cd ..

date
