#!/bin/bash

#SBATCH --job-name=BLAST_identification_of_biosynth_genes.sh
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=arj66@cornell.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=50gb
#SBATCH --time=999:00:00
#SBATCH --output=joblogs/BLAST_identification_of_biosynth_genes.joblog

date

mkdir 28_BLAST_identification_of_biosynth_genes
cd 28_BLAST_identification_of_biosynth_genes

#copy fasta files for genes of interest
cp -r /workdir/arj66/EUPHORBIA_PEPLUS_GENOME/raw_data/genes_of_interest/terpenoid_biosynth_genes_of_interest/* ./

#copy proteins
cp ../22_finalizing_genome_for_upload/E_peplus_chromosomal_assembly/Euphorbia_peplus.aa ./

cd alcohol_dehydrogenase
blastp -query ../Euphorbia_peplus.aa -subject E_lathyris_ADH1.fasta -qcov_hsp_perc 80 \
-evalue 1e-10 -outfmt "6 qcovs qseqid sseqid pident length mismatch gapopen evalue" > putative_ADH1_lathyris_blast.out
sort -k 4nr putative_ADH1_lathyris_blast.out > sorted_putative_ADH1_lathyris_blast.out
tblastn -query ../Euphorbia_peplus.aa -subject E_peplus_ADH1.fasta -qcov_hsp_perc 80 \
-evalue 1e-10 -outfmt "6 qcovs qseqid sseqid pident length mismatch gapopen evalue" > putative_ADH1_peplus_blast.out
sort -k 4nr putative_ADH1_peplus_blast.out > sorted_putative_ADH1_peplus_blast.out
cd ..

cd casbene_synthase
blastp -query ../Euphorbia_peplus.aa -subject Ricinus_casbene_synthase.faa -qcov_hsp_perc 80 \
-evalue 1e-10 -outfmt "6 qcovs qseqid sseqid pident length mismatch gapopen evalue" > putative_casbene_synthase_Ricinus_blast.out
sort -k 4nr putative_casbene_synthase_Ricinus_blast.out > sorted_putative_casbene_synthase_Ricinus_blast.out
blastp -query ../Euphorbia_peplus.aa -subject E_lathyris_casbene_synthase.fasta -qcov_hsp_perc 80 \
-evalue 1e-10 -outfmt "6 qcovs qseqid sseqid pident length mismatch gapopen evalue" > putative_casbene_synthase_E_lathyris_blast.out
sort -k 4nr putative_casbene_synthase_E_lathyris_blast.out > sorted_putative_casbene_synthase_E_lathyris_blast.out
cd ..

cd CYP71D445
blastp -query ../Euphorbia_peplus.aa -subject E_lathyris_CYP71D445.fasta -qcov_hsp_perc 80 \
-evalue 1e-10 -outfmt "6 qcovs qseqid sseqid pident length mismatch gapopen evalue" > putative_CYP71D445_blast.out
sort -k 4nr putative_CYP71D445_blast.out > sorted_putative_CYP71D445_blast.out
cd ..

cd CYP726A27
blastp -query ../Euphorbia_peplus.aa -subject E_lathyris_CYP726A27.fasta -qcov_hsp_perc 80 \
-evalue 1e-10 -outfmt "6 qcovs qseqid sseqid pident length mismatch gapopen evalue" > putative_CYP726A27_blast.out
sort -k 4nr putative_CYP726A27_blast.out > sorted_putative_CYP726A27_blast.out
cd ..

cd FPP_synthase
blastp -query ../Euphorbia_peplus.aa -subject E_pekinensis_FPP_synthase.fasta -qcov_hsp_perc 80 \
-evalue 1e-10 -outfmt "6 qcovs qseqid sseqid pident length mismatch gapopen evalue" > putative_FPP_synthase_blast.out
sort -k 4nr putative_FPP_synthase_blast.out > sorted_putative_FPP_synthase_blast.out
cd ..

cd GGDP_synthase
tblastn -query ../Euphorbia_peplus.aa -subject E_peplus_GGDP_synthase_RNA_1.fasta -qcov_hsp_perc 80 \
-evalue 1e-10 -outfmt "6 qcovs qseqid sseqid pident length mismatch gapopen evalue" > putative_GGDP_synthase_1_blast.out
sort -k 4nr putative_GGDP_synthase_1_blast.out > sorted_putative_GGDP_synthase_1_blast.out
tblastn -query ../Euphorbia_peplus.aa -subject E_peplus_GGDP_synthase_RNA_2.fasta -qcov_hsp_perc 80 \
-evalue 1e-10 -outfmt "6 qcovs qseqid sseqid pident length mismatch gapopen evalue" > putative_GGDP_synthase_2_blast.out
sort -k 4nr putative_GGDP_synthase_2_blast.out > sorted_putative_GGDP_synthase_2_blast.out
cd ..

date
