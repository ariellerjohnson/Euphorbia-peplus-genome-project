# *Euphorbia peplus* genome assembly v1

## Brief explanation

This is the code for Johnson et al 2022: "Chromosome-level genome assembly of *Euphorbia peplus*, a model system for plant latex, reveals that suppression of Ty3 transposons contributed to its small genome size".

Code was run in a working directory on a Cornell BioHPC server.  All .sh code was run through SLURM with the command "sbatch [scriptname].sh" unless otherwise noted.  All .R code was run in RStudio.  Before starting to run code, I made an empty directory called "joblogs" as well as a directory called "raw_data" that contained six subdirectories called "assemblies_for_Kimura", "genes_of_interest", "GENESPACE_gff_files", "HiC_data", "HiFi_CCS_reads", "OrthoFinder_protein_files", and "RNAseq_data", all with the relevant data inside.  The raw data generated for this project, and the finished assembly, can be found at PRJNA837952 “Euphorbia peplus Genome sequencing and assembly”.

## Contact information 
I am a PhD student at Cornell, feel free to contact me through Github if you have any questions about this code.  The co-corresponding authors on the paper are Gaurav Moghe (gdm67@cornell.edu) and Margaret Frank (mhf47@cornell.edu), please contact them with any other questions about the project.  
