#https://bioconductor.org/packages/release/bioc/html/karyoploteR.html

setwd("~/Box Sync/Box Storage/Labwork/Genome paper/KaryoPloteR")

# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("karyoploteR")

#BiocManager::install("genomation")

# browseVignettes("karyoploteR")
# browseVignettes("regioneR")

library(karyoploteR)
library(regioneR)
library(genomation)
library(scales)

#Adding E peplus genome with full chrom length
Ep.genome <- toGRanges(data.frame(chr=c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8"), start=c(1, 1, 1, 1, 1, 1, 1, 1), end=c(29603500, 29994419, 34601726, 37676872, 36807906, 38403033, 31008493, 29100068)))

# Add gene density
# 22_finalizing_genome_for_upload/E_peplus_chromosomal_assembly/Euphorbia_peplus.gff
Ep_genes <- gffToGRanges("Euphorbia_peplus.gff", filter = "gene", zero.based = FALSE, ensembl = FALSE)
#head(Ep_genes)
seqlevels(Ep_genes)<- sub('Chromosome_','chr',seqlevels(Ep_genes))
#head(Ep_genes)

# add centromeric repeat
#21_putative_centromeric_repeat/top_repeat_blast.out.gff
Ep_centro <- gffToGRanges("top_repeat_blast.out.gff", filter = NULL, zero.based = FALSE, ensembl = FALSE)
#head(Ep_centro)
seqlevels(Ep_centro)<- sub('HiC_scaffold_','chr',seqlevels(Ep_centro))
#head(Ep_centro)



###TEs
#31_comparative_KaryoPloteR_setup/*.KaryoPlotR.gff

#all repetitive sequences
Ep_repeats <- gffToGRanges("E_peplus.KaryoPlotR.gff", zero.based = FALSE, ensembl = FALSE)
seqlevels(Ep_repeats)<- sub('Chromosome_','chr',seqlevels(Ep_repeats))

#Ep_Copia
Ep_Copia <- gffToGRanges("E_peplus.KaryoPlotR.gff", filter = "LTR/Copia", zero.based = FALSE, ensembl = FALSE)
seqlevels(Ep_Copia)<- sub('Chromosome_','chr',seqlevels(Ep_Copia))

#Ep_Ty3 (formerly Gypsy)
Ep_Ty3 <- gffToGRanges("E_peplus.KaryoPlotR.gff", filter = "LTR/Gypsy", zero.based = FALSE, ensembl = FALSE)
seqlevels(Ep_Ty3)<- sub('Chromosome_','chr',seqlevels(Ep_Ty3))


#Ep_Harbinger
Ep_Harbinger <- gffToGRanges("E_peplus.KaryoPlotR.gff", filter = "DNA/PIF-Harbinger", zero.based = FALSE, ensembl = FALSE)
seqlevels(Ep_Harbinger)<- sub('Chromosome_','chr',seqlevels(Ep_Harbinger))


#Peplus map
png("peplus_KaryoPloteR_for_contact_map.png", width = 7, height = 3.2, units = 'in', res = 300)

kp <- plotKaryotype(genome=Ep.genome, plot.type = 4, cex=0.8)
kpAddLabels(kp, "genes", r0=0.7, srt=90, label.margin=0.035,font=1.5, cex=0.6)
kpAddLabels(kp, "centromere", r0=0.23, srt=90, label.margin=0.035,font=1.5, cex=0.6)
kpAddLabels(kp, "all repeats", r0=-0.5, srt=90, label.margin=0.035,font=1.5, cex=0.6)
kpAddBaseNumbers(kp, tick.dist = 10000000, tick.len = 2, tick.col="black", cex=0.6)

kpPlotDensity(kp, data=Ep_genes, col = alpha("purple", 0.05), r0=0.69,r1=0.93)
kpAxis(kp, ymax=kp$latest.plot$computed.values$max.density, r0=0.69,r1=0.93,cex=0.6,side=2)
kpPlotDensity(kp, data=Ep_centro, col = alpha("gold", 0.05), r0=0.36,r1=0.60)
kpAxis(kp, ymax=kp$latest.plot$computed.values$max.density, r0=0.36,r1=0.60,cex=0.6,side=2)
kpPlotDensity(kp, data=Ep_repeats, col=alpha("black", 0.05), r0=0.0,r1=0.27)
kpAxis(kp, ymax=kp$latest.plot$computed.values$max.density, r0=0.0,r1=0.27,cex=0.6,side=2)

dev.off()

#maps for peplus vs lathyris 

seqlevels(Ep_Ty3)<- sub('Chromosome_','chr',seqlevels(Ep_Ty3))

png("peplus_KaryoPloteR_Copia_Ty3.png", width = 4, height = 3.2, units = 'in', res = 300)
kp <- plotKaryotype(genome = Ep.genome, plot.type = 1, cex=0.7)
kpAddBaseNumbers(kp, tick.dist = 10000000)
kpPlotDensity(kp, data=Ep_Copia, col = alpha("gold", 1), r0=0, r1=0.8, lty=3)
kpPlotDensity(kp, data=Ep_Ty3, col = alpha("blue3", 0.2), r0=0, r1=0.8, lty=1)
dev.off()

#adding lathyris genome
El.genome <- toGRanges(data.frame(chr=c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10"), start=c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1), end=c(140418633, 138705850, 107888011, 96626626, 95724319, 88300888, 86325489, 79934230, 79259015, 73611708)))

#El_Copia
El_Copia <- gffToGRanges("E_lathyris.KaryoPlotR.gff", filter = "LTR/Copia", zero.based = FALSE, ensembl = FALSE)
seqlevels(El_Copia)<- sub('LG0','chr',seqlevels(El_Copia))
seqlevels(El_Copia)<- sub('LG','chr',seqlevels(El_Copia))

#El_Ty3 (formerly Gypsy)
El_Ty3 <- gffToGRanges("E_lathyris.KaryoPlotR.gff", filter = "LTR/Gypsy", zero.based = FALSE, ensembl = FALSE)
seqlevels(El_Ty3)<- sub('LG0','chr',seqlevels(El_Ty3))
seqlevels(El_Ty3)<- sub('LG','chr',seqlevels(El_Ty3))

png("lathyris_KaryoPloteR_Copia_Ty3.png", width = 4, height = 3.2, units = 'in', res = 300)
kp2 <- plotKaryotype(genome = El.genome, plot.type = 1, cex=0.7)
kpAddBaseNumbers(kp2, tick.dist = 10000000)
kpPlotDensity(kp2, data=El_Copia, col = alpha("gold", 1), r0=0, r1=0.8, lty=3)
kpPlotDensity(kp2, data=El_Ty3, col = alpha("blue3", 0.2), r0=0, r1=0.8, lty=1)
dev.off()
