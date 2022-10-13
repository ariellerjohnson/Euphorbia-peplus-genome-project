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

#loading in peplus data
Ep_chrom_lengths <- read.csv("E_peplus.lengths.karyoploter.txt", sep="\t", header=F)
Ep_repeats <- read.csv("E_peplus.TEs.500kwindows.cov.txt", sep="\t", header=F)
Ep_repeat_length <- read.csv("Ep_TE_lengths.mean.txt", sep="\t", header=F)
Ep_repeats_for_density <- gffToGRanges("E_peplus.KaryoPlotR.gff", zero.based = FALSE, ensembl = FALSE)
Ep_genes <- read.csv("E_peplus.genes.500kwindows.cov.txt", sep="\t", header=F)
Ep_centro <- read.csv("E_peplus.centro.500kwindows.cov.txt", sep="\t", header=F)

#loading in lathyris genome 
El_chrom_lengths <- read.csv("E_lathyris.lengths.karyoploter.txt", sep="\t", header=F)

#loading in TE data for peplus and lathyris 
Ep_Ty3 <- read.csv("E_peplus.Ty3.500kwindows.cov.txt", sep="\t", header=F)
Ep_Copia <- read.csv("E_peplus.Copia.500kwindows.cov.txt", sep="\t", header=F)
El_Ty3 <- read.csv("E_lathyris.Ty3.500kwindows.cov.txt", sep="\t", header=F)
El_Copia <- read.csv("E_lathyris.Copia.500kwindows.cov.txt", sep="\t", header=F)

#Adding E peplus genome with full chrom length
Ep.genome <- toGRanges(Ep_chrom_lengths)
seqlevels(Ep.genome)<- sub('Chromosome_','chr',seqlevels(Ep.genome))

#fix peplus repeat density track
seqlevels(Ep_repeats_for_density)<- sub('Chromosome_','chr',seqlevels(Ep_repeats_for_density))

#add peplus repeat coverage track
Ep_repeats[,1] <- sub('Chromosome_','chr',Ep_repeats[,1])
colnames(Ep_repeats)[7] <- "value"
#plot(Ep_repeats$value[which(Ep_repeats$V1=="chr1")]~Ep_repeats$V2[which(Ep_repeats$V1=="chr1")], type="l")
Ep_repeat_granges <- toGRanges(Ep_repeats)

#add peplus repeat length track 
Ep_repeat_length[,1] <- sub('Chromosome_','chr',Ep_repeat_length[,1])
Ep_repeat_length$V5 <- Ep_repeat_length$V4/max(Ep_repeat_length$V4)

#fix peplus gene coverage track
Ep_genes[,1] <- sub('Chromosome_','chr',Ep_genes[,1])
colnames(Ep_genes)[7] <- "value"

#fix peplus centromere coverage track
Ep_centro[,1] <- sub('Chromosome_','chr',Ep_genes[,1])
colnames(Ep_centro)[7] <- "value"

#fix Ty3 and Copia coverage tracks
Ep_Ty3[,1] <- sub('Chromosome_','chr',Ep_Ty3[,1])
colnames(Ep_Ty3)[7] <- "value"
Ep_Copia[,1] <- sub('Chromosome_','chr',Ep_Copia[,1])
colnames(Ep_Copia)[7] <- "value"
El_Ty3[,1] <- sub('LG','chr',El_Ty3[,1])
colnames(El_Ty3)[7] <- "value"
El_Copia[,1] <- sub('LG','chr',El_Copia[,1])
colnames(El_Copia)[7] <- "value"


#plotting peplus TE lengths to show why coverage and density look so different
kp <- plotKaryotype(genome=Ep.genome, plot.type = 4, cex=0.8)
kpAddLabels(kp, "TE length", r0=0.7, srt=90, label.margin=0.035,font=1.5, cex=0.6)
kpAddLabels(kp, "TE coverage", r0=0.15, srt=90, label.margin=0.035,font=1.5, cex=0.6)
kpAddLabels(kp, "TE density", r0=-0.6, srt=90, label.margin=0.035,font=1.5, cex=0.6)
kpAddBaseNumbers(kp, tick.dist = 10000000, tick.len = 2, tick.col="black", cex=0.6)
kp <- kpLines(kp, data=Ep_repeat_granges, data.panel = 1, col=("gray0"),y=Ep_repeat_length$V5,r0=0.69,r1=0.93)
kp <- kpLines(kp, data=Ep_repeat_granges, data.panel = 1, col=("gray0"),y=Ep_repeats$value,r0=0.36,r1=0.60)
kpPlotDensity(kp, data=Ep_repeats_for_density, col=alpha("black", 0.05), r0=0.0,r1=0.27, window.size=500000)

#plotting peplus genes, centro, and repeats coverage
png("peplus_KaryoPloteR_for_contact_map.png", width = 7, height = 3.2, units = 'in', res = 300)
kp <- plotKaryotype(genome=Ep.genome, plot.type = 4, cex=0.8)
kpAddLabels(kp, "genes", r0=0.7, srt=90, label.margin=0.035,font=1.5, cex=0.6)
kpAddLabels(kp, "centromere", r0=0.15, srt=90, label.margin=0.035,font=1.5, cex=0.6)
kpAddLabels(kp, "all repeats", r0=-0.55, srt=90, label.margin=0.035,font=1.5, cex=0.6)
kpAddBaseNumbers(kp, tick.dist = 10000000, tick.len = 2, tick.col="black", cex=0.6)
kp <- kpArea(kp, data=Ep_repeat_granges, data.panel = 1, col=("thistle2"),y=Ep_genes$value,r0=0.69,r1=0.93, ymax=0.6)
kpAxis(kp, ymax=0.6, r0=0.69,r1=0.93 ,cex=0.6,side=2)
kp <- kpArea(kp, data=Ep_repeat_granges, data.panel = 1, col=("lemonchiffon"),y=Ep_centro$value,r0=0.36,r1=0.60, ymax=0.9)
kpAxis(kp, ymax=0.9,r0=0.36,r1=0.60 ,cex=0.6,side=2)
kp <- kpArea(kp, data=Ep_repeat_granges, data.panel = 1, col=("snow2"),y=Ep_repeats$value,r0=0.0,r1=0.27, ymax=1)
kpAxis(kp, ymax=1, r0=0.0,r1=0.27,cex=0.6,side=2)
dev.off()

#peplus coverage of Copia and Ty3 across the chromosomes
png("peplus_KaryoPloteR_Copia_Ty3.png", width = 3, height = 3.2, units = 'in', res = 300)
kp <- plotKaryotype(genome = Ep.genome, plot.type = 1, cex=0.7)
kpAddBaseNumbers(kp, tick.dist = 10000000)
kp <- kpArea(kp, data=Ep_repeat_granges, col=("gold"),y=Ep_Copia$value,r0=0, r1=0.8, lty=0, ymax=0.7)
kp <- kpArea(kp, data=Ep_repeat_granges, col=("blue3"),y=Ep_Ty3$value,r0=0, r1=0.8, lty=0, ymax=0.7)
kpAxis(kp,cex=0.6,side=2, r0=0, r1=0.8, ymax=0.7, numticks=1)
dev.off()

#Adding E lathyris genome with full chrom length
El.genome <- toGRanges(El_chrom_lengths)
seqlevels(El.genome)<- sub('LG','chr',seqlevels(El.genome))
El_repeat_granges <- toGRanges(El_Copia)

#lathyris coverage of Copia and Ty3 across the chromosomes
png("lathyris_KaryoPloteR_Copia_Ty3.png", width = 6, height = 3.2, units = 'in', res = 300)
kp <- plotKaryotype(genome = El.genome, plot.type = 1, cex=0.7)
kpAddBaseNumbers(kp, tick.dist = 10000000)
kp <- kpArea(kp, data=El_repeat_granges, col=("gold"),y=El_Copia$value,r0=0, r1=0.8, lty=0, ymax=0.7)
kp <- kpArea(kp, data=El_repeat_granges, col=("blue3"),y=El_Ty3$value,r0=0, r1=0.8, lty=0, ymax=0.7)
kpAxis(kp,ymax=0.7,cex=0.6,side=2, r0=0, r1=0.8, numticks = 1)
dev.off()


