setwd("~/Box Sync/Box Storage/Labwork/Genome paper/UpSetR")

#install.packages("UpSetR")
library(UpSetR)

#Input: from results folder: /Orthogroups/Orthogroups.GeneCount.tsv

#import orthogroup gene counts
gene_count <- read.table("Orthogroups.GeneCount.tsv", header = T, sep = "\t")
#head(gene_count)
gene_count <- gene_count[,2:6]
#head(gene_count)
colnames(gene_count) = gsub("_", ". ", colnames(gene_count))
#head(gene_count)

#make binary
gene_count[gene_count > 0] <- 1
#head(gene_count)

#?upset

#make upset plot with numbers on bars
upset(gene_count,  sets = c("H. brasiliensis", "M. esculenta", "R. communis", "E. lathyris", "E. peplus"), keep.order = T, sets.bar.color = "dark blue",
      order.by = "freq", empty.intersections = "on", show.numbers = "yes", mainbar.y.label = "Number of orthologous groups", sets.x.label = "Groups with >0 members from sp.")

#export image without numbers on bars
png(file="Orthogroups_UpSet_plot.png", width=7, height=5, units="in", res=200)
upset(gene_count,  sets = c("H. brasiliensis", "M. esculenta", "R. communis", "E. lathyris", "E. peplus"), keep.order = T, sets.bar.color = "dark blue",
      order.by = "freq", empty.intersections = "on", show.numbers = "no", mainbar.y.label = "Number of orthologous groups", sets.x.label = "Groups with >0 members from sp.")
dev.off()