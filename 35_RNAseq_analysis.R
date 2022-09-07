#Found this vignette useful: https://github.com/CBC-UCONN/RNA-seq-with-reference-genome-and-annotation

# Load the libraries 
library("DESeq2")
library("apeglm")
library("pheatmap")
library("tidyverse")
library(ggpubr)

path_to_reads="~/Box Sync/Box Storage/Labwork/Genome paper/23_htseq_count"
setwd(path_to_reads)

list.files(path_to_reads)

sampleFiles <- c("10_ATCACG_Aligned.out.sorted.counts", 
                 "11_CGATGT_Aligned.out.sorted.counts", 
                 "12_TTAGGC_Aligned.out.sorted.counts",
                 "16_TGACCA_Aligned.out.sorted.counts",
                 "17_CAGATC_Aligned.out.sorted.counts",
                 "18_CTTGTA_Aligned.out.sorted.counts",
                 "13_ACAGTG_Aligned.out.sorted.counts",
                 "14_GCCAAT_Aligned.out.sorted.counts",
                 "15_GGCTAC_Aligned.out.sorted.counts",
                 "4_TGACCA_Aligned.out.sorted.counts",
                 "5_ACAGTG_Aligned.out.sorted.counts",
                 "6_GCCAAT_Aligned.out.sorted.counts",
                 "1_ATCACG_Aligned.out.sorted.counts",
                 "2_CGATGT_Aligned.out.sorted.counts",
                 "3_TTAGGC_Aligned.out.sorted.counts",
                 "7_GGCTAC_Aligned.out.sorted.counts",
                 "8_CTTGTA_Aligned.out.sorted.counts",
                 "9_TAGCTT_Aligned.out.sorted.counts")

# vector of sample names prior to removing bad root outlier sample
sampleNames <- c("young_leaf_1","young_leaf_2","young_leaf_3",
                  "mature_leaf_1", "mature_leaf_2", "mature_leaf_3",
                  "stem_1", "stem_2", "stem_3",
                  "inflorescence_1", "inflorescence_2", "inflorescence_3",
                  "fruit_1", "fruit_2", "fruit_3",
                  "root_1", "root_2", "root_3")


# vector of conditions prior to removing bad root outlier sample
sampleCondition <- c("young_leaf", "young_leaf", "young_leaf",
                     "mature_leaf", "mature_leaf", "mature_leaf",
                     "stem", "stem", "stem",
                     "inflorescence", "inflorescence", "inflorescence",
                     "fruit", "fruit","fruit",
                     "root", "root", "root")

# create data frame of data with bad root sample included
sampleTable <- data.frame(
  sampleName = sampleNames,
  fileName = sampleFiles,
  condition = sampleCondition
)

# create the DESeq data object
ddsHTSeq <- DESeqDataSetFromHTSeqCount(
  sampleTable = sampleTable, 
  directory = path_to_reads, 
  design = ~ condition
)


##PCA PLOT INCLUDING BAD ROOT SAMPLE

# get genes with summed counts greater than 20
sumcounts <- rowSums(counts(ddsHTSeq))
keep <- sumcounts > 20

# keep only the genes for which the vector "keep" is TRUE
ddsHTSeq <- ddsHTSeq[keep,]

dds <- DESeq(ddsHTSeq, fitType = "local")

# normalized, variance-stabilized transformed counts for visualization
vsd <- vst(dds, blind=TRUE)

plotPCA(vsd, intgroup="condition")

# alternatively, using ggplot

dat <- plotPCA(vsd, intgroup="condition",returnData=TRUE)

dat$group <- factor(dat$group,
                     levels = c("young_leaf", "mature_leaf", "stem", "inflorescence", "fruit", "root"))

#shape version
p <- ggplot(dat,aes(x=PC1,y=PC2,shape=group))
p <- p + geom_point()
p
ggsave("all_samples_RNAseq_PCA_shapes.png")

#color version
all_sample_PCA_color <- ggplot(dat,aes(x=PC1,y=PC2,col=group)) + geom_point() + theme_classic()
all_sample_PCA_color
ggsave("all_samples_RNAseq_PCA_colors.png")

#REMOVING FILE 8 (DEGRADED ROOT SAMPLE)
sampleFiles2 <- sampleFiles[! sampleFiles %in% c("8_CTTGTA_Aligned.out.sorted.counts")]

# new vector of sample names
sampleNames2 <- c("young_leaf_1","young_leaf_2","young_leaf_3", 
                 "mature_leaf_1", "mature_leaf_2", "mature_leaf_3",
                 "stem_1", "stem_2", "stem_3", 
                 "inflorescence_1", "inflorescence_2", "inflorescence_3", 
                 "fruit_1", "fruit_2", "fruit_3",
                 "root_1", "root_3")

# new vector of conditions
sampleCondition2 <- c("young_leaf", "young_leaf", "young_leaf",
                     "mature_leaf", "mature_leaf", "mature_leaf",
                     "stem", "stem", "stem",
                     "inflorescence", "inflorescence", "inflorescence", 
                     "fruit", "fruit","fruit",
                     "root", "root")

# now create a data frame from these three vectors. 
sampleTable2 <- data.frame(
  sampleName = sampleNames2,
  fileName = sampleFiles2,
  condition = sampleCondition2
)

# create the DESeq data object
ddsHTSeq2 <- DESeqDataSetFromHTSeqCount(
  sampleTable = sampleTable2, 
  directory = path_to_reads, 
  design = ~ condition
)


##CHANGE ORDER OF TREATMENTS

# change order of treatments:
treatments <- c("young_leaf", "mature_leaf", "stem", "root", "inflorescence", "fruit")

# reset the factor levels:
ddsHTSeq2$condition <- factor(ddsHTSeq2$condition, levels = treatments)

# verify the order
ddsHTSeq2$condition


##PCA PLOT WITH BAD ROOT SAMPLE REMOVED

# get genes with summed counts greater than 20
sumcounts <- rowSums(counts(ddsHTSeq2))
keep <- sumcounts > 20

# keep only the genes for which the vector "keep" is TRUE
ddsHTSeq2_keep <- ddsHTSeq2[keep,]

dds2 <- DESeq(ddsHTSeq2_keep, fitType = "local")

# normalized, variance-stabilized transformed counts for visualization
vsd2 <- vst(dds2, blind=TRUE)

plotPCA(vsd2, intgroup="condition")

# alternatively, using ggplot

dat2 <- plotPCA(vsd2, intgroup="condition",returnData=TRUE)

dat2$group <- factor(dat2$group,
                     levels = c("young_leaf", "mature_leaf", "stem", "inflorescence", "fruit", "root"))


#shape version
p <- ggplot(dat2,aes(x=PC1,y=PC2,shape=group))
p <- p + geom_point()
p
ggsave("final_RNAseq_PCA_shapes.png")

#color version
final_PCA_color <- ggplot(dat2,aes(x=PC1,y=PC2,col=group)) + geom_point()+theme_classic()
final_PCA_color
ggsave("final_RNAseq_PCA_colors.png")

figure <- ggarrange(all_sample_PCA_color, final_PCA_color,
                    labels = c("A", "B"),
                    ncol = 2, nrow = 1, common.legend = TRUE, legend = "bottom")
ggsave("RNAseq_PCA_figure.png", plot=figure)


###CALCULATE TPM 

sequence_counts=as.data.frame(counts(ddsHTSeq2))
sequence_counts[] <- lapply(sequence_counts, as.numeric)
head(sequence_counts)
sequence_counts$gene=row.names(sequence_counts)
head(sequence_counts)

gene_lengths=read.table("Euphorbia_peplus.genelength",sep="\t", header=T)
head(gene_lengths)
gene_lengths=select(gene_lengths,-(median:merged))
names(gene_lengths)[names(gene_lengths) == 'mean'] <- 'length'
gene_lengths$length <- as.numeric(gene_lengths$length)
head(gene_lengths)

counts_and_lengths=merge(gene_lengths, sequence_counts, by="gene")
head(counts_and_lengths)

rpk=counts_and_lengths[,3:ncol(counts_and_lengths)] / (counts_and_lengths[,2]/1000)
head(rpk)
rpk_matrix= as.matrix(rpk)
rpk=cbind(counts_and_lengths$gene, rpk)
names(rpk)[names(rpk) == 'counts_and_lengths$gene'] <- 'gene'
head(rpk)

colSums(rpk_matrix)
tpm=prop.table(rpk_matrix, margin = 2)*1000000
head(tpm)
tpm=as.data.frame(tpm)
tpm=cbind(counts_and_lengths$gene, tpm)
names(tpm)[names(tpm) == 'counts_and_lengths$gene'] <- 'gene'
head(tpm)
write.table(tpm, file='Euphorbia_peplus_tpm.tsv', quote=FALSE, sep='\t',row.names=F)
#quickly manually edited to just chromosomes in excel
#put both tpm tsv files in 22_RNAseq_PCAs_and_TPM on server

##HEATMAPS

# regularized log transformation of counts
# a good explanation: https://hbctraining.github.io/DGE_workshop_salmon_online/lessons/03_DGE_QC_analysis.html
rld <- rlog(dds2, blind=FALSE)

library("stringr") 
rownames(rld)=str_replace_all(rownames(rld), "sc0001", "chr1") 
rownames(rld)=str_replace_all(rownames(rld), "sc0002", "chr2") 
rownames(rld)=str_replace_all(rownames(rld), "sc0003", "chr3") 
rownames(rld)=str_replace_all(rownames(rld), "sc0004", "chr4") 
rownames(rld)=str_replace_all(rownames(rld), "sc0005", "chr5") 
rownames(rld)=str_replace_all(rownames(rld), "sc0006", "chr6") 
rownames(rld)=str_replace_all(rownames(rld), "sc0007", "chr7") 
rownames(rld)=str_replace_all(rownames(rld), "sc0008", "chr8") 

#Get the gene functions csv (used for chromoMap)
gene_functions=read.csv("../35_1_gene_functions.csv", header=F)
head(gene_functions)

my_ing_meb_genes=gene_functions[,1]
rownames(gene_functions)=my_ing_meb_genes
my_gene_functions=gene_functions[,3]


#length(which(rownames(rld) %in% my_ing_meb_genes==TRUE))

pheatmap(
  assay(rld)[rownames(rld) %in% my_ing_meb_genes,], 
  cluster_rows=FALSE, 
  show_rownames=TRUE,
  cluster_cols=FALSE
)

ing_meb_heatmap_table <- assay(rld)[rownames(rld) %in% my_ing_meb_genes,]
ing_meb_heatmap_table_averaged <- data.frame(row.names=rownames(ing_meb_heatmap_table))
ing_meb_heatmap_table_averaged$"young leaf"<- rowMeans(ing_meb_heatmap_table[, grep("young_leaf", colnames(ing_meb_heatmap_table))])
ing_meb_heatmap_table_averaged$"mature leaf"<- rowMeans(ing_meb_heatmap_table[, grep("mature_leaf", colnames(ing_meb_heatmap_table))])
ing_meb_heatmap_table_averaged$"stem"<- rowMeans(ing_meb_heatmap_table[, grep("stem", colnames(ing_meb_heatmap_table))])
ing_meb_heatmap_table_averaged$"inflorescence"<- rowMeans(ing_meb_heatmap_table[, grep("inflorescence", colnames(ing_meb_heatmap_table))])
ing_meb_heatmap_table_averaged$"fruit"<- rowMeans(ing_meb_heatmap_table[, grep("fruit", colnames(ing_meb_heatmap_table))])
ing_meb_heatmap_table_averaged$"root"<- rowMeans(ing_meb_heatmap_table[, grep("root", colnames(ing_meb_heatmap_table))])
averaged_heatmap_table_with_names <- merge(ing_meb_heatmap_table_averaged, gene_functions, by=0)
rownames(averaged_heatmap_table_with_names) <- paste(averaged_heatmap_table_with_names$V1, "--", averaged_heatmap_table_with_names$V3)
averaged_heatmap_table_with_names <- averaged_heatmap_table_with_names %>% select(-(V1:V4)) %>% select(-(Row.names)) 
averaged_heatmap_table_short_names <- merge(ing_meb_heatmap_table_averaged, gene_functions, by=0)
rownames(averaged_heatmap_table_short_names) <- make.unique(averaged_heatmap_table_short_names$V3)
averaged_heatmap_table_short_names <- averaged_heatmap_table_short_names %>% select(-(V1:V4)) %>% select(-(Row.names)) 


pheatmap(
  ing_meb_heatmap_table_averaged, 
  cluster_rows=FALSE, 
  show_rownames=TRUE,
  cluster_cols=TRUE,
  show_colnames = TRUE,
  #color = colorRampPalette(c("light yellow", "gold", "red"))(50), 
  scale = "row"
)


im_heatmap_transpose <- t(ing_meb_heatmap_table_averaged)

pheatmap(
  im_heatmap_transpose, 
  cluster_rows=TRUE, 
  show_rownames=TRUE,
  cluster_cols=FALSE,
  show_colnames = TRUE,
  #color = colorRampPalette(c("light yellow", "gold", "red"))(50), 
  scale = "column"
)

im_heatmap_with_names_transpose <- t(averaged_heatmap_table_with_names)
im_heatmap_short_names_transpose <- t(averaged_heatmap_table_short_names)


pheatmap(
  im_heatmap_with_names_transpose, 
  cluster_rows=TRUE, 
  show_rownames=TRUE,
  cluster_cols=FALSE,
  show_colnames = TRUE,
  #color = colorRampPalette(c("light yellow", "gold", "red"))(50), 
  scale = "column", 
  fontsize = 20
)


png("../heatmap for figure.png", width = 2000, height = 700, units = "px")
pheatmap(
  im_heatmap_short_names_transpose, 
  cluster_rows=TRUE, 
  show_rownames=TRUE,
  cluster_cols=FALSE,
  show_colnames = TRUE,
  #color = colorRampPalette(c("light yellow", "gold", "red"))(50), 
  scale = "column", 
  fontsize = 20
)
dev.off()

pheatmap(
  im_heatmap_transpose, 
  cluster_rows=TRUE, 
  show_rownames=TRUE,
  cluster_cols=TRUE,
  show_colnames = TRUE,
  #color = colorRampPalette(c("light yellow", "gold", "red"))(50), 
  scale = "column"
)
