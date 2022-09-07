##GENESPACE for E. peplus genome
#Arielle Johnson 2/17/22

library("GENESPACE")

gpar <- init_genespace(
  genomeIDs = c("Euphorbia peplus", "Euphorbia lathyris", "Manihot esculenta", "Hevea brasiliensis", "Ricinus communis"),
  speciesIDs = c("EuphorbiaPeplus", "EuphorbiaLathyris", "ManihotEsculenta", "HeveaBrasiliensis", "RicinusCommunis"),
  versionIDs = c( "v", "v", "v", "v", "v"),
  ploidy = c(1,1,1,1,1),
  wd = "/workdir/arj66/GENESPACE/E_peplus_GENESPACE/Euphorbiaceae",
  orthofinderMethod = "fast",
  diamondMode = "fast",
  nCores = 20,
  gffString = "gff3",
  pepString = ".fa",
  path2mcscanx = "/programs/MCScanX",
  rawGenomeDir = "/workdir/arj66/GENESPACE/E_peplus_GENESPACE/genomes/")

parse_phytozome(gsParam = gpar, genomeIDs=c("Manihot esculenta"), overwrite = T)

#parse_ncbi(gsParam = gpar, genomeIDs = "Jatropha curcas", overwrite = T)

parse_annotations(
  gsParam = gpar, genomeIDs = "Hevea brasiliensis", troubleshoot = F, 
  gffEntryType = "CDS", gffIdColumn = "protein_id", headerStripText = ".*=|]", 
  overwrite = T)

parse_annotations(
  genomeIDs = c("Euphorbia lathyris", "Ricinus communis"),
  gsParam = gpar,
  gffEntryType = "CDS",
  headerEntryIndex = 1,
  overwrite = T,
  troubleshoot = T,
  headerSep=" ",
  gffIdColumn = "Parent")

parse_annotations(
  genomeIDs = "Euphorbia peplus",
  gffEntryType = "CDS",
  gsParam = gpar,
  headerEntryIndex = 1,
  overwrite = T,
  troubleshoot = T,
  headerSep=" ",
  gffIdColumn = "transcript_id")


gpar <- run_orthofinder(gsParam = gpar, overwrite = T)
gpar <-  synteny(gsParam=gpar, overwrite = T)

png("Peplus_GENESPACE.png", width = 6, height = 4, units = 'in', res = 300)
ripdat <- plot_riparian(gpar, 
                        #https://colorbrewer2.org/#type=sequential&scheme=YlGnBu&n=8, plus "#F0E442"
                        colByChrs = c("#F0E442", "#edf8b1", "#c7e9b4", "#7fcdbb", "#41b6c4", "#1d91c0", "#225ea8", "#0c2c84"),
                        onlyTheseChrs = c("scaffold_1", "scaffold_2", "scaffold_3", "scaffold_4", "scaffold_5", "scaffold_6", "scaffold_7", "scaffold_8"),
                        blackBg = FALSE)
dev.off()

png("Peplus_chr1_GENESPACE.png", width = 6, height = 4, units = 'in', res = 300)
ripdat <- plot_riparian(gpar, 
                        #https://colorbrewer2.org/#type=sequential&scheme=YlGnBu&n=8, plus "#F0E442"
                        colByChrs = c("#F0E442", "#edf8b1", "#c7e9b4", "#7fcdbb", "#41b6c4", "#1d91c0", "#225ea8", "#0c2c84"),
                        onlyTheseChrs = c("scaffold_1"),
                        blackBg = FALSE)
dev.off()

png("Peplus_chr2_GENESPACE.png", width = 6, height = 4, units = 'in', res = 300)
ripdat <- plot_riparian(gpar, 
                        #https://colorbrewer2.org/#type=sequential&scheme=YlGnBu&n=8, plus "#F0E442"
                        colByChrs = c("#edf8b1", "#c7e9b4", "#7fcdbb", "#41b6c4", "#1d91c0", "#225ea8", "#0c2c84"),
                        onlyTheseChrs = c("scaffold_2"),
                        blackBg = FALSE)
dev.off()

png("Peplus_chr3_GENESPACE.png", width = 6, height = 4, units = 'in', res = 300)
ripdat <- plot_riparian(gpar, 
                        #https://colorbrewer2.org/#type=sequential&scheme=YlGnBu&n=8, plus "#F0E442"
                        colByChrs = c("#c7e9b4", "#7fcdbb", "#41b6c4", "#1d91c0", "#225ea8", "#0c2c84"),
                        onlyTheseChrs = c("scaffold_3"),
                        blackBg = FALSE)
dev.off()

png("Peplus_chr4_GENESPACE.png", width = 6, height = 4, units = 'in', res = 300)
ripdat <- plot_riparian(gpar, 
                        #https://colorbrewer2.org/#type=sequential&scheme=YlGnBu&n=8, plus "#F0E442"
                        colByChrs = c("#7fcdbb", "#41b6c4", "#1d91c0", "#225ea8", "#0c2c84"),
                        onlyTheseChrs = c("scaffold_4"),
                        blackBg = FALSE)
dev.off()

png("Peplus_chr5_GENESPACE.png", width = 6, height = 4, units = 'in', res = 300)
ripdat <- plot_riparian(gpar, 
                        #https://colorbrewer2.org/#type=sequential&scheme=YlGnBu&n=8, plus "#F0E442"
                        colByChrs = c("#41b6c4", "#1d91c0", "#225ea8", "#0c2c84"),
                        onlyTheseChrs = c("scaffold_5"),
                        blackBg = FALSE)
dev.off()

png("Peplus_chr6_GENESPACE.png", width = 6, height = 4, units = 'in', res = 300)
ripdat <- plot_riparian(gpar, 
                        #https://colorbrewer2.org/#type=sequential&scheme=YlGnBu&n=8, plus "#F0E442"
                        colByChrs = c("#1d91c0", "#225ea8", "#0c2c84"),
                        onlyTheseChrs = c("scaffold_6"),
                        blackBg = FALSE)
dev.off()

png("Peplus_chr7_GENESPACE.png", width = 6, height = 4, units = 'in', res = 300)
ripdat <- plot_riparian(gpar, 
                        #https://colorbrewer2.org/#type=sequential&scheme=YlGnBu&n=8, plus "#F0E442"
                        colByChrs = c("#225ea8", "#0c2c84"),
                        onlyTheseChrs = c("scaffold_7"),
                        blackBg = FALSE)
dev.off()

png("Peplus_chr8_GENESPACE.png", width = 6, height = 4, units = 'in', res = 300)
ripdat <- plot_riparian(gpar, 
                        #https://colorbrewer2.org/#type=sequential&scheme=YlGnBu&n=8, plus "#F0E442"
                        colByChrs = c("#0c2c84"),
                        onlyTheseChrs = c("scaffold_8"),
                        blackBg = FALSE)
dev.off()