#!/usr/bin/env Rscript

### PoolSeq: Pool sequencing of spores coming from 4- vs 2-spore asci of Y vs Wa87p
#############################################################################
# To confirm that Spok2 is responsible for killing in crosses between Psk-1
# and Psk-5 a pooled sequencing approach was employed. A cross was conducted
# between Wa87 and Y, and spores from 2-spored (spore killing) and 4-spored
# asci (heterozygous for killers) were collected. We obtained 21 progeny from
# 2-spored asci and 63 progeny from 4-spored asci. Progeny from each ascus
# type were grown on the same plate and extracted for sequencing with Illumina
# Hi-seq X.
# =======================================
# Sandra Lorena Ament Velasquez
# Johannesson Lab, Evolutionary Biology Center, Uppsala University, Sweden
# 2018-08-28
# Version 0
# =======================================

library(vcfR) # browseVignettes("vcfR")
library(poppr)
library(ggplot2)
library(dplyr) # For filter
library(reshape2) # To fix dataframe with melt

# ============================
# Reading the data
# ============================
# Getting all sites for all samples
vcffile <- snakemake@input[[1]]

# READ THE VCF
vcf <- read.vcfR(vcffile, verbose = FALSE)

## Remove the mitochondria
indexchr <- which(vcf@fix[, "CHROM"] != "PaMt_NC_001329.3")  # Get the range of positions excluding that chr
vcf@fix <- vcf@fix[vcf@fix[, "CHROM"] != "PaMt_NC_001329.3", ]  # Get only the relevant chromosomes
vcf@gt <- vcf@gt[indexchr, ] # Get only the relevant genotypes using the range of the chromosome

## Get the depth of coverage 
dp <- extract.gt(vcf, element='DP', as.numeric = TRUE) %>% data.frame()
dp <- tibble::rownames_to_column(dp, "ChrPos") %>% data.frame()

ad <- extract.gt(vcf, element='AD', as.numeric = TRUE) %>% data.frame() # Trick here is that as.numeric forces to take only the coverage of reference
ad <- tibble::rownames_to_column(ad, "ChrPos") %>% data.frame()

# Get genotypes
gt <- extract.gt(vcf, element='GT', return.alleles = TRUE) %>% data.frame()
gt <- tibble::rownames_to_column(gt, "ChrPos") %>% data.frame()


# Reshape
dp <- melt(dp, id = c("ChrPos"))
names(dp) <- c("ChrPos", "Sample", "DP")
ad <- melt(ad, id = c("ChrPos"))
names(ad) <- c("ChrPos", "Sample", "RefAll")
gt <- melt(gt, id = c("ChrPos"))
names(gt) <- c("ChrPos", "Sample", "allele")


## Make the coordinates easier to handle
positions <- strsplit(dp$ChrPos, "chromosome_._") %>% unlist
positions <- positions[positions != ""] %>% as.numeric

# Make a column for each chromosome
chr <- strsplit(dp$ChrPos, "_[0-9]+$", perl = TRUE) %>% unlist

# Add to dataframe
genotypes <- cbind(chr, positions, dp, Ref = ad$RefAll, Alt = abs(dp$DP - ad$RefAll), allele = gt$allele)

# Calculate the Minor allele frequency 
MAFs <- genotypes %>% rowwise() %>% mutate(maf = min(Ref, Alt)/DP, relRef = Ref/DP, relAlt = Alt/DP) %>% data.frame() # rowwise makes the trick to do it per row

### Now, I want to filter out everything that is not "haploid-looking" in the parentals
# Get indexes of bad SNPs per sample
indexpersamp_2sp <- which(MAFs %>% filter(Sample == "PaWa87m") %>% .$maf >= 0.2 )
indexpersamp_4sp <- which(MAFs %>% filter(Sample == "PaYp") %>% .$maf >= 0.2 )
# Merge the two lists of SNPs
badsnps <- c(indexpersamp_2sp, indexpersamp_4sp) %>% sort %>% unique()

# Number of SNPs per sample
nrowPerSample <- MAFs %>% filter(Sample == "Pool2spores") %>% nrow()
albadsnps <- c(badsnps, badsnps + nrowPerSample, badsnps + nrowPerSample*2, badsnps + nrowPerSample*3 ) # Expand the bad snps to all samples
# This filters out all sites that are too polymorphic in the haploid parents, under the assumption that the SNPs come from repeats, probably

## Data frame including only the poolsequencing
pools <- MAFs %>% filter(Sample == "Pool2spores" | Sample == "Pool4spores")


#### --- Points ----
centromeres <- data.frame(chr = c("chromosome_1", "chromosome_2", "chromosome_3", "chromosome_4", "chromosome_5", "chromosome_6", "chromosome_7"),
                          positions = c(4479205, 236468.5, 675115.5, 1236388.5, 2062808.5, 238150.0, 3562141.5), Locus = "Centromere", Sample = c(rep("Pool2spores", 7), rep("Pool4spores", 7)), DP = 0, Ref = 0, Alt = 0, maf = 0)

## Map some interesting genes
markers <- data.frame(chr = c("chromosome_5", "chromosome_1", "chromosome_3", "chromosome_3", "chromosome_4", "chromosome_2", "chromosome_1", "chromosome_2", "chromosome_3"),
                      positions = c(1094363, 7345878.5, 208938.5, 495787.5, 691354.5, 3056777, 3957511.5, 2683360.5, 355881), 
                      Locus = c("Spok2", "MAT", "het-s", "het-c", "het-e", "het-d", "het-Z", "het-r", "Spok3/4"), 
                      Sample = c(rep("Pool2spores", 9), rep("Pool4spores", 9)), 
                      DP = 0, Ref = 0, Alt = 0, maf = 0)

# Put them together
allmarkers <- rbind(centromeres, markers)

### ---- Keep track of the parental ----
filteredMAFs <- MAFs[-albadsnps,] # %>% filter(DP >= 50 & DP <=150 & maf >= 0.01)
alleleWa87 <- filteredMAFs %>% filter(Sample == "PaWa87m") %>% .$allele
alleleYp <- filteredMAFs %>% filter(Sample == "PaYp") %>% .$allele

spores2 <- filteredMAFs %>% filter(Sample == "Pool2spores")
spores4 <- filteredMAFs %>% filter(Sample == "Pool4spores")

# Extract individual alleles from the dikaryons
getalleledik <- function(dikaryondf, alleleWa87, alleleYp){
  allele1 <- strsplit(dikaryondf$allele %>% as.character, "/.") %>% unlist
  allele2 <- strsplit(dikaryondf$allele %>% as.character, "./") %>% unlist
  allele2 <- allele2[allele2 != ""]
  
  all1 <- rep("?", dikaryondf %>% nrow())
  all1[alleleWa87 == allele1] <- "Wa87"
  all1[alleleYp == allele1] <- "Y"
  
  all2 <- rep("?", dikaryondf %>% nrow())
  all2[alleleWa87 == allele2] <- "Wa87"
  all2[alleleYp == allele2] <- "Y"
  
  return(cbind(all1, all2))
}

spores2 <- cbind(spores2, getalleledik(spores2, alleleWa87, alleleYp)) 
spores4 <- cbind(spores4, getalleledik(spores4, alleleWa87, alleleYp))

# Remove the sites where one of the parents is not recognizible (and filter out bad quality sites too)
dikaryons <- rbind(spores2, spores4) %>% filter(all1 != "?" & all2 != "?") %>% filter(DP >= 50 & DP <=150 & maf >= 0.01)

### ---- Re-shape the data-frame ----
allelecol <- dikaryons[,c(1,2,4,5,12,13)] %>% melt(., id = c("chr", "positions", "Sample", "DP"))
allelefreq <- dikaryons[,c(1,2,4,5,10,11)] %>% melt(., id = c("chr", "positions", "Sample", "DP"))

dikaryonsfreq <- cbind(allelefreq, Parent = allelecol$value)
names(dikaryonsfreq) <- c("chr", "positions", "Sample", "DP", "Allele", "Freq", "Parent")

### ---- Finally plot ----

allpoints <- rbind ( cbind(centromeres, Freq = 0), cbind(markers, Freq = 0) )

ggplot(dikaryonsfreq, aes(x = positions, y = Freq, colour = Parent)) + 
  geom_point(alpha = 0.2) + 
  facet_grid(chr ~ Sample, scales = "free") +
  geom_point(data = allpoints, aes(x = positions, shape = Locus, fill = Locus), colour = "black", size = 2) + 
  scale_shape_manual(values=rep(c(19, 18, 13, 23, 25, 15, 9, 3, 4, 24),2)) +
  theme(panel.background = element_rect(fill = "gray98")) +
  xlab("Position (bp)") + ylab("Allele frequency")


ggsave(file = snakemake@output[[1]], width = 10, height = 11, scale = 2.5, units = "cm")
ggsave(file = paste0(snakemake@output[[1]],".png"), width = 10, height = 11, scale = 2.5, units = "cm")
