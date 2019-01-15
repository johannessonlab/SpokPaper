# The 5xS backcrosses of the *P. anserina*'s Psk
#############################################################################
# ===========================================================================
# Sandra Lorena Ament Velasquez
# Johannesson Lab, Evolutionary Biology Center, Uppsala University, Sweden
# 2018/11/27
# ++++++++++++++++++++++++++++++++++++++++++++++

# ============================
# Load the necessary libraries
# ============================
library(vcfR) # browseVignettes("vcfR")
library(poppr) # To manage the SNPs
library(ggplot2)
library(dplyr) # For filter
library(reshape2) # To fix dataframe
# ============================

# library(ggfortify) # https://cran.r-project.org/web/packages/ggfortify/vignettes/plot_pca.html
# http://www.molecularecologist.com/2011/03/should-i-use-fst-gst-or-d-2/
# http://ggplot2.tidyverse.org/reference/geom_tile.html
# http://t-redactyl.io/blog/2015/12/creating-plots-in-r-using-ggplot2-part-1-line-plots.html

# ============================
# Reading the data
# ============================
# Getting all sites for all samples
vcffile <- snakemake@input[[1]]

# ============================
# READ THE VCF
# ============================
vcf <- read.vcfR(vcffile, verbose = FALSE)

# Create the adegenet objects
genlight <- vcfR2genlight(vcf) # That's faster, omitting 289 loci with more than 2 alleles

# Subset the genlight object for only the relevant strains
backcrosssamples <- c("Pa170m", "Pa180p", "Pa130p", "Pa200p", "PaSp", "PaSm", "PaWa53p", "PaWa58m", "PaWa28m", "PaYp") # With parentals
snpclone <- genlight[genlight@ind.names %in% backcrosssamples] %>% as.snpclone() # Make snpclone object

# ============================
# Shape the snpclone into a dataframe
# ============================
# First rename the individuals to look like in the paper (including parentals)
snpclone@ind.names <- c("Psk2xS5p","Psk1xS5m","Psk7xS5p","Psk5xS5p","Sm","Sp", "PaWa28m", "PaWa53p", "PaWa58m", "PaYp")

### Make a nice data frame
genotypes_matrix <- t(as.matrix(snpclone))
genotypes <- data.frame(site = row.names(genotypes_matrix), 
                        pos = snpclone@position,
                        genotypes_matrix, 
                        chr = snpclone@chromosome)

genotypes <- genotypes %>% melt(id = c("site", "pos", "chr"))
names(genotypes)[4:5] <- c("Strain", "Allele")

# Rearrange labels to be in pretty order
# genotypes$Strain <- factor(genotypes$Strain, levels = rev(c("Sm", "Sp", "Psk1xS5m", "Psk2xS5p", "Psk5xS5p", "Psk7xS5p", "PaWa53p","PaWa28m", "PaYp", "PaWa58m", character(1))))
genotypes$Strain <- factor(genotypes$Strain, levels = rev(c("Sm", "Sp", "Psk1xS5m", "PaWa53p", "Psk2xS5p", "PaWa28m", "Psk5xS5p", "PaYp", "Psk7xS5p", "PaWa58m", character(1))))

# ============================
# Extra data frames with marker's information
# ============================

# Make a dataframe with the lengths of the Podan2 chromosomes
chrlines <- data.frame(pos =c(0), endchr = c(8813524, 5165621, 4137471, 3808395, 4734292, 4264132, 4087160) , chr = c("chromosome_1", "chromosome_2", "chromosome_3", "chromosome_4", "chromosome_5", "chromosome_6", "chromosome_7"), Strain = character(1))

# The genes
genes <- data.frame(site = c("Spok2", "Spok3", "Spok3", "Spok3", "Spok3", "Spok4", "Spok4", "Spok4",
                             "MAT", rep("Centromere", 7), 
                             "Spok3", "Spok3", "Spok3", "Spok3", "Spok4", "Spok4", "Spok4"), 
                    pos =c(1094363, 355881, 355881, 3328395.5, 900387.5, 355881, 355881, 900387.5,
                           7345878.5, 4479205, 236468.5, 675115.5, 1236388.5, 2062808.5, 238150, 3562141.5,
                           355881, 355881, 3328395.5, 900387.5, 355881, 355881, 900387.5), 
                    chr = c("chromosome_5", "chromosome_3", "chromosome_3", "chromosome_5", "chromosome_5", "chromosome_3", "chromosome_3", "chromosome_5",
                            "chromosome_1","chromosome_1", "chromosome_2", "chromosome_3", "chromosome_4", "chromosome_5", "chromosome_6", "chromosome_7",
                            "chromosome_3", "chromosome_3", "chromosome_5", "chromosome_5", "chromosome_3", "chromosome_3", "chromosome_5"), 
                    Strain = c(character(1), "Psk1xS5m", "Psk5xS5p", "Psk2xS5p", "Psk7xS5p", "Psk1xS5m", "Psk5xS5p", "Psk7xS5p",
                               rep(character(1), 8), 
                               "PaYp", "PaWa53p", "PaWa28m", "PaWa58m", "PaWa53p", "PaYp", "PaWa58m"))

# Het genes for these strains with alleles
hets <- data.frame(site = c("het-S", "het-S", "het-S", "het-S", "het-S", "het-S", "het-S", "het-S", "het-S", "het-S",
                            "het-c", "het-c", "het-c", "het-c", "het-c", "het-c", "het-c", "het-c", "het-c", "het-c"), 
                   pos =c(rep(208504, 10),
                          rep(495416, 10)), 
                   chr = c("chromosome_3"), 
                   Strain = rep(c("Sm", "Sp", "Psk1xS5m", "Psk2xS5p", "Psk5xS5p", "Psk7xS5p", 
                                  "PaWa53p","PaWa28m", "PaYp", "PaWa58m"), 2), 
                   Allele = c("S", "S", "s", "s", "s", "S", 
                              "s", "S", "s", "s", 
                              "c1", "c1", "c1", "c2", "c3", "c1",
                              "c1", "c6", "c3", "c3"))

# ============================
# Plot
# ============================
# Make base plot
allchrs <- ggplot(genotypes %>% filter(Allele == 1), aes(x = pos, y = Strain)) + xlab("Chromosome position (bp)") +
  stat_bin2d(binwidth = c(20000, 1)) + scale_fill_gradientn(colours = terrain.colors(5)) +
  theme(panel.background = element_rect(fill = "gray98") ) +
  facet_grid(chr ~ ., scales = "free") #+ coord_flip()

# Plot all chromosomes on top
allchrs <- allchrs + geom_segment(data = chrlines, 
                                  aes(x = pos, y = Strain, xend = endchr, yend = 1), 
                                  colour = "gray",
                                  size = 2) + 
  theme(axis.ticks=element_blank()) + # Remove ticks of strains
  scale_y_discrete(drop=FALSE) # Do not drop unused factors to restart your own
# I had to add scale_y_discrete(drop=FALSE) because adding layers into a plot with a new data frame defaults to drop the unused layers, fuse them, and have a new set of factors order alphabetically. See https://github.com/tidyverse/ggplot2/issues/577

# Add markers
allchrs <- allchrs + geom_point(data = genes, aes(x = pos, shape = site), size = 1.9) + # 
  scale_shape_manual(values=c(19, 17, 15, 3, 4)) + scale_y_discrete(drop=FALSE) # Do not drop unused factors to restart your own (but it conflicts with Y)

# Add het genes
modDarjeeling1 <- c("#FF0000", "#F2AD00", "#00F296", "#336EDE", "#C895E3", "#ED4FAF") # Nice colors
allchrs <- allchrs + geom_point(data = hets, aes(x = pos, shape = site, colour = Allele), size = 1.9) + 
  scale_shape_manual(values=c(19, 18, 15, 3, 4, 1,2)) + scale_colour_manual(values = modDarjeeling1)

#Save the final plot
ggsave(plot = allchrs, file = snakemake@output[[1]], width = 9, height = 10, scale = 2.5, units = "cm")

