#!/usr/bin/env Rscript

# Calculating Pairwise Nucleotide Differences of Spoks
# =======================================
# Sandra Lorena Ament Velasquez
# Johannesson Lab, Evolutionary Biology Center, Uppsala University, Sweden
# 2018-10-31
# =======================================

library(ggplot2)
library(dplyr) # For filter
library(ape)

# Read alignment
spoks <- read.dna("SpoksPaper_CDS_al.fas", format = "fasta", skip = 0, nlines = 0, comment.char = "#", as.character = FALSE, as.matrix = NULL)

# A function to calculate pairwise nucleotide diversity by moving windows of fixed size from an aligment
# See discussion in http://seqanswers.com/forums/showthread.php?t=33639
winpiseq_step <- function(fasta, windowsize = 100, winstep = 20){
  alglen <- dim(fasta)[2] # The length of the alignment
  seqsIDs <- attributes(fasta)$dimnames[[1]] # How many individual sequences there are
  
  starts <- seq(1, alglen-winstep, by = winstep) # The starts of the windows
  n <- length(starts)    # Find the length of the vector "starts"
  
  windowpis <- c() # The vectors to save all the PIs per window
  windoweffectiven <- c()
  
  # Loop through each window and get pi
  for (win in 1:n){
    
    if (starts[win] + winstep > alglen - windowsize){ # If the next step moves into the last window
      winseq <- fasta[,starts[win]:alglen] %>% del.colgapsonly(threshold = 1/length(seqsIDs))
    } else {
      winseq <- fasta[,starts[win]:(starts[win] + windowsize)] %>% del.colgapsonly(threshold = 1/length(seqsIDs))
    }   
    
    windoweffectiven <- c(windoweffectiven, dim(winseq)[2])
    
    # Pi = (n/(n-1))*2*[(sum_i=1)^n (sum_j=1)^(i-1) (x_i*x_j*pi_ij) ] 
    # Nei and Li in 1979 plus sample size correction, and see https://en.wikipedia.org/wiki/Nucleotide_diversity
    
    # Make a list with pi per site for that window
    pis_window <- c()
    
    for (i in winseq %>% seg.sites){ # For each site in the window, using the segregating sites
      sitefreqs <- winseq[,i] %>% base.freq() # Extract nucleotide frequencies
      sitefreqs <- sitefreqs[sitefreqs > 0] # Keep only nucleotides with counts
      nu_bases_site <- winseq[,i] %>% dim() %>% .[1] # Sample size n
      
      # Calculate pi
      if (length(sitefreqs) <= 2){ # Use only biallelic sites
        if (sitefreqs[1] == 1){ # If site is not segregating
          pi_site <- 0
        } else { # Segregating sites
          pi_site <- prod(sitefreqs)*1*2*(nu_bases_site/(nu_bases_site-1))
        }
        pis_window <- c(pis_window, pi_site) # record in the final list
      }
    }
    
    thepi_window <- sum(pis_window)/dim(winseq)[2]  # The sum of all pis per site divided by the number of sites in window
    windowpis <- c(windowpis, thepi_window)
  }
  
  coords <- starts + (windowsize/2)
  ends <- c(starts[1:(length(starts)-1)] + windowsize, alglen)
  result <- data.frame(coords = coords, starts = starts, ends = ends, pi = windowpis, n = windoweffectiven)
  
  return(result)
}

# Get the IDs of the sequences
seqsIDs <- attributes(spoks)$dimnames[[1]]

# Subset to match only the representative of each spok homolog
references <- c("Spok2_PaSp", "Spok3_PaWa87m", "Spok4_PaWa87m", "Spok1_PcTp")
spoksref <- spoks[match(references, seqsIDs),]

# Subset by Spok
spok2ids <- c("Spok2_PaWa21m", "Spok2_PaWa28m", "Spok2_PaWa53m", "Spok2_PaWa58m", "Spok2_PaWa63p", "Spok2_PaWa87m", "Spok2_PaWa100p", "Spok2_PaSp")
spok3ids <- c("Spok3_PaWa21m", "Spok3_PaWa28m", "Spok3_PaWa58m", "Spok3_PaWa87m", "Spok3_PaYp", "Spok3_CBS237.71m", "Spok3_PaWa53m", "Spok3_PaTgp")
spok4ids <- c("Spok4_PaWa53m", "Spok4_PaWa58m", "Spok4_PaWa87m", "Spok4_PaWa100p", "Spok4_PaYp", "Spok4_PaTgp", "Spok4_CBS237.71m")

spok2 <- spoks[match(spok2ids, seqsIDs),]
spok3 <- spoks[match(spok3ids, seqsIDs),]
spok4 <- spoks[match(spok4ids, seqsIDs),]

# Make a nice data frame
spokwindow <- 100
spokstep <- 20

spokdf <- cbind(winpiseq_step(spoksref, spokwindow, spokstep), Spok = "All")
spokdf <- rbind(spokdf, cbind(winpiseq_step(spok2, spokwindow, spokstep), Spok = "Spok2"))
spokdf <- rbind(spokdf, cbind(winpiseq_step(spok3, spokwindow, spokstep), Spok = "Spok3"))
spokdf <- rbind(spokdf, cbind(winpiseq_step(spok4, spokwindow, spokstep), Spok = "Spok4"))

ggplot(spokdf, aes(x = coords, y = pi, colour = Spok)) + geom_line(alpha = 0.8) + #geom_point(alpha = 0.3) +
  xlab("Nucleotide position (bp)") + ylab("Pairwise nucleotide differences") + 
  scale_colour_manual(values = c("black", "chocolate1", "chartreuse4", "brown")) +
  theme_bw()

outputpath <- "PaperSpoks/Figures/"
ggsave(file = paste0(outputpath, "Pi_spoks.pdf"), width = 8, height = 1.5, scale = 3, units = "cm")
