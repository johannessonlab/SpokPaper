#!/usr/bin/env Rscript

# Dotplot of the haplotype of Pseudospok1 vs the Spok block
# =======================================
# NUCmer was run like such:
# $ nucmer -b 2000 -c 20 -p pseudospok.delta Psk1_PaWa87p.fa Pseudospok1_PaWa87p.fa -t 1 --maxmatch
# $ show-coords -r -B pseudospok.delta > pseudospok.tab
# =======================================
# Sandra Lorena Ament Velasquez
# Johannesson Lab, Evolutionary Biology Center, Uppsala University, Sweden
# 2018-11-09
# =======================================
Version <- 1.0
# ============================
# Load the necessary libraries
# ============================
library(ggplot2, quietly = TRUE)
library(grid)
library(dplyr, warn.conflicts = FALSE)

file <- "pseudospok.tab"
coords <- read.delim(file, header = FALSE)
names(coords) <- c("QUERY","DATE","LEN_Q","PROGRAM","REF_FILE","REF","S2","E2","S1","E1","IDY","SYM","LEN_2","V14","V15","V16","V17","SENSE","LEN_R","V20","V21")

# ============================
# Processing
# ============================
# Prepare a function to produce the DotPlots from a given data set

PaWa87p <- data.frame(Gene = c("Spok3", "Spok4"),
                      S1 = c(82737, 40494),
                      end = c(84974, 42770),
                      LEN_Q = 0)

pseudo <- data.frame(Gene = c("Pseudospok1", "Pseudospok1"),
                     S1 = c(0, 0),
                     end = c(7716, 10582),
                     LEN_Q = c(7198, 9581))

dotplotter <- function(data, multiref = FALSE){
  # Make an empty plot, with panels of height equal to lenght of queries
  if (multiref){
    p <- ggplot(data, aes(x=S1, y=LEN_Q)) + geom_blank() + facet_grid(QUERY ~ REF, scales="free", space = "free") + expand_limits(y=0) +
      theme(panel.margin = unit(0.1, "lines")) # To put the facets closer together
  } else{
    p <- ggplot(data, aes(x=S1, y=LEN_Q)) + geom_blank() + facet_grid(QUERY ~ ., scales="free_y", space = "free_y") + expand_limits(y=0) # expand_limits(y=0) is to fix the lower bound of the y axis to 0
  }
  
  # Label axes
  if ( length(unique(data$QUERY)) == 1){ # If there is only one query
    query <- as.character(data$QUERY[1])
    p <- p + ylab(query) + ylim(0, max(data$LEN_Q)) # Make axis the length of that query sequence
  } else {
    p <- p + ylab("Query sequences")
  }
  if ( length(unique(data$REF)) == 1){
    reference <- as.character(data$REF[1])
    p <- p + xlab(reference) + xlim(0, max(data$LEN_R)) # Make axis the length of that reference sequence
  } else{
    p <- p + xlab("Reference sequences")
  }
  
  # Plot the alignments 
  for (i in 1:dim(data)[1]){
    currentrow <- data[i,]
    if (currentrow$S2 > currentrow$E2){ # Reverse
      colorsito <- "cyan3"
    }
    else{ # Forward
      colorsito <- "coral2"
    }
    p <- p + geom_segment(data = currentrow, aes(x = S1, xend = E1, y = S2, yend = E2), colour = colorsito) +
      geom_point(data = currentrow, aes(x=S1, y = S2), colour = colorsito, size = 1) + geom_point(data = currentrow, aes(x=E1, y = E2), colour = colorsito, size = 1)
  }
  return(p)
}

# ============================
# Plotting
# ============================

dotplotter(coords) + theme_bw() + # White background
  theme(strip.text.y = element_blank(), panel.border=element_blank()) + # Remove the gray label on the side and border of plot
  geom_hline(yintercept=0) + geom_vline(xintercept=0) + 
  geom_segment(data = PaWa87p, aes(x = S1, xend = end, yend = 1, colour = Gene), size = 5) +
  geom_segment(data = pseudo, aes(x = S1, xend = 1, y = LEN_Q, yend = end, colour = Gene), size = 5) +
  scale_color_manual(values=c("azure4", "chartreuse4", "brown2")) 

ggsave(file = "Figure_S4_PseudospokBlock.pdf", width = 8, height = 5, scale = 2.5, units = "cm")

