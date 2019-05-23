library(tidyverse)
library(stringr)
library(cowplot)
library(broom)
library(ggseqlogo)

# read in data
alignment <-read_lines("data/simulated/results_1B4T_A_evolved.txt") 
partial_align <- substring(alignment, 1, 10) 

ggseqlogo(partial_align, method = 'prob') +
  theme(legend.position = "none")
