library(tidyverse)
library(cowplot)
theme_set(theme_cowplot())

sim_plot <- plot_grid(plot1, plot2, plot3, plot4, plot5, plot6, plot7, plot8, NULL, plot9, plot10, NULL, nrow = 3)

save_plot(
  "sim_plot.pdf", sim_plot,
  ncol = 4, 
  nrow = 2, 
  base_asp = 1.35 
)
