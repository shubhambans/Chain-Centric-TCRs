#hello

library(ggplot2)
library(dplyr)
library(tidyr)
librayr(tibble)
library(igraph)
#library(stringDist)
library(ggsignif)
library(Seurat)
library(monocle3)
library(ggraph2)

theme_set(
  theme_bw(base_size = 36) + 
  theme(
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    legend.key = element_blank(),        # Remove the legend background
    plot.title = element_text(face = "bold", size = 40, hjust = 0.5),  # Bold and larger plot title, centered
    axis.title.x = element_text(face = "bold", size = 38),  # Bold and larger x-axis title
    axis.title.y = element_text(face = "bold", size = 38),  # Bold and larger y-axis title
    axis.text.x = element_text(size = 34),  # Increase size of x-axis labels
    axis.text.y = element_text(size = 34)   # Increase size of y-axis labels
  )
)

