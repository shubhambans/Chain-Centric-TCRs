#hello
#This is basic setup
#Adttional figure setup code can be found in other Set up files

library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)
library(igraph)
#library(stringDist)
library(ggsignif)
library(Seurat)
library(monocle3)
library(ggplot2)

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




#load in files (pre-existing)
#files can be found in the GIT folder in BOX
load("C:/Users/sbansal/Downloads/vdjdb09122024ALLHuman.RData")
load("C:/Users/sbansal/Downloads/finalfiguresData.RData")
#load("C:/Users/super/OneDrive/Desktop/Research-Surface/vdjdb09122024ALLHuman.RData")
#load("C:/Users/super/OneDrive/Desktop/Research-Surface/finalfiguresData.RData")
