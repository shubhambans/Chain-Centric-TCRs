#hello
#This is basic setup
#Additional figure setup code can be found in other Set up files
#
# -----------------------------------------------------------------------------
# Workspace Objects Key (loaded from the workspace/RData files below)
# -----------------------------------------------------------------------------
# CD4merge
#   Cell level: all CD4 T cells (master dataframe). Each row is one CD4+ cell
#   with junctions, V/D/J genes, expansion tag, donor, study-of-origin, etc.
# CD4merge_filtered
#   Cell level: expanded CD4 cells only (subset of CD4merge). Not used in figures.
# CD4merge_unexpanded
#   Cell level: unexpanded CD4 cells (complement of filtered). Not used in figures.
# data_TRA / data_TRB
#   Cell level: alpha-/beta-chain slices of CD4merge with identical metadata.
# geneseuratcd4PossibleMistake
#   Seurat object: raw Su et al. CD4 import (~4.8GB).
# indices_combined
#   Donor-level CD4 metrics aggregated from CD4merge (one row per donor).
# indices_combined_filtered / indices_combined_unexpanded
#   Donor-level metrics from filtered/unexpanded CD4 subsets. Not used in figures.
# indices_TRA / indices_TRB
#   Donor-level metrics for TRA-only / TRB-only slices of indices_combined.
# merged_df
#   All T-cell subsets across studies (CD4, CD8, γδ, etc.) for cross-subset comparisons.
# percent_public_df
#   Per-donor % public TCRs summary (needs verification). Not used in figures.
# RealBatchCorrectedBacherSeurat
#   Seurat object: Harmony-corrected Bacher et al. dataset (~7.2GB).
# vdjdb09122024ALLHuman
#   VDJdb human antigen-specific TCRs (downloaded 09/12/2024).
# MetadataForRealBacherSeurat
#   Metadata for the batch-corrected Bacher Seurat object.
# combined_seurat_object_bacher
#   Non-batch-corrected Bacher Seurat object.
# cds_full_monocle3_object_bacher
#   Monocle3 object derived from Bacher data (used to calculate pseudotime).
#
# -----------------------------------------------------------------------------
# Data set sources used in this study (Table 1)
# -----------------------------------------------------------------------------
# Set 1: VDJdb (human antigen-specific TCRs)
# Set 2: Bacher et al. (SARS-CoV-2–enriched T cells)
# Set 3: Su et al. (healthy controls)
# Set 4: Su et al. (COVID-19 + healthy controls)
# Set 5: Combined datasets (COVID-19 + healthy controls)
# -----------------------------------------------------------------------------

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




# ---------------------------------------------------------------------------
# Load data files
# ---------------------------------------------------------------------------
# Download the data bundle from Box:
# https://bri.box.com/s/ny65o2w9dukpt5ejvmc7nl1hl26ubicl
# Then update the paths below to point at the downloaded `.RData` files.
# ---------------------------------------------------------------------------
#load in files (pre-existing)
#files can be found in the GIT folder in BOX
load("C:/Users/sbansal/Downloads/vdjdb09122024ALLHuman.RData")
load("C:/Users/sbansal/Downloads/finalfiguresData.RData")
#load("C:/Users/super/OneDrive/Desktop/Research-Surface/vdjdb09122024ALLHuman.RData")
#load("C:/Users/super/OneDrive/Desktop/Research-Surface/finalfiguresData.RData")
