###############################################################################
# END-TO-END: Su et al CD4 — Fraction of metaclones by disease stage (NHC)
# - Removes MAIT / iNKT by libid (TRA heuristics)
# - Filters cohort: study == "su"
# - Defines GLOBAL metaclones from paired TRA/TRB within cohort
# - Restricts to Collection_Day ∈ {BL, single_tp}
# - Computes fraction = n_metaclone / n_total per donorVisit × NHC × type
# - Boxplot + jitter + per-facet stats (ggpubr/rstatix)
# - Saves PDF + PNG to OUT_DIR
###############################################################################

rm(list = ls())

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(ggplot2)
  library(ggpubr)
  library(rstatix)
  library(PropCIs)
  library(confintr)
})

###############################################################################
# 0) OUTPUT DIRECTORY + FILENAMES
###############################################################################
OUT_DIR <- file.path(getwd(), "Figure_Images")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

OUT_STEM <- "Fig.3E_Su_CD4_fraction_metaclone_by_stage"
OUT_PDF  <- file.path(OUT_DIR, paste0(OUT_STEM, ".pdf"))
OUT_PNG  <- file.path(OUT_DIR, paste0(OUT_STEM, ".png"))

###############################################################################
# 1) GLOBAL THEME (edit sizes here)
###############################################################################
BASE_SIZE  <- 28
TITLE_SIZE <- 32
AXIS_TTL   <- 30
AXIS_TXT   <- 26

theme_set(
  theme_bw(base_size = BASE_SIZE) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.key       = element_blank(),
      plot.title       = element_text(face = "bold", size = TITLE_SIZE, hjust = 0.5),
      axis.title.x     = element_text(face = "bold", size = AXIS_TTL),
      axis.title.y     = element_text(face = "bold", size = AXIS_TTL),
      axis.text.x      = element_text(size = AXIS_TXT),
      axis.text.y      = element_text(size = AXIS_TXT)
    )
)

###############################################################################
# 2) LOAD DATA
###############################################################################
CD4merge <- readRDS("~/Desktop/Shubhams_paper_code/data/CD4merge07212025.rds")

# ---- Required columns (fail fast) ----
req_cols <- c("libid","junction","v_gene","j_gene","chain","study",
              "donorVisit","NHC_Scalefactor","Collection_Day")
miss <- setdiff(req_cols, names(CD4merge))
if (length(miss)) stop("Missing required columns in CD4merge: ", paste(miss, collapse = ", "))

###############################################################################
# 3) REMOVE iNKT / MAIT (by libid; TRA heuristics)
###############################################################################
iNkt_libs <- CD4merge %>%
  filter(junction == "CVVSDRGSTLGRLYF") %>%
  distinct(libid) %>%
  pull(libid)

mait_libs <- CD4merge %>%
  filter(v_gene == "TRAV1-2", j_gene %in% c("TRAJ33","TRAJ20","TRAJ12")) %>%
  distinct(libid) %>%
  pull(libid)

plot_base <- CD4merge %>%
  filter(!(libid %in% union(iNkt_libs, mait_libs)))

###############################################################################
# 4) COHORT FILTER
###############################################################################
plot_base <- plot_base %>%
  filter(study == "su")

###############################################################################
# 5) GLOBAL METACLONE CLASSIFICATION (within cohort)
# alphaCentric: TRA junction pairs with >1 distinct TRB junction
# betaCentric : TRB junction pairs with >1 distinct TRA junction
###############################################################################
tra_all <- plot_base %>% filter(chain == "TRA")
trb_all <- plot_base %>% filter(chain == "TRB")

paired_all <- inner_join(
  tra_all, trb_all, by = "libid", suffix = c(".tra", ".trb")
)

alpha_meta <- paired_all %>%
  group_by(junction.tra) %>%
  filter(n_distinct(junction.trb) > 1) %>%
  ungroup() %>%
  mutate(type = "alphaCentric")

beta_meta <- paired_all %>%
  group_by(junction.trb) %>%
  filter(n_distinct(junction.tra) > 1) %>%
  ungroup() %>%
  mutate(type = "betaCentric")

metaclones_global <- bind_rows(alpha_meta, beta_meta) %>%
  distinct(libid, type, .keep_all = TRUE) %>%
  transmute(
    libid,
    donorVisit = donorVisit.tra,
    NHC        = NHC_Scalefactor.tra,
    type
  )

###############################################################################
# 6) RESTRICT TO SINGLE VISIT (BL, single_tp)
###############################################################################
plot_base <- plot_base %>%
  filter(Collection_Day %in% c("BL", "single_tp"))

###############################################################################
# 7) COUNT METACLONES + TOTAL CELLS AND COMPUTE FRACTIONS
###############################################################################
plot_df <- metaclones_global %>%
  semi_join(plot_base, by = "libid") %>%
  count(donorVisit, NHC, type, name = "n_metaclone")

total_cells <- plot_base %>%
  count(donorVisit, NHC = NHC_Scalefactor, name = "n_total")

plot_df <- plot_df %>%
  complete(donorVisit, NHC, type = c("alphaCentric", "betaCentric"),
           fill = list(n_metaclone = 0)) %>%
  left_join(total_cells, by = c("donorVisit", "NHC")) %>%
  mutate(fraction = n_metaclone / n_total) %>%
  filter(!is.na(fraction))

###############################################################################
# 8) ORDERING: Healthy→Mild→Moderate→Severe; within each, sort by |α-β|
###############################################################################
clin_levels <- c("Healthy", "Mild", "Moderate", "Severe")
plot_df$NHC <- factor(plot_df$NHC, levels = clin_levels)

ordering_df <- plot_df %>%
  pivot_wider(names_from = type, values_from = fraction,
              values_fill = list(fraction = 0)) %>%
  mutate(diff = abs(alphaCentric - betaCentric)) %>%
  arrange(NHC, desc(diff)) %>%
  distinct(donorVisit, .keep_all = TRUE)

plot_df$donorVisit <- factor(plot_df$donorVisit, levels = ordering_df$donorVisit)

###############################################################################
# 9) OPTIONAL: CI COMPUTATION (example shown for Healthy alpha/beta)
###############################################################################
sub_hc  <- plot_df %>% filter(NHC == "Healthy")
sub_hcA <- sub_hc %>% filter(type == "alphaCentric")
sub_hcB <- sub_hc %>% filter(type == "betaCentric")

# Bootstrap CI of the mean (can be slow at 99999; keep as-is if you want)
ciAlpha <- confintr::ci_mean(
  sub_hcA$fraction,
  probs = c(0.0005, 0.9959),
  type  = "bootstrap",
  R     = 99999L,
  seed  = 42
)

ciBeta <- confintr::ci_mean(
  sub_hcB$fraction,
  probs = c(0.0005, 0.9995),
  type  = "bootstrap",
  R     = 99999L,
  seed  = 42
)

cat("\nMean fractions (all):\n")
cat("alphaCentric mean =", mean(plot_df$fraction[plot_df$type=="alphaCentric"]), "\n")
cat("betaCentric  mean =", mean(plot_df$fraction[plot_df$type=="betaCentric"]), "\n")
cat("\nMean fractions (Healthy):\n")
cat("alphaCentric mean =", mean(sub_hcA$fraction), "\n")
cat("betaCentric  mean =", mean(sub_hcB$fraction), "\n")
cat("\nCI Alpha (Healthy):\n"); print(ciAlpha)
cat("\nCI Beta  (Healthy):\n"); print(ciBeta)

###############################################################################
# 10) STATS FOR ANNOTATION (per type facet)
###############################################################################
plot_df$NHC <- factor(plot_df$NHC, levels = clin_levels)
comps <- list(c("Healthy","Moderate"))  # edit/add more comparisons if desired

stat.test <- plot_df %>%
  group_by(type) %>%
  rstatix::t_test(fraction ~ NHC, comparisons = comps) %>%
  add_significance("p") %>%
  add_xy_position(x = "NHC", dodge = 0.8, step.increase = 0.06)

###############################################################################
# 11) PLOT
###############################################################################
p <- ggplot(plot_df, aes(x = NHC, y = fraction)) +
  geom_boxplot(outlier.shape = NA, fill = "gray") +
  facet_wrap(~ type) +
  geom_jitter(color = "black", width = 0.2, size = 2, alpha = 0.5) +
  ggpubr::stat_pvalue_manual(
    stat.test,
    label = "p.signif",
    hide.ns = TRUE,
    tip.length = 0.01,
    size = 10
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.2))) +
  coord_cartesian(clip = "off") +
  labs(x = "\nDisease stage", y = "Fraction metaclone\n") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(p)

###############################################################################
# 12) SAVE (PDF + PNG)
###############################################################################
ggsave(
  filename = OUT_PDF,
  plot     = p,
  width    = 10,
  height   = 8,
  units    = "in",
  device   = cairo_pdf,
  bg       = "white"
)

ggsave(
  filename = OUT_PNG,
  plot     = p,
  width    = 10,
  height   = 8,
  units    = "in",
  dpi      = 300,
  bg       = "white"
)

cat("\nSaved files:\n", OUT_PDF, "\n", OUT_PNG, "\n", sep = "")

###############################################################################
# END
###############################################################################
