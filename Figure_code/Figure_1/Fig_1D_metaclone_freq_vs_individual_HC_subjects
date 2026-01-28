rm(list = ls())
# ── Packages ────────────────────────────────────────────────────────
library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(ggsignif)
library(PropCIs)
library(confintr)

# ── Theming ─────────────────────────────────────────────────────────
rm(list = ls())
theme_set(
    theme_bw(base_size = 28) +
        theme(
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            legend.key = element_blank(),
            plot.title = element_text(face = "bold", size = 32, hjust = 0.5),
            axis.title.x = element_text(face = "bold", size = 30),
            axis.title.y = element_text(face = "bold", size = 30),
            axis.text.x = element_text(size = 26),
            axis.text.y = element_text(size = 26)
			        )
			)

# ── load Su et al data ────────────────────────────────────────────────────────
CD4merge <- readRDS("~/Desktop/Shubhams_paper_code/data/CD4merge07212025.rds") #
OUTPUT_DIR   <- "~/Desktop/Shubhams_paper_code/Figure_PDFs"

## remove iNKT and MAIT cell sequencces

iNkt1 = subset(CD4merge, junction == "CVVSDRGSTLGRLYF")
iNkt2 = subset(CD4merge, libid %in% iNkt1$libid)
mait1 = subset(CD4merge, v_gene == "TRAV1-2" & (j_gene == "TRAJ33" | j_gene == "TRAJ20" | j_gene == "TRAJ12")) # 18
mait2 = subset(CD4merge, libid %in% mait1$libid) # 45

tcrsCombSub3 = subset(CD4merge, !junction %in% iNkt2$junction) # 
tcrsCombSub3 = subset(tcrsCombSub3, !junction %in% mait2$junction) # 7633

# ── 0. Cohort filter (optional) ────────────────────────────────────
plot_base <- CD4merge %>% 
  dplyr::filter(study == "su") # ← adjust / remove as needed
  
# ── 1. GLOBAL metaclone classification on full CD4merge ────────────
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
        NHC = NHC_Scalefactor.tra,
        type
    )
# ── 2. Count metaclone & total cells in the chosen cohort ──────────
plot_df <- metaclones_global %>% 
    semi_join(plot_base, by = "libid") %>% 
    dplyr::count(donorVisit, NHC, type, name = "n_metaclone")
total_cells <- plot_base %>% 
    dplyr::count(donorVisit, NHC = NHC_Scalefactor, name = "n_total")
# ── 3. Fill missing α/β combos with 0 & compute fractions ──────────
plot_df <- plot_df %>% 
    complete(donorVisit, NHC, type = c("alphaCentric", "betaCentric"),
             fill = list(n_metaclone = 0)) %>% 
    left_join(total_cells, by = c("donorVisit", "NHC")) %>% 
    mutate(fraction = n_metaclone / n_total)
# ── 4. Build x-axis order: Healthy → Mild → Moderate → Severe, │
# then descending |α – β| difference within each group. ──
clin_levels <- c("Healthy", "Mild", "Moderate", "Severe")
plot_df <- plot_df[!is.na(plot_df$fraction), ]
ordering_df <- plot_df %>% 
    mutate(NHC = factor(NHC, levels = clin_levels)) %>% 
    pivot_wider(
        names_from = type,
        values_from = fraction,
        values_fill = list(fraction = 0) # ← named list, not bare 0
    ) %>% 
    mutate(diff = abs(alphaCentric - betaCentric)) %>% 
    arrange(NHC, desc(diff)) %>% 
    distinct(donorVisit, .keep_all = TRUE)
plot_df$donorVisit <- factor(plot_df$donorVisit,
                             levels = ordering_df$donorVisit)
plot_df <- plot_df[!is.na(plot_df$fraction), ]

# ── subset and calculations ──────────────────────────────────────────────

sub1 = subset(plot_df, NHC == "Healthy")
sub2 = subset(plot_df, NHC == "Mild")
sub3 = subset(plot_df, NHC == "Moderate")
sub4 = subset(plot_df, NHC == "Severe")

aSub = subset(plot_df, type == "alphaCentric")
bSub = subset(plot_df, type == "betaCentric")

sub1a = subset(sub1, type =="alphaCentric")
sub1b = subset(sub1, type =="betaCentric")

mean(aSub$fraction) # 0.1469291
mean(bSub$fraction) # 0.03599285
# 0.1469291/0.03599285
#[1] 4.082175
# NB, ~4x more alpha centric than beta centric

mean(sub1a$fraction) # 0.1402561
mean(sub1b$fraction) # 0.02787892
# 0.1402561/0.02787892
#[1] 5.030901
# NB, ~5x more alpha centric than beta centric in HC

# model
#lm_model <- lm(fraction ~ NHC*type, data = plot_df)
#anova(lm_model)

# Compute confidence intervals 

ciAlpha = ci_mean(
sub1a$fraction,
probs = c(0.0005, 0.9959),
type = c("bootstrap"),
R = 99999L,
seed = 42,
)

ciBeta = ci_mean(
sub1b$fraction,
probs = c(0.0005, 0.9995),
type = c("bootstrap"),
R = 99999L,
seed = 42,
)

	
# ── 7. dot-plot HC only ──────────────────────────────────────────────

if (dev.cur() > 1) dev.off()
quartz(height = 8, width = 14, dpi = 72)

# Set default theme and point size
theme_set(
  theme_bw(36) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.key = element_blank()
    )
)
update_geom_defaults("point", aes(size = 8))

cbPalette = c('#66c2a5','#fc8d62', 'gray')

# Plot
ggplot(subset(plot_df, NHC == "Healthy"), aes(donorVisit, fraction, colour = type)) +
  geom_point(size = 6) +
  scale_y_continuous(limits = c(0, 0.2), breaks = seq(0, 0.20, 0.05)) +
  scale_x_discrete(labels = seq_along(levels(plot_df$donorVisit))) +
  ylab("Fraction TCRs in metaclone\n") +
  xlab("\nSample number") +
  theme(axis.text.x = element_text(angle = 0, hjust = 1)) +
  scale_color_manual(values = cbPalette) 
	
## save plot

setwd(OUTPUT_DIR)

toSave <- subset(plot_df, NHC == "Healthy")
filename <- "Fig1D_metaclone_frequencies_HC.pdf"
ggsave(filename)
