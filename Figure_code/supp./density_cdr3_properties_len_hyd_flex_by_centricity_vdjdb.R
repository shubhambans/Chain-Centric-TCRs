###############################################################################
# END-TO-END SCRIPT — Metaclone Contrast × Chain (3×2 panel)
# Metrics: Length, Hydrophobicity, Flexibility
# Facets: metric (rows) × chain (columns)
###############################################################################

rm(list = ls())
suppressPackageStartupMessages({
    library(dplyr)
    library(stringr)
    library(tidyr)
    library(ggplot2)
    library(ggforce)
    library(Peptides)
    library(effsize)
	library(ggh4x) 
})

###############################################################################
# USER SETTINGS
###############################################################################

DATA_FILE    <- "~/Desktop/Shubhams_paper_code/data/hu_with_metaclone_20250826_105038.csv"

OUT_PDF <- "Metaclone_AA_properties.pdf"
setwd("~/Desktop/Shubhams_paper_code/Figure_PDFs")

# Contrast toggle (choose one)
PVAL_CONTRAST <- "alphaCentric_vs_none"
# options:
#   "alphaCentric_vs_none"
#   "betaCentric_vs_none"
#   "alphaCentric_vs_betaCentric"

###############################################################################
# LOAD + STANDARDIZE HU
###############################################################################

hu <- read.csv(DATA_FILE, check.names = FALSE, stringsAsFactors = FALSE)

hu <- hu %>%
    mutate(
        libid   = complex.id,
        junction = cdr3,
        v_gene_raw = v.segm,
        j_gene_raw = j.segm,
        Type = metaclone
    )

###############################################################################
# CLEAN METACLONE TYPE LABELS
###############################################################################

hu <- hu %>%
    mutate(
        Type = case_when(
            str_detect(tolower(Type), "alpha") ~ "alphaCentric",
            str_detect(tolower(Type), "beta")  ~ "betaCentric",
            TRUE                               ~ "none"
        )
    )

###############################################################################
# MOTIF FILTER (C…F/W)
###############################################################################

hu <- hu %>%
    filter(str_detect(junction, "^C")) %>%
    filter(str_detect(junction, "[FW]$"))

###############################################################################
# CHAIN INFERENCE
###############################################################################

strip_allele <- function(x) toupper(sub("\\*.*$", "", trimws(x)))

hu <- hu %>%
    mutate(
        v_gene = strip_allele(v_gene_raw),
        j_gene = strip_allele(j_gene_raw),
        chain = case_when(
            grepl("^TRAV", v_gene) | grepl("^TRAJ", j_gene) ~ "TRA",
            grepl("^TRBV", v_gene) | grepl("^TRBJ", j_gene) ~ "TRB",
            TRUE ~ NA_character_
        )
    ) %>%
    filter(chain %in% c("TRA", "TRB"))

###############################################################################
# REMOVE MAIT + iNKT + paired TRB
###############################################################################

mait_inkt_libs <- hu %>%
    filter(
        chain == "TRA" &
        (
            junction == "CVVSDRGSTLGRLYF" | 
            (v_gene == "TRAV1-2" & j_gene %in% c("TRAJ33","TRAJ20","TRAJ12"))
        )
    ) %>% pull(libid)

hu <- hu %>% filter(!(chain == "TRA" & libid %in% mait_inkt_libs))
hu <- hu %>% filter(!(chain == "TRB" & libid %in% mait_inkt_libs))

###############################################################################
# REQUIRE PAIRED TRA+TRB
###############################################################################

paired_ids <- hu %>%
    distinct(libid, chain) %>%
    count(libid, name="n_chain") %>%
    filter(n_chain >= 2) %>%
    pull(libid)

hu <- hu %>% filter(libid %in% paired_ids)

###############################################################################
# DEDUPLICATE AA JUNCTIONS PER LIBRARY
###############################################################################

df <- hu %>%
    distinct(libid, chain, junction, Type, .keep_all = TRUE)

###############################################################################
# COMPUTE AA-BASED METRICS
###############################################################################

# Junction length
df <- df %>% mutate(len = nchar(junction))

# Hydrophobicity (Eisenberg)
df <- df %>% mutate(hyd = Peptides::hydrophobicity(junction, scale="Eisenberg"))

# Flexibility (Vihinen VINM940101)
aa_flex <- c(
  A=0.984, R=1.008, N=1.048, D=1.068, C=0.906,
  Q=1.037, E=1.094, G=1.031, H=0.950, I=0.927,
  L=0.935, K=1.102, M=0.952, F=0.915, P=1.049,
  S=1.046, T=0.997, W=0.904, Y=0.929, V=0.931
)

vihinen_flex <- function(seq){
    aa <- strsplit(seq,"")[[1]]
    vals <- aa_flex[aa]
    if (any(is.na(vals))) return(NA_real_)
    mean(vals)
}

df$flex <- vapply(df$junction, vihinen_flex, numeric(1))

df <- df %>% filter(!is.na(hyd), !is.na(flex))

###############################################################################
# LONG-FORM PANEL
###############################################################################

panel_df <- df %>%
    select(chain, Type, len, hyd, flex) %>%
    pivot_longer(
        cols = c(len, hyd, flex),
        names_to = "metric",
        values_to = "value"
    )

panel_df$metric <- factor(
    panel_df$metric,
    levels = c("len","hyd","flex"),
    labels = c("Length","Hydrophobicity","Flexibility")
)

###############################################################################
# APPLY CONTRAST FILTER — keep exactly 2 metaclones
###############################################################################

contrast_pair <- switch(
    PVAL_CONTRAST,
    alphaCentric_vs_none        = c("alphaCentric", "none"),
    betaCentric_vs_none         = c("betaCentric",  "none"),
    alphaCentric_vs_betaCentric = c("alphaCentric", "betaCentric")
)

panel_df <- panel_df %>% filter(Type %in% contrast_pair)
df       <- df       %>% filter(Type %in% contrast_pair)

###############################################################################
# STATS
###############################################################################

stats_df <- panel_df %>%
    group_by(metric, chain) %>%
    summarise(
        g1 = list(value[Type == contrast_pair[1]]),
        g2 = list(value[Type == contrast_pair[2]]),
        .groups="drop"
    ) %>%
    mutate(
        ks_p = mapply(function(a,b)
            if (length(a)<2 || length(b)<2) NA_real_
            else ks.test(a,b)$p.value,
            g1, g2
        ),
        cliff = mapply(function(a,b)
            effsize::cliff.delta(a,b)$estimate,
            g1, g2
        ),
        p_adj = p.adjust(ks_p, "BH")
    )

stars <- function(p){
    if (is.na(p)) return("")
    if (p < 1e-4) "****"
    else if (p < 1e-3) "***"
    else if (p < 1e-2) "**"
    else if (p < 0.05) "*"
    else "NS"
}

stats_df$label <- paste0(
    sapply(stats_df$p_adj, stars),
    "\nD=", sprintf("%.2f", stats_df$cliff)
)

###############################################################################
# MEDIANS
###############################################################################

med_df <- panel_df %>%
    group_by(metric, chain, Type) %>%
    summarise(med = median(value), .groups="drop")

###############################################################################
# THEME: LARGE FONTS
###############################################################################

theme_big <- theme_bw(base_size = 24) +
    theme(
        legend.position = "top",
        legend.title    = element_text(size = 22, face="bold"),
        legend.text     = element_text(size = 20),
        axis.text       = element_text(size = 18),
        axis.title      = element_text(size = 24, face="bold"),

        strip.text.x    = element_text(size = 24, face="bold"),
        strip.text.y    = element_text(size = 24, face="bold"),
        strip.background = element_blank(),   # <<< remove gray shading

        plot.title      = element_text(size = 30, face="bold", hjust=0.5),
        panel.grid      = element_blank()
    )

###############################################################################
# PLOT
###############################################################################

p <- ggplot(panel_df, aes(x=value, fill=Type)) +
    geom_density(alpha=0.45) +

    facet_grid2(
        metric ~ chain,
        scales = "free",
        independent = "x"
    ) +

    geom_vline(
        data=med_df,
        aes(xintercept=med, color=Type),
        size=0.7,
        show.legend=FALSE
    ) +

    geom_text(
        data=stats_df,
        aes(label=label, x = -Inf, y = Inf),
        hjust = -0.1, vjust = 1.1,
        size = 6,
        inherit.aes = FALSE
    ) +

    scale_fill_manual(values=c(
        alphaCentric="#66c2a5",
        betaCentric="gray40",
        none="#fc8d62"
    )) +
    scale_color_manual(values=c(
        alphaCentric="#66c2a5",
        betaCentric="gray40",
        none="#fc8d62"
    )) +

    labs(
        x = "",
        y = "Density",
        title = "TRA vs TRB — Length / Hydrophobicity / Flexibility"
    ) +
    theme_big

ggsave(OUT_PDF, p, width=12, height=12, units="in")
print(p)

###############################################################################
# END
###############################################################################
