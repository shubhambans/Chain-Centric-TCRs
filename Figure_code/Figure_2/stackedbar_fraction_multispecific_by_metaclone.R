rm(list = ls())

## =========================
## Fraction of multi-epitope by metaclone
## - Toggle UNIT: "libid" (default) or "cdr3"
## - "multiple" is defined on the CHAIN_MODE CDR3 globally (>1 unique epitopes)
## - Counting/plotting/p-values are done per UNIT
## - PAIRED_ONLY filters libids before building UNITs
## - PLOT_THIRD_GROUP enforces which of {"none","bothCentric"} is kept as 3rd bar
## - Saves BOTH PNG + PDF to OUT_DIR
## =========================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(scales)
  library(ggpubr)
  library(rlang)
})

## ---------- OUTPUT DIRECTORY ----------
OUT_DIR <- file.path(getwd(), "Figure_Images")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

## ---------- GLOBAL FONT / SIZE CONTROLS ----------
BASE_SIZE_PT   <- 34
STAR_SIZE_PT   <- 32
TEXT_SIZE_PT   <- 12
BRACKET_LW     <- 0.6
BAR_WIDTH      <- 0.85

AXIS_X_ANGLE   <- 60
PLOT_W_IN      <- 10
PLOT_H_IN      <- 8
#PNG_DPI        <- 300

## ---------- TOGGLES ----------
UNIT <- "cdr3"                 # "libid" or "cdr3"
CHAIN_MODE <- "TRA"             # "TRA" or "TRB"
PAIRED_ONLY <- TRUE             # require both chains per libid (applies before UNIT build)
PLOT_THIRD_GROUP <- "none"      # "none" or "bothCentric"

PVAL_FROM_PRE_DOWNSAMPLED <- TRUE
PVAL_USE_STARS <- TRUE
PVAL_SHOW_NS   <- TRUE
ALPHA_NS       <- 0.05

Y_START  <- 1.04
Y_STEP   <- 0.08
Y_MARGIN <- 0.02
EXTRA_Y_HEADROOM <- 0.10

X_ORDER   <- c("alphaCentric","betaCentric","none","bothCentric")
LVL_ORDER <- c("alphaCentric","betaCentric","bothCentric","none")

set.seed(42)

## ---------- INPUT ----------
in_csv <- "~/Desktop/Shubhams_paper_code/data/hu_with_metaclone_20250826_095644.csv"
tot <- read.csv(in_csv, check.names = FALSE, stringsAsFactors = FALSE)

## ---------- OUTPUT FILENAMES ----------
stem <- sprintf("fraction_multiple_%s_%s%s_%s_plotThird-%s",
                UNIT, CHAIN_MODE, if (PAIRED_ONLY) "_pairedonly" else "",
                if (PVAL_USE_STARS) "stars" else "ptext",
                PLOT_THIRD_GROUP)
#OUTFILE_PNG <- file.path(OUT_DIR, paste0(stem, ".png"))
OUTFILE_PDF <- file.path(OUT_DIR, paste0(stem, ".pdf"))

## ---------- Ensure libid exists ----------
if ("complex.id" %in% names(tot)) tot$libid <- tot$complex.id
if (!"libid" %in% names(tot)) stop("Need 'complex.id' (or precomputed 'libid').")

## ---------- Detect CDR3 column ----------
cdr3_candidates <- c("cdr3","junction","junction_aa","cdr3_aa")
cdr3_col <- cdr3_candidates[cdr3_candidates %in% names(tot)][1]
if (is.na(cdr3_col)) stop("No CDR3 column found among: ", paste(cdr3_candidates, collapse = ", "))

## ---------- Column checks ----------
if (!("metaclone" %in% names(tot))) stop("Missing 'metaclone' column.")

epitope_col <- intersect(c("antigen.epitope","epitope","epitope.species"), names(tot))
if (!length(epitope_col)) stop("Need an epitope column: 'antigen.epitope'/'epitope'/'epitope.species'.")
epitope_col <- epitope_col[1]

v_col <- intersect(c("v.segm","v_gene","V.GENE.and.allele","v_gene_raw"), names(tot))[1]
if (is.na(v_col)) stop("No V-gene column found.")

## Normalize metaclone spelling (typo guard)
tot$metaclone <- gsub("betaCEntric","betaCentric", tot$metaclone, fixed = TRUE)

## ---------- Chain tagging ----------
tot_all <- tot %>%
  mutate(
    v_clean = toupper(sub("^([A-Z0-9-]+).*", "\\1", .data[[v_col]])),
    chain   = case_when(
      grepl("^TRB", v_clean) ~ "TRB",
      grepl("^TRA", v_clean) ~ "TRA",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(chain))

## ---------- PAIRED filter (before building UNITs) ----------
if (PAIRED_ONLY) {
  paired_ids <- tot_all %>%
    distinct(libid, chain) %>%
    count(libid, name = "n_chain") %>%
    filter(n_chain >= 2) %>%
    pull(libid)
  tot_all <- tot_all %>% filter(libid %in% paired_ids)
}

## ---------- GLOBAL CDR3 -> EPITOPE COUNTS (CHAIN_MODE only) ----------
cdr3_sym <- rlang::sym(cdr3_col)
epi_sym  <- rlang::sym(epitope_col)

cdr3_ep_counts <- tot_all %>%
  filter(chain == CHAIN_MODE) %>%
  filter(!is.na(!!epi_sym), !!epi_sym != "") %>%
  distinct(!!cdr3_sym, !!epi_sym) %>%
  count(!!cdr3_sym, name = "n_epitopes_cdr3")

## ---------- Helper: enforce third group choice ----------
apply_third_choice <- function(df, third_choice = c("none","bothCentric")) {
  third_choice <- match.arg(third_choice)
  m <- as.character(df$metaclone)
  avail <- unique(m)

  if (all(c("none","bothCentric") %in% avail)) {
    keep_mask <- m %in% c("alphaCentric","betaCentric", third_choice)
    df <- df[keep_mask, , drop = FALSE]
    m <- as.character(df$metaclone)
  }

  existing <- intersect(c("none","bothCentric"), unique(m))
  if (length(existing) == 1 && existing != third_choice) {
    m[m == existing] <- third_choice
  }

  df$metaclone <- factor(m,
                         levels = intersect(c("alphaCentric","betaCentric", third_choice), unique(m)),
                         ordered = TRUE)
  df
}

## ---------- Build UNIT table ----------
if (UNIT == "libid") {

  libid_map <- tot_all %>%
    filter(chain == CHAIN_MODE) %>%
    distinct(libid, !!cdr3_sym, metaclone)

  dups <- libid_map %>% count(libid) %>% filter(n > 1)
  if (nrow(dups)) {
    message("Warning: multiple CHAIN_MODE rows for some libids; taking first.")
    libid_map <- libid_map %>% group_by(libid) %>% slice_head(n = 1) %>% ungroup()
  }

  unit_df <- libid_map %>%
    left_join(cdr3_ep_counts, by = setNames(cdr3_col, cdr3_col)) %>%
    mutate(n_epitopes_cdr3 = tidyr::replace_na(n_epitopes_cdr3, 0L)) %>%
    filter(n_epitopes_cdr3 >= 1) %>%
    transmute(
      unit_id   = libid,
      metaclone = as.character(metaclone),
      multiple  = ifelse(n_epitopes_cdr3 > 1, "multiple", "one")
    )

  y_label <- "Fraction of libids (by CDR3 multi-epitope)"

} else if (UNIT == "cdr3") {

  cdr3_meta <- tot_all %>%
    filter(chain == CHAIN_MODE) %>%
    count(!!cdr3_sym, metaclone, name = "n") %>%
    group_by(!!cdr3_sym) %>%
    slice_max(n, with_ties = FALSE) %>%
    ungroup() %>%
    select(!!cdr3_sym, metaclone)

  unit_df <- cdr3_meta %>%
    left_join(cdr3_ep_counts, by = setNames(cdr3_col, cdr3_col)) %>%
    mutate(n_epitopes_cdr3 = tidyr::replace_na(n_epitopes_cdr3, 0L)) %>%
    filter(n_epitopes_cdr3 >= 1) %>%
    transmute(
      unit_id   = !!cdr3_sym,
      metaclone = as.character(metaclone),
      multiple  = ifelse(n_epitopes_cdr3 > 1, "multiple", "one")
    )

  y_label <- "Fraction of junctions\n"

} else {
  stop("UNIT must be 'libid' or 'cdr3'.")
}

## ---------- Enforce third-group choice & factor levels ----------
unit_df <- apply_third_choice(unit_df, third_choice = PLOT_THIRD_GROUP)
present_levels <- levels(unit_df$metaclone)
if (!length(present_levels)) stop("No metaclone groups present after filtering.")
if (nrow(unit_df) == 0) stop("No rows remain after classification and third-group enforcement.")

analysis_stats <- unit_df

## ---------- Down-sample to smallest group (plot only) ----------
sizes    <- unit_df %>% count(metaclone, name = "n_group")
target_n <- min(sizes$n_group, na.rm = TRUE)

analysis_ds <- unit_df %>%
  group_by(metaclone) %>%
  group_modify(~ dplyr::slice_sample(.x, n = min(target_n, nrow(.x)), replace = FALSE)) %>%
  ungroup()

## ---------- Plot data ----------
plot_df <- analysis_ds %>%
  mutate(multiple = factor(multiple, levels = c("multiple","one"))) %>%
  count(metaclone, multiple, name = "n_unit") %>%
  tidyr::complete(metaclone, multiple, fill = list(n_unit = 0)) %>%
  mutate(metaclone = factor(as.character(metaclone), levels = present_levels, ordered = TRUE))

plot_levels_display <- levels(plot_df$metaclone)

## ---------- Helpers ----------
to_wide <- function(df_long) {
  df_long %>%
    mutate(
      multiple  = factor(multiple, levels = c("multiple","one")),
      metaclone = as.character(metaclone)
    ) %>%
    count(metaclone, multiple, name = "n_unit") %>%
    tidyr::complete(metaclone, multiple, fill = list(n_unit = 0)) %>%
    tidyr::pivot_wider(names_from = multiple, values_from = n_unit, values_fill = 0) %>%
    mutate(total = multiple + one)
}

fisher_pval <- function(g1, g2, W) {
  cats <- unique(W$metaclone)
  map_alias <- function(g) {
    if (g %in% cats) return(g)
    if (g == "bothCentric" && "none" %in% cats) return("none")
    if (g == "none" && "bothCentric" %in% cats) return("bothCentric")
    NA_character_
  }
  g1a <- map_alias(g1); g2a <- map_alias(g2)
  if (is.na(g1a) || is.na(g2a)) return(NA_real_)
  r1 <- W %>% filter(metaclone == g1a)
  r2 <- W %>% filter(metaclone == g2a)
  if (!nrow(r1) || !nrow(r2)) return(NA_real_)
  mat <- matrix(c(r1$multiple, r1$one, r2$multiple, r2$one), nrow = 2, byrow = TRUE)
  stats::fisher.test(mat)$p.value
}

p_to_stars <- function(p, show_ns = TRUE, alpha_ns = 0.05) {
  if (is.na(p)) return("n/a")
  if (p < 1e-4) "****" else if (p < 1e-3) "***" else if (p < 1e-2) "**"
  else if (p < alpha_ns) "*" else if (show_ns) "ns" else ""
}
p_to_text <- function(p) if (is.na(p)) "n/a" else if (p < 1e-4) "p<1e-4" else paste0("p = ", format(p, digits = 2))

## ---------- P-values ----------
wide_stats <- to_wide(analysis_stats)
wide_plot  <- to_wide(analysis_ds)
wide <- if (PVAL_FROM_PRE_DOWNSAMPLED) wide_stats else wide_plot

third_label_for_plot <- if (PLOT_THIRD_GROUP %in% c("none","bothCentric")) PLOT_THIRD_GROUP else "none"
all_pairs <- list(
  c("alphaCentric","betaCentric"),
  c("alphaCentric", third_label_for_plot),
  c("betaCentric",  third_label_for_plot)
)
pairs_present <- Filter(function(p) all(p %in% plot_levels_display), all_pairs)

pvals_df <- tibble::tibble(
  group1 = vapply(pairs_present, `[[`, character(1), 1),
  group2 = vapply(pairs_present, `[[`, character(1), 2)
) %>%
  rowwise() %>%
  mutate(
    p     = fisher_pval(group1, group2, wide),
    label = if (PVAL_USE_STARS) p_to_stars(p, show_ns = PVAL_SHOW_NS, alpha_ns = ALPHA_NS) else p_to_text(p)
  ) %>%
  ungroup() %>%
  mutate(
    y.position = Y_START + Y_STEP * (row_number() - 1),
    xmin = factor(group1, levels = plot_levels_display),
    xmax = factor(group2, levels = plot_levels_display)
  )

y_top <- if (nrow(pvals_df)) max(1, max(pvals_df$y.position) + Y_MARGIN) else 1
y_top <- y_top + EXTRA_Y_HEADROOM

## ---------- PLOT ----------
theme_set(theme_bw(BASE_SIZE_PT))
cbPalette <- c('#66c2a5','#fc8d62')  # multiple / one

p <- ggplot(plot_df, aes(x = metaclone, y = n_unit, fill = multiple)) +
  geom_col(position = "fill", width = BAR_WIDTH) +
  scale_x_discrete(limits = plot_levels_display, drop = FALSE) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.25),
                     labels = label_number(accuracy = 0.01),
                     expand = c(0, 0)) +
  scale_fill_manual(values = c(multiple = cbPalette[1], one = cbPalette[2]), drop = FALSE) +
  coord_cartesian(ylim = c(0, y_top), clip = "off") +
  labs(
    x = "Metaclone",
    y = y_label,
    fill = "No. epitopes"
  ) +
  theme(
    panel.grid   = element_blank(),
    legend.key   = element_blank(),
    axis.text.x  = element_text(angle = AXIS_X_ANGLE, hjust = 1),
    axis.title.x = element_text(margin = margin(t = 14)),
    plot.margin  = margin(t = 10, r = 16, b = 36, l = 16),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.3),
    axis.line    = element_line(linewidth = 0.3),
    axis.ticks   = element_line(linewidth = 0.3)
  )

if (nrow(pvals_df)) {
  p <- p + ggpubr::stat_pvalue_manual(
    data = pvals_df,
    label = "label",
    xmin = "xmin", xmax = "xmax",
    y.position = "y.position",
    tip.length = 0,
    step.increase = 0,
    bracket.size  = BRACKET_LW,
    size   = if (PVAL_USE_STARS) STAR_SIZE_PT / ggplot2::.pt else TEXT_SIZE_PT / ggplot2::.pt,
    family = "sans",
    inherit.aes = FALSE
  )
}

print(p)

## ---------- SAVE (PDF ) ----------
# PDF (vector)
# PDF (vector, Cairo)
ggsave(
  filename = OUTFILE_PDF,
  plot     = p,
  width    = PLOT_W_IN,
  height   = PLOT_H_IN,
  units    = "in",
  device   = cairo_pdf,
  bg       = "white"
)

cat("\nSaved:\n", OUTFILE_PDF, "\n", sep = "")

## ---------- Diagnostics ----------
cat("\n=== Distribution of global CDR3 epitope counts (CHAIN_MODE only) ===\n")
print(cdr3_ep_counts %>% count(n_epitopes_cdr3) %>% arrange(desc(n)))

cat("\n=== Classification (pre-downsample) ===\n")
print(analysis_stats %>% count(metaclone, multiple) %>% arrange(metaclone, multiple))

cat("\n=== Classification (post-downsample) ===\n")
print(analysis_ds %>% count(metaclone, multiple) %>% arrange(metaclone, multiple))

by_meta_before <- analysis_stats %>%
  distinct(metaclone = as.character(metaclone), unit_id) %>%
  count(metaclone, name = "n_units_before")

by_meta_after <- analysis_ds %>%
  distinct(metaclone = as.character(metaclone), unit_id) %>%
  count(metaclone, name = "n_units_after")

by_meta_both <- dplyr::full_join(by_meta_before, by_meta_after, by = "metaclone") %>%
  mutate(
    metaclone = factor(metaclone, levels = LVL_ORDER, ordered = TRUE),
    fraction_kept = n_units_after / n_units_before
  ) %>%
  arrange(metaclone)

cat("\n=== Units by metaclone (pre vs post downsample) ===\n")
print(by_meta_both)
