
# ────────────────────────────────────────────────────────────────────────────────
# End-to-end: chain-centric classification + J-region length densities + KS table
# ────────────────────────────────────────────────────────────────────────────────

## ===== End-to-end: MHCI vs MHCII across alpha/beta/both/neither; BH + Fisher-safe =====
rm(list = ls())

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(stringr); library(tibble); library(readr)
  library(ggplot2); library(scales)
})

## --------------------------- TOGGLES ---------------------------------------
DATA_SOURCE       <- "file"   # "file" | "object"
DATA_PATH         <- "~/Desktop/Shubhams_paper_code/data/hu_with_metaclone_20250826_105038.csv"
DATA_OBJECT_NAME  <- "hu"

CHAIN_FOR_GROUP   <- "TRA"    # which chain defines the lib set: "TRA" or "TRB"
FILTER_INKT_MAIT  <- TRUE
REQUIRE_MOTIF     <- TRUE

## Plot-only downsampling (set STATS_MATCH_PLOT=TRUE to run stats on DS too)
DOWNSAMPLE_GROUPS <- TRUE
DOWNSAMPLE_SCOPE  <- "all"        # "alpha_beta" or "all"
DOWNSAMPLE_SEED   <- 2025
STATS_MATCH_PLOT  <- TRUE         # TRUE = stats use downsampled counts

## Fisher/Stats safety knobs (fix FEXACT error 6 for big tables)
EXACT_TEST_MODE   <- "auto"       # "auto" | "exact" | "simulate" | "chi2"
FISHER_WORKSPACE  <- 2e8
FISHER_WORKSPACE_MAX <- 1e9
MC_B              <- 2e5
MC_SEED           <- 123
CHI2_CORRECT      <- FALSE

## Multiple testing for pairwise group contrasts
P_ADJUST_METHOD   <- "BH"

## Display knobs
SHOW_LEGEND           <- TRUE
SHOW_TOTALS_IN_LABELS <- FALSE
DRAW_BRACKET          <- TRUE
STARS_SIZE            <- 12
STARS_Y_BASE          <- 1.02
STARS_Y_STEP          <- 0.06

## Bar order (all four bars)
X_LEVELS <- c("alphaCentric","betaCentric","bothCentric","neitherCentric")

set.seed(42)

## --------------------------- Helpers ---------------------------------------
ensure_column <- function(df, synonyms, as_name) {
  if (as_name %in% names(df)) return(df)
  hit <- intersect(synonyms, names(df)); if (!length(hit)) return(df)
  dplyr::rename(df, !!as_name := dplyr::all_of(hit[1]))
}
strip_allele <- function(s) {
  s <- if (is.factor(s)) as.character(s) else s
  toupper(trimws(sub("^([A-Z0-9-]+).*", "\\1", s)))
}
read_robust <- function(path) {
  suppressWarnings({
    a <- try(readr::read_csv(path, show_col_types = FALSE, progress = FALSE), silent = TRUE)
    b <- try(readr::read_tsv(path, show_col_types = FALSE, progress = FALSE), silent = TRUE)
  })
  if (inherits(a, "data.frame") && ncol(a) > 1) return(a)
  if (inherits(b, "data.frame") && ncol(b) > 1) return(b)
  stop("Could not read file as CSV/TSV (or only one column).")
}
norm_class <- function(x) {
  s <- toupper(trimws(as.character(x)))
  is_i  <- grepl("\\b(MHC\\s*)?CLASS\\s*I\\b|\\bMHC\\s*I\\b|^I$",  s)
  is_ii <- grepl("\\b(MHC\\s*)?CLASS\\s*II\\b|\\bMHC\\s*II\\b|^II$", s)
  dplyr::case_when(
    is_i  & !is_ii ~ "I",
    is_ii & !is_i  ~ "II",
    is_i  &  is_ii ~ "both",
    TRUE           ~ "unknown"
  )
}
ensure_bool_cols <- function(df, cols = c("I","II","both")) {
  for (nm in cols) if (!nm %in% names(df)) df[[nm]] <- FALSE
  df
}
## Vector-safe stars for (adjusted) p-values
p_to_stars <- function(p) {
  p <- as.numeric(p)
  out <- rep("ns", length(p))
  out[is.na(p)] <- "n/a"
  out[p < 0.05]  <- "*"
  out[p < 0.01]  <- "**"
  out[p < 1e-3]  <- "***"
  out[p < 1e-4]  <- "****"
  out
}
## Fisher wrappers (RxC-safe with fallbacks)
fisher_safe <- function(M,
                        mode = EXACT_TEST_MODE,
                        workspace = FISHER_WORKSPACE,
                        workspace_max = FISHER_WORKSPACE_MAX,
                        B = MC_B, seed = MC_SEED, chi2_correct = CHI2_CORRECT) {
  M <- as.matrix(M)
  if (mode == "chi2") {
    cs <- suppressWarnings(chisq.test(M, correct = chi2_correct))
    return(list(p.value = cs$p.value, method = "Chi-square (chisq.test)"))
  }
  if (mode == "simulate") {
    set.seed(seed); fs <- suppressWarnings(fisher.test(M, simulate.p.value = TRUE, B = B))
    return(list(p.value = fs$p.value, method = sprintf("Fisher (Monte Carlo, B=%d)", B)))
  }
  if (mode == "exact") {
    fs <- try(suppressWarnings(fisher.test(M, hybrid = TRUE, workspace = workspace)), silent = TRUE)
    if (!inherits(fs, "try-error")) return(list(p.value = fs$p.value,
                                                method = sprintf("Fisher exact (workspace=%.0f, hybrid=TRUE)", workspace)))
    stop(attr(fs, "condition")$message)
  }
  ## auto
  ws <- workspace
  repeat {
    fs <- try(suppressWarnings(fisher.test(M, hybrid = TRUE, workspace = ws)), silent = TRUE)
    if (!inherits(fs, "try-error")) {
      return(list(p.value = fs$p.value,
                  method = sprintf("Fisher exact (workspace=%.0f, hybrid=TRUE)", ws)))
    }
    if (ws >= workspace_max) break
    ws <- min(workspace_max, ceiling(ws * 2))
  }
  set.seed(seed)
  fs2 <- try(suppressWarnings(fisher.test(M, simulate.p.value = TRUE, B = B)), silent = TRUE)
  if (!inherits(fs2, "try-error")) {
    return(list(p.value = fs2$p.value, method = sprintf("Fisher (Monte Carlo, B=%d)", B)))
  }
  cs <- suppressWarnings(chisq.test(M, correct = chi2_correct))
  list(p.value = cs$p.value, method = "Chi-square (fallback)")
}

## ---------------------- Load & standardize ---------------------------------
hu <- if (DATA_SOURCE == "file") read_robust(DATA_PATH) else {
  if (!exists(DATA_OBJECT_NAME)) stop(sprintf("Object '%s' not found.", DATA_OBJECT_NAME))
  get(DATA_OBJECT_NAME)
}
hu0_n <- nrow(hu)

hu <- hu %>%
  ensure_column(c("junction","cdr3","cdr3aa","junction_aa","aa.cdr3","AACDR3"), "junction") %>%
  ensure_column(c("libid","complex.id","pair_id","pair.id","cell_id"),         "libid") %>%
  ensure_column(c("v.segm","v_gene","V.GENE.and.allele","v_gene_raw","V.gene","VGene"), "v_gene_raw") %>%
  ensure_column(c("j.segm","j_gene","J.GENE.and.allele","j_gene_raw","J.gene","JGene"), "j_gene_raw") %>%
  ensure_column(c("MHC.class","mhc.class","MHC_class","mhc_class","MHCClass","mhcClass","MHC.CLASS"), "mhc_class")

req <- c("junction","libid","v_gene_raw","j_gene_raw","mhc_class")
missing <- setdiff(req, names(hu)); if (length(missing)) stop("Missing columns: ", paste(missing, collapse=", "))

if ("species" %in% names(hu)) {
  hu <- hu %>%
    mutate(species_norm = str_replace_all(.data$species, "\\s+", "")) %>%
    filter(is.na(species_norm) | species_norm == "HomoSapiens") %>%
    select(-species_norm)
}

hu <- hu %>%
  filter(!is.na(junction), junction != "",
         !is.na(v_gene_raw), v_gene_raw != "",
         !is.na(j_gene_raw), j_gene_raw != "")
if (REQUIRE_MOTIF) hu <- hu %>% filter(str_detect(junction, "^[Cc].*[FfWw]$"))
hu <- hu %>% mutate(v_gene = strip_allele(v_gene_raw),
                    j_gene = strip_allele(j_gene_raw))

chaincol <- intersect(c("chain","gene","Gene","tcr.chain"), names(hu))
hu <- hu %>%
  mutate(
    chain_guess = if (length(chaincol)) toupper(.data[[chaincol[1]]]) else NA_character_,
    chain = ifelse(is.na(chain_guess) | chain_guess == "",
                   ifelse(grepl("^TRAV", v_gene) | grepl("^TRAJ", j_gene), "TRA",
                   ifelse(grepl("^TRBV", v_gene) | grepl("^TRBJ", j_gene), "TRB", NA_character_)),
                   chain_guess)
  ) %>%
  filter(chain %in% c("TRA","TRB")) %>%
  select(-chain_guess)

if (isTRUE(FILTER_INKT_MAIT)) {
  hu <- hu %>%
    filter(!(junction == "CVVSDRGSTLGRLYF")) %>%
    filter(!(v_gene == "TRAV1-2" & j_gene %in% c("TRAJ33","TRAJ20","TRAJ12")))
}

## --------------- Paired libs & metaclone groups (lib-level) ----------------
paired_ids <- hu %>% distinct(libid, chain) %>% count(libid, name="nchain") %>%
  filter(nchain >= 2) %>% pull(libid)
hu <- hu %>% filter(libid %in% paired_ids)

tra_all <- hu %>% filter(chain == "TRA")
trb_all <- hu %>% filter(chain == "TRB")

paired_all <- inner_join(
  tra_all %>% select(libid, junction.tra = junction) %>% distinct(),
  trb_all %>% select(libid, junction.trb = junction) %>% distinct(),
  by = "libid"
)

alpha_ids <- paired_all %>% group_by(junction.tra) %>%
  filter(dplyr::n_distinct(junction.trb) > 1) %>% pull(libid) %>% unique()
beta_ids  <- paired_all %>% group_by(junction.trb) %>%
  filter(dplyr::n_distinct(junction.tra) > 1) %>% pull(libid) %>% unique()
both_ids     <- intersect(alpha_ids, beta_ids)
neither_ids  <- setdiff(paired_all$libid, union(alpha_ids, beta_ids))

chain_df <- if (toupper(CHAIN_FOR_GROUP) == "TRA") tra_all else trb_all

## Assign ONE of four groups to each lib (priority: both > alpha > beta > neither)
groups_stats <- chain_df %>%
  distinct(libid) %>%
  mutate(group = case_when(
    libid %in% both_ids    ~ "bothCentric",
    libid %in% alpha_ids   ~ "alphaCentric",
    libid %in% beta_ids    ~ "betaCentric",
    libid %in% neither_ids ~ "neitherCentric",
    TRUE                   ~ NA_character_
  )) %>%
  filter(!is.na(group)) %>%
  mutate(group = factor(group, levels = X_LEVELS))
stopifnot(nrow(groups_stats) > 0)

## --------------- Per-lib restriction (I/II/both/unknown) -------------------
lib_cls_wide <- hu %>%
  transmute(libid, cls = norm_class(mhc_class)) %>%
  distinct(libid, cls) %>%
  mutate(seen = TRUE) %>%
  tidyr::pivot_wider(names_from = cls, values_from = seen, values_fill = FALSE)
lib_cls_wide <- ensure_bool_cols(lib_cls_wide, c("I","II","both"))

lib_labels <- lib_cls_wide %>%
  mutate(
    has_I  = I,
    has_II = II,
    has_both_label = both,
    mhc_cat = case_when(
      (has_I & has_II) | has_both_label ~ "both",
      has_I  ~ "I",
      has_II ~ "II",
      TRUE   ~ "unknown"
    )
  ) %>% select(libid, mhc_cat)

## ---------------------- Optional downsampling (PLOT ONLY) ------------------
groups_plot <- groups_stats
if (isTRUE(DOWNSAMPLE_GROUPS)) {
  set.seed(DOWNSAMPLE_SEED)
  scope <- if (DOWNSAMPLE_SCOPE == "alpha_beta")
             intersect(levels(groups_plot$group), c("alphaCentric","betaCentric"))
           else levels(groups_plot$group)
  scope <- scope[scope %in% unique(groups_plot$group)]
  if (length(scope) > 1) {
    size_tbl <- groups_plot %>% filter(group %in% scope) %>% count(group, name = "n")
    target_n <- min(size_tbl$n, na.rm = TRUE)
    sampled <- groups_plot %>%
      filter(group %in% scope) %>%
      group_by(group) %>% group_split()
    sampled <- lapply(sampled, function(df) df %>% slice_sample(n = min(target_n, nrow(df))))
    kept_scoped <- bind_rows(sampled)
    others <- groups_plot %>% filter(!group %in% scope)
    groups_plot <- bind_rows(kept_scoped, others) %>%
      mutate(group = factor(group, levels = X_LEVELS))
  }
}
counts_source <- if (isTRUE(STATS_MATCH_PLOT)) groups_plot else groups_stats

## ---------------- Build table used for BOTH stats & plot -------------------
tab_raw <- counts_source %>%
  left_join(lib_labels, by = "libid") %>%
  filter(group %in% X_LEVELS, mhc_cat %in% c("I","II")) %>%
  mutate(
    mhc_cat = recode(mhc_cat, "I" = "MHCI", "II" = "MHCII"),
    group   = factor(group, levels = X_LEVELS)
  ) %>%
  count(group, mhc_cat, name = "n") %>%
  tidyr::complete(
    group   = factor(X_LEVELS, levels = X_LEVELS),
    mhc_cat = factor(c("MHCI","MHCII"), levels = c("MHCI","MHCII")),
    fill = list(n = 0)
  ) %>%
  group_by(group) %>%
  mutate(tot = sum(n), frac = ifelse(tot > 0, n/tot, 0),
         pct = percent(frac, accuracy = 0.1),
         lab = if (isTRUE(SHOW_TOTALS_IN_LABELS)) sprintf("%s\n(%d/%d)", pct, n, tot) else pct) %>%
  ungroup()

cat("\n[Table used for BOTH stats and plot]\n")
print(tab_raw %>% select(group, mhc_cat, n, tot, pct), n = Inf)

## Save raw tables for provenance
write_csv(tab_raw, "mhc_counts_alpha_beta_both_neither.csv")

## ---------------- Statistics (Fisher-safe; BH pairwise) --------------------
mat_global <- tab_raw %>%
  select(group, mhc_cat, n) %>%
  tidyr::pivot_wider(names_from = mhc_cat, values_from = n) %>%
  arrange(group)
M <- as.matrix(mat_global[, c("MHCI","MHCII")]); rownames(M) <- mat_global$group

gl <- fisher_safe(M)
cat(sprintf("\nGlobal test method: %s\n", gl$method))
cat(sprintf("Global p-value: %s\n", format.pval(gl$p.value, digits = 3, eps = 1e-4)))

pairs <- if (length(levels(tab_raw$group)) >= 2) combn(levels(tab_raw$group), 2, simplify = FALSE) else list()
pairwise <- tibble(
  g1 = if (length(pairs)) vapply(pairs, `[`, character(1), 1) else character(0),
  g2 = if (length(pairs)) vapply(pairs, `[`, character(1), 2) else character(0)
)
pair_p <- function(g1, g2) {
  if (!all(c(g1,g2) %in% rownames(M))) return(list(p=NA_real_, method="n/a"))
  mm <- M[c(g1,g2), , drop = FALSE]
  fr <- fisher_safe(mm)
  list(p = fr$p.value, method = fr$method)
}
if (nrow(pairwise)) {
  pw <- lapply(seq_len(nrow(pairwise)), function(i) pair_p(pairwise$g1[i], pairwise$g2[i]))
  pairwise$p <- vapply(pw, function(x) x$p, numeric(1))
  pairwise$method <- vapply(pw, function(x) x$method, character(1))
  pairwise <- pairwise %>%
    mutate(p_adj = p.adjust(p, method = P_ADJUST_METHOD),
           stars = p_to_stars(p_adj))
  write_csv(pairwise, "mhc_pairwise_pvalues_BH.csv")
  cat("\nPairwise tests (BH-adjusted):\n")
  print(pairwise %>% mutate(p_adj = format.pval(p_adj, digits = 3, eps = 1e-4)), n = Inf)
} else {
  cat("\n[No pairwise tests: fewer than 2 groups present]\n")
}

## ---------------- Plot ------------------------------------------------------
if (dev.cur() > 1) dev.off()
if (Sys.info()[["sysname"]] == "Darwin")
  quartz(height = 14, width = 9 + 1*(length(levels(tab_raw$group))-2), dpi = 72)

theme_set(
  theme_bw(36) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.key       = element_blank(),
      axis.text.x      = element_text(angle = 45, hjust = 1),
      plot.margin      = margin(t = 24, r = 32, b = 36, l = 32),
      panel.border     = element_rect(linewidth = 0.3, colour = "black"),
      axis.ticks       = element_line(linewidth = 0.3),
      axis.line        = element_line(linewidth = 0.3)
    )
)

pal <- c(MHCI = "#66c2a5", MHCII = "#fc8d62")

gg <- ggplot(tab_raw, aes(x = group, y = n, fill = mhc_cat)) +
  geom_col(position = "fill", width = 0.7) +
  geom_text(aes(label = lab),
            position = position_fill(vjust = 0.5),
            lineheight = 0.95, size = 8, color = "black") +
  scale_fill_manual(values = pal, name = "Restriction",
                    guide = if (SHOW_LEGEND) "legend" else "none") +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(x = NULL, y = "Fraction of TCRs\n")

add_bracket <- function(g, x1, x2, y, lbl) {
  if (!DRAW_BRACKET) return(g + annotate("text", x = mean(c(x1,x2)), y = y, label = lbl,
                                         size = STARS_SIZE, vjust = 0, fontface = "bold"))
  g +
    annotate("segment", x = x1, xend = x2, y = y, yend = y, linewidth = 0.3) +
    annotate("segment", x = x1, xend = x1, y = y, yend = y - 0.02, linewidth = 0.3) +
    annotate("segment", x = x2, xend = x2, y = y, yend = y - 0.02, linewidth = 0.3) +
    annotate("text", x = mean(c(x1,x2)), y = y, label = lbl, size = STARS_SIZE,
             vjust = 0, fontface = "bold")
}
if (exists("pairwise") && nrow(pairwise)) {
  gg <- gg + coord_cartesian(ylim = c(0, STARS_Y_BASE + (nrow(pairwise)+1)*STARS_Y_STEP + 0.04),
                             clip = "off")
  xpos <- setNames(seq_along(levels(tab_raw$group)), levels(tab_raw$group))
  y0 <- STARS_Y_BASE
  for (i in seq_len(nrow(pairwise))) {
    g1 <- pairwise$g1[i]; g2 <- pairwise$g2[i]
    if (is.na(pairwise$p_adj[i])) next
    gg <- add_bracket(gg, xpos[g1], xpos[g2], y0 + (i-1)*STARS_Y_STEP, pairwise$stars[i])
  }
}
print(gg)

## ---------------- Diagnostics / provenance ---------------------------------
cat("\n[Row counts]\n")
cat(sprintf("hu (input): %d rows\n", hu0_n))
cat(sprintf("hu (after filters, paired libs): %d rows\n", nrow(hu)))
cat(sprintf("Unique libs in groups_stats: %d\n", n_distinct(groups_stats$libid)))
cat(sprintf("Unique libs in counts_source: %d\n", n_distinct(counts_source$libid)))
