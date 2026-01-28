## =========================
## End-to-end enrichment (V/J or MHCI/MHCII) + PDF output
## MAKE RESULTS MATCH "FIRST SCRIPT" ON CD4merge07212025.rds:
## - If hu contains columns study + Status, subset to: study=="su" & Status=="Healthy"
##   (this reproduces: plot_base <- CD4merge %>% filter(study=="su", Status=="Healthy"))
## - Everything else keeps your toggles + file locations
## =========================

rm(list = ls())

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(tibble)
  library(ggplot2)
  library(stringr)
  library(rlang)
  library(readr)
  library(scales)
})

set.seed(1)

## ---------- TOGGLES ----------
CHAIN_TARGET           <- "TRA"      # "TRA" or "TRB" (chain used for enrichment rows)
FEATURE_MODE           <- "mhcII"        # "v" | "j" | "mhcI" | "mhcII"
KEEP_ALLELES           <- FALSE      # TRUE = keep allele strings; FALSE = strip alleles
DATA_SOURCE            <- "file"     # "file" | "object"
#DATA_PATH              <- "~/Desktop/Shubhams_paper_code/data/CD4merge07212025.rds"  # .rds or .csv/.tsv
DATA_PATH         <- "~/Desktop/Shubhams_paper_code/data/hu_with_metaclone_20250826_105038.csv" # use this for MHC

DATA_OBJECT_NAME       <- "myTCRs"   # used if DATA_SOURCE == "object"
FILTER_INKT_MAIT       <- TRUE

DOWNSAMPLE_JUNCTIONS   <- FALSE
SEED_DOWNSAMPLE        <- 2025

RANDOMIZE_SAMPLES      <- FALSE
SEED_RANDOM            <- 42

MIN_GROUP_COUNT        <- 10
REQUIRE_BOTH_GROUPS    <- FALSE

PLOT_KEEP_NS_TOP       <- 5

## ---------- OUTPUT (timestamped) ----------
OUTPUT_DIR <- "~/Desktop/Shubhams_paper_code/Figure_PDFs"
RUN_TIMESTAMP <- format(Sys.time(), "%Y%m%d_%H%M%S")

OUTPUT_TAG <- paste(
  CHAIN_TARGET,
  FEATURE_MODE,
  if (KEEP_ALLELES) "allele" else "gene",
  RUN_TIMESTAMP,
  sep = "_"
)

if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)

out_pdf <- function(stem) file.path(OUTPUT_DIR, paste0(stem, ".pdf"))
out_csv <- function(stem) file.path(OUTPUT_DIR, paste0(stem, ".csv"))
out_rds <- function(stem) file.path(OUTPUT_DIR, paste0(stem, ".rds"))
out_txt <- function(stem) file.path(OUTPUT_DIR, paste0(stem, ".txt"))

## ---------- HELPERS ----------
ensure_column <- function(df, synonyms, as_name) {
  if (as_name %in% names(df)) return(df)
  hit <- intersect(synonyms, names(df))
  if (!length(hit)) return(df)
  dplyr::rename(df, !!as_name := dplyr::all_of(hit[1]))
}

strip_allele <- function(s) {
  s <- if (is.factor(s)) as.character(s) else s
  s <- trimws(toupper(s))
  sub("^([A-Z0-9-]+).*", "\\1", s, perl = TRUE)
}

read_robust_table <- function(path) {
  suppressWarnings({
    csv_try <- try(readr::read_csv(path, show_col_types = FALSE, progress = FALSE), silent = TRUE)
    tsv_try <- try(readr::read_tsv(path, show_col_types = FALSE, progress = FALSE), silent = TRUE)
  })
  cand <- list(csv = csv_try, tsv = tsv_try)
  cand <- cand[vapply(cand, function(x) inherits(x, "data.frame"), logical(1))]
  if (!length(cand)) stop("Failed to read file as CSV or TSV: ", path)
  pick <- cand[[which.max(vapply(cand, ncol, integer(1)))]]
  if (ncol(pick) == 1L) stop("File parsed into a single column. Likely wrong delimiter.")
  pick
}

read_any <- function(path) {
  ext <- tolower(tools::file_ext(path))
  if (ext == "rds") {
    obj <- readRDS(path)
    if (!inherits(obj, "data.frame")) stop("RDS did not contain a data.frame: ", path)
    return(obj)
  }
  if (ext %in% c("csv","tsv","txt")) return(read_robust_table(path))
  read_robust_table(path)
}

split_mhc_any <- function(x) {
  if (is.null(x)) return(character(0))
  v <- unlist(strsplit(as.character(x), "[,;+/\\|& ]+"))
  v <- toupper(trimws(v))
  v <- gsub("\u00A0", "", v)
  v <- v[nzchar(v)]
  unique(v)
}

guess_mhc_class <- function(tok) {
  t <- toupper(trimws(tok))
  t <- sub("^HLA[-:]", "", t)
  if (grepl("^(A|B|C|E|F|G)\\b", t)) return("I")
  if (grepl("^(DRB[1-5]|DQA1|DQB1|DPA1|DPB1|DRA)\\b", t)) return("II")
  "unknown"
}

find_mhc_columns <- function(df) {
  nm <- names(df)
  hits <- nm[grepl("^(mhc|hla)", nm, ignore.case = TRUE)]
  setdiff(hits, nm[grepl("class", nm, ignore.case = TRUE)])
}

run_fisher_table <- function(a, b, c, d) {
  any_zero <- any(c(a,b,c,d) == 0)
  a0 <- ifelse(any_zero, a + 0.5, a)
  b0 <- ifelse(any_zero, b + 0.5, b)
  c0 <- ifelse(any_zero, c + 0.5, c)
  d0 <- ifelse(any_zero, d + 0.5, d)
  list(
    p  = fisher.test(matrix(c(a0,b0,c0,d0), nrow = 2, byrow = TRUE))$p.value,
    OR = (a0 / b0) / (c0 / d0)
  )
}

## ---------- LOAD DATA ----------
if (DATA_SOURCE == "file") {
  hu <- read_any(DATA_PATH)
} else if (DATA_SOURCE == "object") {
  if (!exists(DATA_OBJECT_NAME)) stop(sprintf("Object '%s' not found.", DATA_OBJECT_NAME))
  hu <- get(DATA_OBJECT_NAME)
  if (!inherits(hu, "data.frame")) stop("DATA_OBJECT_NAME did not resolve to a data.frame.")
} else stop('DATA_SOURCE must be "file" or "object".')

## ---------- CRITICAL: MATCH FIRST SCRIPT'S OBJECT PREP ----------
## First script does:
## plot_base <- CD4merge %>% filter(study == "su", Status == "Healthy")
## myTCRs <- plot_base
## We replicate that *here* if those columns exist, without changing any toggles/paths.
if (all(c("study", "Status") %in% names(hu))) {
  hu <- hu %>% filter(.data$study == "su", .data$Status == "Healthy")
  message("Applied cohort filter to match first script: study=='su' & Status=='Healthy'")
} else {
  message("No (study, Status) columns found; skipping cohort filter.")
}

## ---------- STANDARDIZE REQUIRED COLUMNS ----------
hu <- hu %>%
  ensure_column(c("junction","cdr3","aa.cdr3","AACDR3","cdr3aa","junction_aa"), "junction") %>%
  ensure_column(c("libid","complex.id","pair.id","pair_id","cell_id","barcode"), "libid") %>%
  ensure_column(c("v_gene_raw","v_gene","v.segm","V.GENE.and.allele","V.gene"), "v_gene_raw") %>%
  ensure_column(c("j_gene_raw","j_gene","j.segm","J.GENE.and.allele","J.gene"), "j_gene_raw")

if (!"libid" %in% names(hu)) stop("Missing 'libid' (or synonym).")
if (!"junction" %in% names(hu)) stop("Missing 'junction' (CDR3 AA).")
if (!all(c("v_gene_raw","j_gene_raw") %in% names(hu))) stop("Missing V/J genes (need 'v_gene_raw' and 'j_gene_raw').")

hu <- hu %>%
  filter(!is.na(junction), junction != "",
         !is.na(.data$v_gene_raw), .data$v_gene_raw != "",
         !is.na(.data$j_gene_raw), .data$j_gene_raw != "",
         str_detect(junction, "^[Cc].*[FfWw]$")) %>%
  mutate(
    v_gene = if (KEEP_ALLELES) toupper(trimws(.data$v_gene_raw)) else strip_allele(.data$v_gene_raw),
    j_gene = if (KEEP_ALLELES) toupper(trimws(.data$j_gene_raw)) else strip_allele(.data$j_gene_raw)
  )

chaincol <- intersect(c("chain","gene","Gene","tcr.chain","chainType"), names(hu))
hu <- hu %>%
  mutate(
    chain_guess = if (length(chaincol)) toupper(as.character(.data[[chaincol[1]]])) else NA_character_,
    chain = ifelse(
      is.na(chain_guess) | chain_guess == "",
      ifelse(grepl("^TRAV", v_gene) | grepl("^TRAJ", j_gene), "TRA",
             ifelse(grepl("^TRBV", v_gene) | grepl("^TRBJ", j_gene), "TRB", NA_character_)),
      chain_guess
    )
  ) %>%
  filter(chain %in% c("TRA","TRB")) %>%
  select(-chain_guess)

## ---------- REMOVE iNKT/MAIT (optional) ----------
if (isTRUE(FILTER_INKT_MAIT)) {
  iNkt_libs <- hu %>% filter(junction == "CVVSDRGSTLGRLYF") %>% pull(libid)
  mait_libs <- hu %>% filter(v_gene == "TRAV1-2", j_gene %in% c("TRAJ33","TRAJ20","TRAJ12")) %>% pull(libid)
  hu <- hu %>% filter(!(libid %in% union(iNkt_libs, mait_libs)))
}

## ---------- KEEP PAIRED LIBRARIES (must have both TRA & TRB) ----------
paired_ids <- hu %>%
  distinct(libid, chain) %>%
  count(libid, name = "n_chain") %>%
  filter(n_chain >= 2) %>%
  pull(libid)
hu <- hu %>% filter(libid %in% paired_ids)

## ---------- BUILD METACLONES FROM PAIRED TRA/TRB ----------
tra_all <- hu %>% filter(chain == "TRA")
trb_all <- hu %>% filter(chain == "TRB")

paired_all <- inner_join(
  tra_all %>% select(libid, junction.tra = junction, v.tra = v_gene, j.tra = j_gene),
  trb_all %>% select(libid, junction.trb = junction, v.trb = v_gene, j.trb = j_gene),
  by = "libid"
)

alpha_meta <- paired_all %>% group_by(junction.tra) %>% filter(n_distinct(junction.trb) > 1) %>% ungroup()
beta_meta  <- paired_all %>% group_by(junction.trb) %>% filter(n_distinct(junction.tra) > 1) %>% ungroup()

alpha_ids <- unique(alpha_meta$libid)
beta_ids  <- unique(beta_meta$libid)
overlap_ids <- intersect(alpha_ids, beta_ids)
if (length(overlap_ids)) message("Dropping ", length(overlap_ids), " libids both alpha- and beta-centric.")
alpha_only_ids <- setdiff(alpha_ids, overlap_ids)
beta_only_ids  <- setdiff(beta_ids,  overlap_ids)

## ---------- SELECT CHAIN TO TEST; TAG GROUPS ----------
chain_df <- if (toupper(CHAIN_TARGET) == "TRA") tra_all else trb_all
hu_chain <- chain_df %>%
  filter(libid %in% c(alpha_only_ids, beta_only_ids)) %>%
  mutate(group = ifelse(libid %in% alpha_only_ids, "alphaCentric", "betaCentric"))

## ---------- Optional randomization, then down-sampling by unique junctions ----------
if (isTRUE(RANDOMIZE_SAMPLES)) {
  set.seed(SEED_RANDOM)
  lib_groups <- hu_chain %>% distinct(libid, group)
  n_alpha <- sum(lib_groups$group == "alphaCentric")
  n_beta  <- sum(lib_groups$group == "betaCentric")
  libs <- sample(lib_groups$libid)
  remap <- tibble(
    libid = libs,
    group = c(rep("alphaCentric", n_alpha), rep("betaCentric", n_beta))
  )
  hu_chain <- hu_chain %>% select(-group) %>% left_join(remap, by = "libid")
}

if (isTRUE(DOWNSAMPLE_JUNCTIONS)) {
  set.seed(SEED_DOWNSAMPLE)
  jcounts <- hu_chain %>% distinct(group, junction) %>% count(group, name = "n_j")
  target_n <- min(jcounts$n_j, na.rm = TRUE)

  hu_chain <- hu_chain %>%
    group_by(group) %>%
    group_modify(function(df, key) {
      jset <- unique(df$junction)
      keep_j <- if (length(jset) > target_n) sample(jset, target_n) else jset
      semi_join(df, tibble(junction = keep_j), by = "junction")
    }) %>%
    ungroup()
}

alpha <- hu_chain %>% filter(group == "alphaCentric")
beta  <- hu_chain %>% filter(group == "betaCentric")
if (nrow(alpha) == 0 || nrow(beta) == 0) stop("Empty alpha or beta group after filtering.")

## ---------- ENRICHMENT ----------
if (tolower(FEATURE_MODE) %in% c("v","j")) {

  feature_col   <- if (tolower(FEATURE_MODE) == "v") "v_gene" else "j_gene"
  feature_label <- paste0(CHAIN_TARGET, " ", toupper(FEATURE_MODE), " ", if (KEEP_ALLELES) "allele" else "gene")

  alpha_tot <- nrow(alpha)
  beta_tot  <- nrow(beta)

  alpha_counts <- alpha %>% count(.data[[feature_col]], name = "alpha_n")
  beta_counts  <- beta  %>% count(.data[[feature_col]],  name = "beta_n")

  DF_pre <- full_join(beta_counts, alpha_counts, by = feature_col) %>%
    rename(Feature = !!sym(feature_col)) %>%
    mutate(
      alpha_n = coalesce(alpha_n, 0L),
      beta_n  = coalesce(beta_n,  0L),
      a = alpha_n,
      b = alpha_tot - alpha_n,
      c = beta_n,
      d = beta_tot  - beta_n,
      alpha_prop = a / (a + b),
      beta_prop  = c / (c + d)
    )

  if (!is.null(MIN_GROUP_COUNT) && MIN_GROUP_COUNT > 0) {
    DF_pre <- if (REQUIRE_BOTH_GROUPS) {
      DF_pre %>% filter(alpha_n >= MIN_GROUP_COUNT & beta_n >= MIN_GROUP_COUNT)
    } else {
      DF_pre %>% filter(alpha_n >= MIN_GROUP_COUNT | beta_n >= MIN_GROUP_COUNT)
    }
    if (nrow(DF_pre) == 0) stop("All features removed by MIN_GROUP_COUNT in V/J branch.")
  }

  eps <- 1e-12
  DF <- DF_pre %>%
    rowwise() %>%
    mutate(
      res  = list(run_fisher_table(a,b,c,d)),
      pVal = res$p,
      OR   = res$OR,
      dir_raw = ifelse(alpha_prop > beta_prop + eps, "Up in alphaCentric",
                       ifelse(alpha_prop < beta_prop - eps, "Up in betaCentric", "ns_base"))
    ) %>%
    ungroup() %>%
    select(-res) %>%
    mutate(
      fdr = p.adjust(pVal, method = "BH"),
      direction = case_when(
        fdr < 0.05 & dir_raw == "Up in alphaCentric" ~ "Up in alphaCentric",
        fdr < 0.05 & dir_raw == "Up in betaCentric"  ~ "Up in betaCentric",
        TRUE ~ "NS"
      ),
      logFDR = -log10(fdr)
    ) %>%
    arrange(fdr, desc(abs(logFDR)))

} else if (tolower(FEATURE_MODE) %in% c("mhci","mhcii")) {

  mhc_cols_all <- find_mhc_columns(hu)
  if (!length(mhc_cols_all)) stop("No MHC/HLA-like columns found (names starting with 'mhc' or 'HLA').")

  extract_alleles <- function(df_rows) {
    libids <- unique(df_rows$libid)
    base   <- hu %>% filter(libid %in% libids)
    long <- purrr::map_dfr(mhc_cols_all, function(col) {
      tibble(libid = base$libid, value = base[[col]]) %>%
        mutate(tokens = purrr::map(value, split_mhc_any)) %>%
        select(libid, tokens) %>%
        tidyr::unnest(tokens, keep_empty = TRUE) %>%
        transmute(libid, allele = na_if(tokens, ""))
    }) %>%
      distinct() %>%
      filter(!is.na(allele)) %>%
      mutate(mhc_class = vapply(allele, guess_mhc_class, character(1)))
    long
  }

  alpha_alleles <- extract_alleles(alpha) %>% mutate(group = "alphaCentric")
  beta_alleles  <- extract_alleles(beta)  %>% mutate(group = "betaCentric")
  all_alleles   <- bind_rows(alpha_alleles, beta_alleles)

  target_class <- if (tolower(FEATURE_MODE) == "mhci") "I" else "II"
  alleles_test <- all_alleles %>%
    filter(mhc_class == target_class) %>%
    distinct(libid, allele, group)

  alpha_tot_libs <- n_distinct(all_alleles$libid[all_alleles$group == "alphaCentric"])
  beta_tot_libs  <- n_distinct(all_alleles$libid[all_alleles$group == "betaCentric"])

  alpha_counts <- alleles_test %>% filter(group == "alphaCentric") %>% count(allele, name = "alpha_n")
  beta_counts  <- alleles_test %>% filter(group == "betaCentric")  %>% count(allele, name = "beta_n")

  DF_pre <- full_join(beta_counts, alpha_counts, by = "allele") %>%
    rename(Feature = allele) %>%
    mutate(
      alpha_n = coalesce(alpha_n, 0L),
      beta_n  = coalesce(beta_n,  0L),
      a = alpha_n,
      b = alpha_tot_libs - alpha_n,
      c = beta_n,
      d = beta_tot_libs  - beta_n,
      alpha_prop = a / (a + b),
      beta_prop  = c / (c + d)
    )

  if (!is.null(MIN_GROUP_COUNT) && MIN_GROUP_COUNT > 0) {
    DF_pre <- if (REQUIRE_BOTH_GROUPS) {
      DF_pre %>% filter(alpha_n >= MIN_GROUP_COUNT & beta_n >= MIN_GROUP_COUNT)
    } else {
      DF_pre %>% filter(alpha_n >= MIN_GROUP_COUNT | beta_n >= MIN_GROUP_COUNT)
    }
    if (nrow(DF_pre) == 0) stop("All alleles removed by MIN_GROUP_COUNT in MHC branch.")
  }

  eps <- 1e-12
  DF <- DF_pre %>%
    rowwise() %>%
    mutate(
      res  = list(run_fisher_table(a,b,c,d)),
      pVal = res$p,
      OR   = res$OR,
      dir_raw = ifelse(alpha_prop > beta_prop + eps, "Up in alphaCentric",
                       ifelse(alpha_prop < beta_prop - eps, "Up in betaCentric", "ns_base"))
    ) %>%
    ungroup() %>%
    select(-res) %>%
    mutate(
      fdr = p.adjust(pVal, method = "BH"),
      direction = case_when(
        fdr < 0.05 & dir_raw == "Up in alphaCentric" ~ "Up in alphaCentric",
        fdr < 0.05 & dir_raw == "Up in betaCentric"  ~ "Up in betaCentric",
        TRUE ~ "NS"
      ),
      logFDR = -log10(fdr)
    ) %>%
    arrange(fdr, desc(abs(logFDR)))

  feature_label <- paste0("MHC Class ", target_class, " allele")

} else {
  stop('FEATURE_MODE must be one of: "v","j","mhcI","mhcII"')
}

## ---------- PLOT ----------
DF_plot <- DF %>%
  group_by(direction) %>%
  arrange(desc(logFDR), .by_group = TRUE) %>%
  filter(direction != "NS" | dplyr::row_number() <= PLOT_KEEP_NS_TOP) %>%
  ungroup()

pal <- c(
  "Up in alphaCentric" = "#66c2a5",
  "Up in betaCentric"  = "#fc8d62",
  "NS"                 = "gray70"
)

theme_set(
  theme_bw(28) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.key = element_blank(),
      axis.title.x = element_text(margin = margin(t = 14)),
      axis.title.y = element_text(margin = margin(r = 16)),
      plot.margin = margin(t = 14, r = 20, b = 20, l = 28)
    )
)

gg <- ggplot(DF_plot, aes(x = reorder(Feature, logFDR), y = logFDR, fill = direction)) +
  geom_col(width = 0.8) +
  coord_flip(clip = "off") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  scale_fill_manual(values = pal, drop = FALSE, name = "Direction") +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.06))) +
  scale_x_discrete(guide = guide_axis(n.dodge = 1)) +
  labs(
    x = paste0(feature_label, "\n"),
    y = expression("\n-log"[10]*"(FDR)")
  ) +
  theme(plot.title = element_text(size = 16, face = "plain", hjust = 0.5))

print(gg)

## ---------- SAVE ARTIFACTS (PDF) ----------
plot_stem    <- paste0("enrichment_", OUTPUT_TAG)
table_stem   <- paste0("enrichment_table_", OUTPUT_TAG)
summary_stem <- paste0("summary_", OUTPUT_TAG)
objects_stem <- paste0("run_objects_", OUTPUT_TAG)

if (capabilities("cairo")) {
  grDevices::cairo_pdf(out_pdf(plot_stem), width = 12, height = 9)
  print(gg)
  grDevices::dev.off()
} else {
  grDevices::pdf(out_pdf(plot_stem), width = 12, height = 9, useDingbats = FALSE)
  print(gg)
  grDevices::dev.off()
}

readr::write_csv(DF, out_csv(table_stem))
saveRDS(DF, out_rds(table_stem))

saveRDS(
  list(
    toggles = list(
      CHAIN_TARGET = CHAIN_TARGET,
      FEATURE_MODE = FEATURE_MODE,
      KEEP_ALLELES = KEEP_ALLELES,
      DATA_SOURCE = DATA_SOURCE,
      DATA_PATH = DATA_PATH,
      DATA_OBJECT_NAME = DATA_OBJECT_NAME,
      FILTER_INKT_MAIT = FILTER_INKT_MAIT,
      DOWNSAMPLE_JUNCTIONS = DOWNSAMPLE_JUNCTIONS,
      SEED_DOWNSAMPLE = SEED_DOWNSAMPLE,
      RANDOMIZE_SAMPLES = RANDOMIZE_SAMPLES,
      SEED_RANDOM = SEED_RANDOM,
      MIN_GROUP_COUNT = MIN_GROUP_COUNT,
      REQUIRE_BOTH_GROUPS = REQUIRE_BOTH_GROUPS,
      PLOT_KEEP_NS_TOP = PLOT_KEEP_NS_TOP,
      OUTPUT_DIR = OUTPUT_DIR,
      OUTPUT_TAG = OUTPUT_TAG,
      RUN_TIMESTAMP = RUN_TIMESTAMP
    ),
    hu_chain = hu_chain,
    alpha_only_ids = alpha_only_ids,
    beta_only_ids = beta_only_ids
  ),
  out_rds(objects_stem)
)

sink(out_txt(summary_stem))
cat("OUTPUT_DIR:            ", OUTPUT_DIR, "\n")
cat("RUN_TIMESTAMP:         ", RUN_TIMESTAMP, "\n")
cat("OUTPUT_TAG:            ", OUTPUT_TAG, "\n\n")

cat("Chain target:          ", CHAIN_TARGET, "\n")
cat("Feature mode:          ", FEATURE_MODE, "\n")
cat("Keep alleles:          ", KEEP_ALLELES, "\n")
cat("Downsample junctions:  ", DOWNSAMPLE_JUNCTIONS, "\n")
cat("Randomize samples:     ", RANDOMIZE_SAMPLES, "\n")
cat("Min group count:       ", MIN_GROUP_COUNT, "\n")
cat("Require both groups:   ", REQUIRE_BOTH_GROUPS, "\n\n")

cat("Unique junctions per group (post-filter/post-downsample):\n")
print(hu_chain %>% distinct(group, junction) %>% count(group, name = "n_unique_junctions"))

cat("\n#Features tested:", nrow(DF), "\n")

cat("\nTop features (FDR < 0.05):\n")
print(
  DF %>%
    filter(fdr < 0.05) %>%
    select(Feature, OR, alpha_prop, beta_prop, fdr) %>%
    head(30),
  row.names = FALSE
)
sink()

cat("\nSaved:\n",
    "  Plot PDF:  ", out_pdf(plot_stem), "\n",
    "  Table:     ", out_csv(table_stem), "\n",
    "  Table RDS: ", out_rds(table_stem), "\n",
    "  Summary:   ", out_txt(summary_stem), "\n",
    "  Objects:   ", out_rds(objects_stem), "\n",
    sep = "")
