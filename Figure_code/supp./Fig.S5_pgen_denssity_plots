

## =========================
## VDJdb → metaclones → pgen density plot (3 rows x 2 cols) + vlines
## Rows: alphaCentric, betaCentric, non-centric
## Cols: TRA, TRB
## X axis: Generational probability, -log10(pgen)
## Palette: cbPalette = c('#66c2a5','#fc8d62','#7570b3')
## Vlines:
##  - dotted  = median non-centric TRA  (-log10 pgen)
##  - solid   = median non-centric TRB  (-log10 pgen)
## =========================

rm(list = ls())

suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(readr)
  library(ggplot2)
})

msg <- function(...) message(sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), paste0(...)))

## ---------- USER SETTINGS ----------
VDJDB_PATH <- "/Users/peterlinsley/Dropbox/Mac/Desktop/Shubhams_paper_code/data/vdjdb.txt"
PGEN_PATH  <- "/Users/peterlinsley/Dropbox/Mac/Desktop/Shubhams_paper_code/data/igor_output_combined.csv"
OUT_DIR    <- "/Users/peterlinsley/Dropbox/Mac/Desktop/Shubhams_paper_code/Figure_PDFs"

SEED <- 42
N_CELLS_PER_TYPE <- 300          # cells = unique libid per type
USE_CANONICAL_CDR3 <- TRUE
X_LIMS <- c(5, 25)               # x-axis limits for -log10(pgen)
## ----------------------------------

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

###############################################################################
# 1) LOAD + STANDARDIZE VDJdb
###############################################################################
stopifnot(file.exists(VDJDB_PATH))
db <- read.delim(VDJDB_PATH, stringsAsFactors = FALSE)

names(db) <- gsub("^complex\\.id$", "libid", names(db))
names(db) <- gsub("^v\\.segm$", "v_gene", names(db))
names(db) <- gsub("^j\\.segm$", "j_gene", names(db))
names(db) <- gsub("^antigen\\.species$", "study_group", names(db))

if (!"junction" %in% names(db)) {
  if ("cdr3" %in% names(db)) db$junction <- db$cdr3
  else if ("cdr3fix" %in% names(db)) db$junction <- db$cdr3fix
  else stop("Could not find junction column: expected 'junction' or 'cdr3' or 'cdr3fix'.")
}

need <- c("libid","v_gene","j_gene","junction","species")
miss <- setdiff(need, names(db))
if (length(miss)) stop("Missing required columns in VDJdb: ", paste(miss, collapse=", "))

hu <- db %>%
  dplyr::filter(species == "HomoSapiens") %>%
  dplyr::mutate(
    libid = as.character(libid),
    gene = dplyr::case_when(
      stringr::str_detect(v_gene, "TRA") ~ "TRA",
      stringr::str_detect(v_gene, "TRB") ~ "TRB",
      TRUE ~ NA_character_
    )
  ) %>%
  dplyr::filter(
    libid != "0",
    gene %in% c("TRA","TRB"),
    !is.na(junction), junction != ""
  )

if (USE_CANONICAL_CDR3) {
  hu <- hu %>% dplyr::filter(stringr::str_detect(junction, "^C.*[FW]$"))
}

msg("After base filters: rows=", nrow(hu), " libids=", dplyr::n_distinct(hu$libid))

###############################################################################
# 2) REQUIRE EXACTLY 2 ROWS PER LIBID, BOTH CHAINS, TWO DISTINCT JUNCTIONS
###############################################################################
counts <- hu %>% dplyr::count(libid, name = "n")
twoChains <- counts %>% dplyr::filter(n == 2) %>% dplyr::pull(libid)
hu2 <- hu %>% dplyr::filter(libid %in% twoChains)

ids_both <- hu2 %>%
  dplyr::group_by(libid) %>%
  dplyr::summarize(has_TRA = any(gene=="TRA"), has_TRB = any(gene=="TRB"), .groups="drop") %>%
  dplyr::filter(has_TRA & has_TRB) %>%
  dplyr::pull(libid)
hu2 <- hu2 %>% dplyr::filter(libid %in% ids_both)

twoJuncs <- hu2 %>%
  dplyr::group_by(libid) %>%
  dplyr::summarize(nj = dplyr::n_distinct(junction), .groups="drop") %>%
  dplyr::filter(nj == 2) %>%
  dplyr::pull(libid)
hu2 <- hu2 %>% dplyr::filter(libid %in% twoJuncs)

msg("After pairing constraints: rows=", nrow(hu2), " libids=", dplyr::n_distinct(hu2$libid))

###############################################################################
# 3) BUILD TRA/TRB PAIRS & CALL METACLONES
###############################################################################
tra <- hu2 %>% dplyr::filter(gene=="TRA") %>% dplyr::transmute(libid, junction_TRA = junction)
trb <- hu2 %>% dplyr::filter(gene=="TRB") %>% dplyr::transmute(libid, junction_TRB = junction)
join1 <- dplyr::inner_join(tra, trb, by="libid")

alpha_cent <- join1 %>%
  dplyr::distinct(junction_TRA, junction_TRB) %>%
  dplyr::count(junction_TRA, name="nTRB") %>%
  dplyr::filter(nTRB > 1)

beta_cent <- join1 %>%
  dplyr::distinct(junction_TRB, junction_TRA) %>%
  dplyr::count(junction_TRB, name="nTRA") %>%
  dplyr::filter(nTRA > 1)

alpha_ids <- join1 %>% dplyr::semi_join(alpha_cent, by="junction_TRA") %>% dplyr::pull(libid) %>% unique()
beta_ids  <- join1 %>% dplyr::semi_join(beta_cent,  by="junction_TRB") %>% dplyr::pull(libid) %>% unique()

msg("alphaCentric libids: ", length(alpha_ids))
msg("betaCentric  libids: ", length(beta_ids))

###############################################################################
# 4) LOAD PGEN + BUILD comb2 (your plotting frame)
###############################################################################
stopifnot(file.exists(PGEN_PATH))
values <- read.csv(PGEN_PATH, stringsAsFactors = FALSE)

if (!"TCR_name" %in% names(values)) stop("PGEN file missing 'TCR_name'.")
values$libid <- gsub("^TCR_", "", values$TCR_name)

need_pgen <- c("libid","TRA_pgen","TRB_pgen")
miss_pgen <- setdiff(need_pgen, names(values))
if (length(miss_pgen)) stop("PGEN file missing columns: ", paste(miss_pgen, collapse=", "))

# Add chain-specific pgen
tra_pgen <- hu2 %>% dplyr::filter(gene=="TRA") %>%
  dplyr::mutate(pgen = values$TRA_pgen[match(libid, values$libid)])
trb_pgen <- hu2 %>% dplyr::filter(gene=="TRB") %>%
  dplyr::mutate(pgen = values$TRB_pgen[match(libid, values$libid)])

comb <- dplyr::bind_rows(tra_pgen, trb_pgen) %>%
  dplyr::filter(!is.na(pgen)) %>%
  dplyr::mutate(
    type = dplyr::case_when(
      libid %in% alpha_ids ~ "alphaCentric",
      libid %in% beta_ids  ~ "betaCentric",
      TRUE ~ "non-centric"
    ),
    chainType = ifelse(gene=="TRA", "TRA chain", "TRB chain")
  )

###############################################################################
# 5) DOWNSAMPLE TO 300 CELLS PER TYPE (LIBID), KEEP BOTH CHAINS
###############################################################################
set.seed(SEED)

sample_ids <- function(ids, n) {
  ids <- unique(ids)
  if (!length(ids)) return(character(0))
  sample(ids, size = min(n, length(ids)), replace = FALSE)
}

ids_alpha <- sample_ids(comb %>% dplyr::filter(type=="alphaCentric") %>% dplyr::pull(libid) %>% unique(), N_CELLS_PER_TYPE)
ids_beta  <- sample_ids(comb %>% dplyr::filter(type=="betaCentric")  %>% dplyr::pull(libid) %>% unique(), N_CELLS_PER_TYPE)
ids_non   <- sample_ids(comb %>% dplyr::filter(type=="non-centric")  %>% dplyr::pull(libid) %>% unique(), N_CELLS_PER_TYPE)

keep_ids <- unique(c(ids_alpha, ids_beta, ids_non))

comb2 <- comb %>% dplyr::filter(libid %in% keep_ids)

# keep only libids with both chains after downsampling
both_ids <- comb2 %>%
  dplyr::group_by(libid) %>%
  dplyr::summarize(has_TRA = any(gene=="TRA"), has_TRB = any(gene=="TRB"), .groups="drop") %>%
  dplyr::filter(has_TRA & has_TRB) %>%
  dplyr::pull(libid)

comb2 <- comb2 %>% dplyr::filter(libid %in% both_ids)

msg("Downsampled libids per type:")
print(comb2 %>% dplyr::distinct(libid, type) %>% dplyr::count(type))

###############################################################################
# 6) VLINES EXACTLY LIKE YOUR EXAMPLE (computed from comb2)
###############################################################################
abSub1 <- dplyr::filter(comb2, type == "non-centric" & gene == "TRA")
abSub2 <- dplyr::filter(comb2, type == "non-centric" & gene == "TRB")

vline_TRA_non <- stats::median(-log10(abSub1$pgen), na.rm = TRUE)  # dotted
vline_TRB_non <- stats::median(-log10(abSub2$pgen), na.rm = TRUE)  # solid

###############################################################################
# 7) PLOT (3x2 FACETS) + STYLE
###############################################################################
comb2$type <- factor(comb2$type, levels = c("alphaCentric", "betaCentric", "non-centric"))
comb2$chainType <- factor(comb2$chainType, levels = c("TRA chain", "TRB chain"))

cbPalette <- c("#66c2a5", "#fc8d62", "#7570b3")
names(cbPalette) <- c("alphaCentric", "betaCentric", "non-centric")

p <- ggplot2::ggplot(comb2, ggplot2::aes(x = -log10(pgen))) +
  ggplot2::geom_density(ggplot2::aes(fill = type), alpha = 0.5) +
  ggplot2::labs(
    x = "\nGenerational probability,\n-log10(pgen)",
    y = "Density\n",
    fill = "type"
  ) +
  ggplot2::scale_fill_manual(values = cbPalette) +
  ggplot2::facet_grid(type ~ chainType) +
  ggplot2::geom_vline(xintercept = vline_TRA_non, linetype = "dotted", linewidth = 1.5) +
  ggplot2::geom_vline(xintercept = vline_TRB_non, linetype = "solid",  linewidth = 1.5) +
  ggplot2::theme_bw(base_size = 24) +
  ggplot2::theme(
    panel.grid.major = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),
    legend.key = ggplot2::element_blank()
  )

if (!is.null(X_LIMS)) {
  p <- p + ggplot2::scale_x_continuous(limits = X_LIMS)
}

stamp <- gsub("-", "", as.character(Sys.Date()))
out_pdf <- file.path(OUT_DIR, paste0("pgen_density_3x2_vlines_", stamp, ".pdf"))
ggplot2::ggsave(out_pdf, p, width = 10, height = 8, units = "in", dpi = 300)
msg("Saved: ", out_pdf)

# Save plot data
out_csv <- file.path(OUT_DIR, paste0("pgen_density_3x2_input_", stamp, ".csv"))
readr::write_csv(comb2 %>% dplyr::select(libid, gene, chainType, type, pgen), out_csv)
msg("Saved: ", out_csv)

msg("Done.")
