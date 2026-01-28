###############################################################################
# END-TO-END: Su_et_al_filtered_TCRs.csv
#
# Metaclone definition toggles:
#   DESIGNATION_UNIVERSE: "ALL" / "HEALTHY" / "COVID"   (TCRs used to DEFINE centricity)
#   METACLONE_SCOPE     : "GLOBAL" / "LOCAL"            (define centricity pooled vs within CD4/CD8)
#   MIN_CELLS_FOR_CENTRIC: require centric junction seen in >= this many cells (barcode/libid)
#   CELL_ID_COL         : "barcode" or "libid" (must exist)
#
# Downstream analysis:
#   - Always evaluates fractions in HEALTHY donors (as in your original script)
#   - Computes per-donor fractions (alphaCentric, betaCentric) in CD4 and CD8
#   - Computes:
#       (A) CD4 vs CD8 p-values (paired by donor) within each type + effect sizes
#       (B) alpha vs beta p-values (paired by donor) within each cellType
#   - Merges (A) into a single summary table
#   - Adds asterisks on plot for (B): alpha vs beta within CD4 facet, and within CD8 facet
#   - Saves plot PDF + PNG
###############################################################################

rm(list = ls())

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(stringr)
})

###############################################################################
# 0) USER TOGGLES
###############################################################################
VISIT_TOGGLE <- "ALL"     # "ALL" or e.g. "BL", "single_tp"

# What set of TCRs to use to DEFINE metaclones:
# "ALL" / "HEALTHY" / "COVID"
DESIGNATION_UNIVERSE <- "ALL"

# "GLOBAL" = define centricity pooled across CD4+CD8
# "LOCAL"  = define centricity separately within CD4 vs CD8
METACLONE_SCOPE <- "LOCAL"

# Require centric junction to be observed in >1 cell
MIN_CELLS_FOR_CENTRIC <- 2     # ">1 cell" => 2
CELL_ID_COL <- "barcode"       # "barcode" (default) or "libid" if present

# P-values: paired tests by donor
PVAL_TEST_BETWEEN_CELLTYPES <- "wilcox"  # "wilcox" or "t"
PVAL_TEST_WITHIN_CELLTYPE   <- "wilcox"  # "wilcox" or "t"

DATA_DIR <- "/Users/peterlinsley/Desktop/Shubhams_paper_code/data"
CSV_FILE <- "Su_et_al_filtered_TCRs.csv"

OUT_DIR <- "~/Desktop/Shubhams_paper_code/Figure_PDFs"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

OUT_STEM <- paste0(
  "Su_defineMeta_", DESIGNATION_UNIVERSE,
  "_scope_", METACLONE_SCOPE,
  "_minCells", MIN_CELLS_FOR_CENTRIC, "_cellID_", CELL_ID_COL,
  "_thenHealthy_CD4_CD8",
  ifelse(VISIT_TOGGLE == "ALL", "_allVisits", paste0("_", VISIT_TOGGLE))
)

OUT_PDF <- file.path(OUT_DIR, paste0(OUT_STEM, ".pdf"))
OUT_PNG <- file.path(OUT_DIR, paste0(OUT_STEM, ".png"))

BASE_SIZE <- 22
PLOT_W    <- 8
PLOT_H    <- 5

###############################################################################
# 1) LOAD DATA
###############################################################################
setwd(DATA_DIR)
su.full <- read.csv(CSV_FILE, stringsAsFactors = FALSE)

req <- c("barcode","junction","donor","cellType","visit","chainType")
if (CELL_ID_COL != "barcode") req <- unique(c(req, CELL_ID_COL))
miss <- setdiff(req, names(su.full))
if (length(miss)) stop("Missing required columns: ", paste(miss, collapse = ", "))

if (!DESIGNATION_UNIVERSE %in% c("ALL","HEALTHY","COVID")) {
  stop("DESIGNATION_UNIVERSE must be one of: ALL, HEALTHY, COVID")
}
if (!METACLONE_SCOPE %in% c("GLOBAL","LOCAL")) {
  stop("METACLONE_SCOPE must be 'GLOBAL' or 'LOCAL'")
}
if (!is.numeric(MIN_CELLS_FOR_CENTRIC) || MIN_CELLS_FOR_CENTRIC < 1) {
  stop("MIN_CELLS_FOR_CENTRIC must be a positive integer (>=1).")
}
if (!PVAL_TEST_BETWEEN_CELLTYPES %in% c("wilcox","t")) stop("Bad PVAL_TEST_BETWEEN_CELLTYPES")
if (!PVAL_TEST_WITHIN_CELLTYPE %in% c("wilcox","t")) stop("Bad PVAL_TEST_WITHIN_CELLTYPE")

su.full <- su.full %>%
  mutate(
    donor     = as.character(donor),
    cellType  = tolower(as.character(cellType)),
    visit     = as.character(visit),
    chainType = toupper(as.character(chainType)),
    barcode   = as.character(barcode),
    junction  = as.character(junction),
    cell_id   = as.character(.data[[CELL_ID_COL]])
  ) %>%
  filter(chainType %in% c("TRA","TRB")) %>%
  filter(cellType %in% c("cd4","cd8"))

###############################################################################
# 2) PAIRING (ALL CELLS FIRST; then we subset paired table for designation universe)
###############################################################################
paired_keys <- su.full %>%
  distinct(barcode, donor, visit, cellType, cell_id, chainType) %>%
  count(barcode, donor, visit, cellType, cell_id, name = "n_chain") %>%
  filter(n_chain >= 2) %>%
  select(barcode, donor, visit, cellType, cell_id)

su.paired <- su.full %>%
  semi_join(paired_keys, by = c("barcode","donor","visit","cellType","cell_id"))

tra_all <- su.paired %>%
  filter(chainType == "TRA") %>%
  transmute(barcode, donor, visit, cellType, cell_id, junction.tra = junction)

trb_all <- su.paired %>%
  filter(chainType == "TRB") %>%
  transmute(barcode, donor, visit, cellType, cell_id, junction.trb = junction)

paired_all <- inner_join(
  tra_all, trb_all,
  by = c("barcode","donor","visit","cellType","cell_id")
)

if (nrow(paired_all) == 0) stop("No paired TRA/TRB rows found after pairing step.")

###############################################################################
# 2.5) SUBSET PAIRED TABLE FOR METACLONE DESIGNATION UNIVERSE
###############################################################################
paired_for_designation <- paired_all
if (DESIGNATION_UNIVERSE == "HEALTHY") {
  paired_for_designation <- paired_for_designation %>%
    filter(grepl("Healthy", donor, ignore.case = TRUE))
} else if (DESIGNATION_UNIVERSE == "COVID") {
  paired_for_designation <- paired_for_designation %>%
    filter(grepl("COVID", donor, ignore.case = TRUE))
}
if (nrow(paired_for_designation) == 0) {
  stop("No paired rows available for DESIGNATION_UNIVERSE='", DESIGNATION_UNIVERSE,
       "'. Check donor labels (grep patterns).")
}

###############################################################################
# 3) METACLONE DESIGNATION (GLOBAL or LOCAL) USING paired_for_designation
###############################################################################
if (METACLONE_SCOPE == "GLOBAL") {

  alpha_ids <- paired_for_designation %>%
    group_by(junction.tra) %>%
    filter(
      n_distinct(junction.trb) > 1,
      n_distinct(cell_id) >= MIN_CELLS_FOR_CENTRIC
    ) %>%
    ungroup() %>%
    distinct(barcode) %>%
    pull(barcode)

  beta_ids <- paired_for_designation %>%
    group_by(junction.trb) %>%
    filter(
      n_distinct(junction.tra) > 1,
      n_distinct(cell_id) >= MIN_CELLS_FOR_CENTRIC
    ) %>%
    ungroup() %>%
    distinct(barcode) %>%
    pull(barcode)

  metaclones_defined <- bind_rows(
    tibble(barcode = alpha_ids, type = "alphaCentric"),
    tibble(barcode = beta_ids,  type = "betaCentric")
  ) %>%
    distinct(barcode, type) %>%
    left_join(paired_all %>% distinct(barcode, donor, visit, cellType),
              by = "barcode")

} else { # LOCAL within cellType

  alpha_local <- paired_for_designation %>%
    group_by(cellType, junction.tra) %>%
    filter(
      n_distinct(junction.trb) > 1,
      n_distinct(cell_id) >= MIN_CELLS_FOR_CENTRIC
    ) %>%
    ungroup() %>%
    distinct(barcode, cellType) %>%
    mutate(type = "alphaCentric")

  beta_local <- paired_for_designation %>%
    group_by(cellType, junction.trb) %>%
    filter(
      n_distinct(junction.tra) > 1,
      n_distinct(cell_id) >= MIN_CELLS_FOR_CENTRIC
    ) %>%
    ungroup() %>%
    distinct(barcode, cellType) %>%
    mutate(type = "betaCentric")

  metaclones_defined <- bind_rows(alpha_local, beta_local) %>%
    distinct(barcode, cellType, type) %>%
    left_join(paired_all %>% distinct(barcode, donor, visit, cellType),
              by = c("barcode","cellType"))
}

metaclones_allcells <- metaclones_defined

###############################################################################
# 4) SUBSET TO HEALTHY (and optional VISIT_TOGGLE) FOR DOWNSTREAM FRACTIONS
###############################################################################
paired_hc <- paired_all %>%
  filter(grepl("Healthy", donor, ignore.case = TRUE))

metaclones_hc <- metaclones_allcells %>%
  filter(grepl("Healthy", donor, ignore.case = TRUE))

if (VISIT_TOGGLE != "ALL") {
  paired_hc     <- paired_hc %>% filter(visit == VISIT_TOGGLE)
  metaclones_hc <- metaclones_hc %>% filter(visit == VISIT_TOGGLE)
}
if (nrow(paired_hc) == 0) stop("No paired rows left after Healthy/visit filtering.")

###############################################################################
# 5) FRACTIONS PER DONOR × cellType × type
###############################################################################
total_cells <- paired_hc %>%
  distinct(barcode, donor, cellType) %>%
  count(donor, cellType, name = "n_total_paired")

plot_df <- metaclones_hc %>%
  count(donor, cellType, type, name = "n_metaclone") %>%
  complete(donor, cellType, type = c("alphaCentric","betaCentric"),
           fill = list(n_metaclone = 0)) %>%
  left_join(total_cells, by = c("donor","cellType")) %>%
  mutate(fraction = n_metaclone / n_total_paired) %>%
  filter(!is.na(fraction))

plot_df$type <- factor(plot_df$type, levels = c("alphaCentric","betaCentric"))
plot_df$cellType <- factor(plot_df$cellType, levels = c("cd4","cd8"),
                           labels = c("CD4+","CD8+"))

###############################################################################
# 6) STATISTICS
# 6A) Between cell types: CD4 vs CD8 within each type (paired by donor)
###############################################################################
wide_cd4_cd8 <- plot_df %>%
  mutate(donor = as.character(donor),
         type = as.character(type),
         cellType = as.character(cellType)) %>%
  select(donor, type, cellType, fraction) %>%
  tidyr::pivot_wider(names_from = cellType, values_from = fraction) %>%
  filter(is.finite(`CD4+`) & is.finite(`CD8+`)) %>%
  mutate(delta_cd4_minus_cd8 = `CD4+` - `CD8+`)

paired_test_cd4_cd8 <- function(df_one_type) {
  x <- df_one_type$`CD4+`
  y <- df_one_type$`CD8+`
  if (length(x) < 3) {
    return(tibble(n_paired_donors = length(x), method = "Insufficient paired donors",
                  statistic = NA_real_, p_value = NA_real_))
  }
  if (PVAL_TEST_BETWEEN_CELLTYPES == "wilcox") {
    tt <- suppressWarnings(stats::wilcox.test(x, y, paired = TRUE, exact = FALSE))
    tibble(n_paired_donors = length(x), method = tt$method,
           statistic = unname(tt$statistic), p_value = tt$p.value)
  } else {
    tt <- stats::t.test(x, y, paired = TRUE)
    tibble(n_paired_donors = length(x), method = tt$method,
           statistic = unname(tt$statistic), p_value = tt$p.value)
  }
}

# Effect sizes for CD4 vs CD8 within each type
effect_sizes_cd4_cd8 <- wide_cd4_cd8 %>%
  group_by(type) %>%
  summarize(
    n_paired_donors = n(),
    median_CD4_minus_CD8 = median(delta_cd4_minus_cd8),
    mean_CD4_minus_CD8   = mean(delta_cd4_minus_cd8),
    cohens_dz = mean(delta_cd4_minus_cd8) / sd(delta_cd4_minus_cd8),
    cliffs_delta = {
      s <- sign(delta_cd4_minus_cd8)
      (sum(s > 0) - sum(s < 0)) / length(s)
    },
    .groups = "drop"
  )

pvals_cd4_cd8 <- wide_cd4_cd8 %>%
  group_by(type) %>%
  group_modify(~ paired_test_cd4_cd8(.x)) %>%
  ungroup()

# MERGED SINGLE TABLE (requested): CD4 vs CD8 p-values + effect sizes
between_celltype_table <- pvals_cd4_cd8 %>%
  left_join(effect_sizes_cd4_cd8, by = "type") %>%
  mutate(
    p_adj_BH = p.adjust(p_value, method = "BH")
  ) %>%
  arrange(factor(type, levels = c("alphaCentric","betaCentric")))

cat("\n=== CD4 vs CD8 (paired by donor) within each metaclone type ===\n")
print(between_celltype_table)

###############################################################################
# 6B) Within each cell type: alphaCentric vs betaCentric (paired by donor)
#      These are the comparisons you want to annotate on the plot facets.
###############################################################################
wide_alpha_beta <- plot_df %>%
  mutate(donor = as.character(donor),
         type = as.character(type),
         cellType = as.character(cellType)) %>%
  select(donor, cellType, type, fraction) %>%
  tidyr::pivot_wider(names_from = type, values_from = fraction) %>%
  filter(is.finite(alphaCentric) & is.finite(betaCentric)) %>%
  mutate(delta_alpha_minus_beta = alphaCentric - betaCentric)

paired_test_alpha_beta <- function(df_one_celltype) {
  x <- df_one_celltype$alphaCentric
  y <- df_one_celltype$betaCentric
  if (length(x) < 3) {
    return(tibble(n_paired_donors = length(x), method = "Insufficient paired donors",
                  statistic = NA_real_, p_value = NA_real_))
  }
  if (PVAL_TEST_WITHIN_CELLTYPE == "wilcox") {
    tt <- suppressWarnings(stats::wilcox.test(x, y, paired = TRUE, exact = FALSE))
    tibble(n_paired_donors = length(x), method = tt$method,
           statistic = unname(tt$statistic), p_value = tt$p.value)
  } else {
    tt <- stats::t.test(x, y, paired = TRUE)
    tibble(n_paired_donors = length(x), method = tt$method,
           statistic = unname(tt$statistic), p_value = tt$p.value)
  }
}

within_celltype_table <- wide_alpha_beta %>%
  group_by(cellType) %>%
  group_modify(~ paired_test_alpha_beta(.x)) %>%
  ungroup() %>%
  mutate(p_adj_BH = p.adjust(p_value, method = "BH"))

cat("\n=== Within each cell type: alphaCentric vs betaCentric (paired by donor) ===\n")
print(within_celltype_table)

# p-value -> asterisks
p_to_stars <- function(p) {
  if (is.na(p)) return("n/a")
  if (p < 0.0001) return("****")
  if (p < 0.001)  return("***")
  if (p < 0.01)   return("**")
  if (p < 0.05)   return("*")
  "ns"
}

###############################################################################
# 7) PLOT + annotate alpha vs beta within each facet using asterisks
###############################################################################
# Base plot
p <- ggplot(plot_df, aes(x = type, y = fraction)) +
  geom_boxplot(outlier.shape = NA, fill = "grey85", linewidth = 0.6) +
  geom_jitter(width = 0.12, size = 2, alpha = 0.6) +
  facet_wrap(~ cellType, nrow = 1) +
  labs(
    x = "",
    y = "Fraction TCRs in metaclones\n",
    title = ""
  ) +
  theme_bw(base_size = BASE_SIZE) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5)
  )

# Build annotation dataframe for within-cell-type alpha vs beta
# We'll put a bracket spanning x=1..2 in each facet at y slightly above the max
ymax_by_facet <- plot_df %>%
  group_by(cellType) %>%
  summarize(ymax = max(fraction, na.rm = TRUE), .groups = "drop")

ann_within <- within_celltype_table %>%
  left_join(ymax_by_facet, by = "cellType") %>%
  mutate(
    stars = vapply(p_value, p_to_stars, character(1)),
    x1 = 1, x2 = 2,
    y  = ymax * 1.08,
    y_text = ymax * 1.10
  )

# Add bracket + stars (per facet)
p2 <- p +
  geom_segment(
    data = ann_within,
    aes(x = x1, xend = x2, y = y, yend = y),
    inherit.aes = FALSE,
    linewidth = 0.7
  ) +
  geom_segment(
    data = ann_within,
    aes(x = x1, xend = x1, y = y, yend = y * 0.985),
    inherit.aes = FALSE,
    linewidth = 0.7
  ) +
  geom_segment(
    data = ann_within,
    aes(x = x2, xend = x2, y = y, yend = y * 0.985),
    inherit.aes = FALSE,
    linewidth = 0.7
  ) +
  geom_text(
    data = ann_within,
    aes(x = 1.5, y = y_text, label = stars),
    inherit.aes = FALSE,
    size = 7
  ) +
  coord_cartesian(clip = "off") +
  theme(plot.margin = margin(10, 10, 20, 10))

print(p2)

###############################################################################
# 8) SAVE
###############################################################################
pdf_device <- if (capabilities("cairo")) cairo_pdf else pdf

ggsave(
  filename = OUT_PDF,
  plot     = p2,
  width    = PLOT_W,
  height   = PLOT_H,
  units    = "in",
  device   = pdf_device,
  bg       = "white"
)

ggsave(
  filename = OUT_PNG,
  plot     = p2,
  width    = PLOT_W,
  height   = PLOT_H,
  units    = "in",
  dpi      = 300,
  bg       = "white"
)

cat("\nSaved:\n", OUT_PDF, "\n", OUT_PNG, "\n", sep = "")

###############################################################################
# 9) SANITY CHECKS
###############################################################################
cat("\nSanity:\n")
cat("VISIT_TOGGLE:", VISIT_TOGGLE, "\n")
cat("DESIGNATION_UNIVERSE:", DESIGNATION_UNIVERSE, "\n")
cat("METACLONE_SCOPE:", METACLONE_SCOPE, "\n")
cat("CELL_ID_COL:", CELL_ID_COL, "\n")
cat("MIN_CELLS_FOR_CENTRIC:", MIN_CELLS_FOR_CENTRIC, "\n")
cat("All paired rows:", nrow(paired_all), "\n")
cat("Paired rows used for designation:", nrow(paired_for_designation), "\n")
cat("Healthy paired barcodes (after visit filter):", n_distinct(paired_hc$barcode), "\n")
cat("Healthy donors (after visit filter):", n_distinct(paired_hc$donor), "\n")

# donor drop counts for CD4 vs CD8 (within each type)
drop_diag <- wide_cd4_cd8 %>%
  group_by(type) %>%
  summarize(
    donors_used_paired = n(),
    .groups = "drop"
  )
cat("\nPaired donors used for CD4 vs CD8 tests (per type):\n")
print(drop_diag)
