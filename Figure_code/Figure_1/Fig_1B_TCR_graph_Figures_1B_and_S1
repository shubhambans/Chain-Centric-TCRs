## =========================
## VDJdb → TCRgraph + clean CSV summaries
## (≥2 chains, must have α & β; global metaclone detection; per-species graphs)
## =========================

rm(list = ls())

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(stringr); library(readr)
  library(ggplot2); library(ggraph); library(igraph); library(viridis); library(tcrGraph)
})

`%||%` <- function(x, y) if (is.null(x) || (length(x)==1 && is.na(x))) y else x

## ---------- USER TOGGLES ----------
VDJDB_PATH   <- "/Users/peterlinsley/Dropbox/Mac/Desktop/Shubhams_paper_code/data/vdjdb.txt"
OUTPUT_DIR   <- "~/Desktop/Shubhams_paper_code/Figure_PDFs"
SPECIES_SET  <- "SARS-CoV-2"                  # or c("SARS-CoV-2","CMV", "InfluenzaA", "HomoSapiens", "EBV", "HCV")
DOWNSAMPLE_N <- 350                    # set NA to disable downsampling
SEED         <- 42
USE_CANONICAL_CDR3 <- TRUE            # TRUE enforces ^C.*[FW]$
PREFER_CDR3FIX     <- FALSE             # prefer cdr3fix over cdr3 when available
SAVE_INPUT_TABLE   <- TRUE             # also save the downsampled input used to draw graph
SAVE_SUMMARIES     <- TRUE             # write summary + top metaclone CSVs
TOP_N_METACLONES   <- 25               # top N α/β-centric junctions per species
## ----------------------------------

dir.create(path.expand(OUTPUT_DIR), showWarnings = FALSE, recursive = TRUE)

## ---------- Helpers ----------
read_vdjdb <- function(path){
  tryCatch(suppressWarnings(readr::read_tsv(path, show_col_types = FALSE)),
           error = function(e) read.delim(path, stringsAsFactors = FALSE, sep = "\t"))
}

std_cols <- function(df){
  nm <- names(df)
  need <- c("complex.id","v.segm","j.segm","species")
  miss <- setdiff(need, nm)
  if (length(miss)) stop(sprintf("Missing required columns: %s", paste(miss, collapse=", ")))

  junction_col <- if (PREFER_CDR3FIX && "cdr3fix" %in% nm) "cdr3fix"
  else if ("cdr3" %in% nm) "cdr3"
  else if ("cdr3fix" %in% nm) "cdr3fix" else NA_character_
  if (is.na(junction_col)) stop("Neither 'cdr3fix' nor 'cdr3' present for junctions.")

  df %>%
    dplyr::rename(
      libid           = `complex.id`,
      v_gene          = `v.segm`,
      j_gene          = `j.segm`,
      junction        = !!sym(junction_col),
      antigen_species = `antigen.species`,
      antigen_gene    = `antigen.gene`
    )
}

# Extract plain CDR3 strings if junction is JSON-like; otherwise return as-is
clean_cdr3 <- function(x) {
  x <- as.character(x)
  looks_json <- !is.na(x) & grepl("^\\s*\\{", x)
  out <- x
  if (any(looks_json)) {
    if (!requireNamespace("jsonlite", quietly = TRUE)) {
      stop("Please install 'jsonlite' (install.packages('jsonlite')) to clean JSON-like CDR3 fields.")
    }
    out[looks_json] <- vapply(x[looks_json], function(s) {
      val <- tryCatch(jsonlite::fromJSON(s), error = function(e) NULL)
      if (is.null(val)) s else as.character(val$cdr3 %||% s)
    }, character(1))
  }
  trimws(out)
}

clean_alleles <- function(x){
  gsub("//\\*\\d+$","", gsub("\\*\\d+$","", x))
}

write_csv_safe <- function(df, path){
  # always write, even if empty; ensures headers exist
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  readr::write_csv(df, path, na = "")
  message("   saved: ", path)
}

## ---------- Load & normalize ----------
db <- read_vdjdb(VDJDB_PATH) %>% as.data.frame(stringsAsFactors = FALSE) %>% std_cols()
message("Loaded cols: ", paste(names(db), collapse=", "))

# Keep human; normalize label
hu <- db %>% filter(species %in% c("HomoSapiens","Homo Sapiens","HSapiens"))
hu$species <- "HomoSapiens"

# Infer chain if 'gene' missing
if (!"gene" %in% names(hu)) {
  hu <- hu %>% mutate(gene = case_when(
    str_detect(v_gene, "TRA") ~ "TRA",
    str_detect(v_gene, "TRB") ~ "TRB",
    TRUE ~ NA_character_
  ))
}

hu <- hu %>%
  filter(!is.na(libid), libid != 0,
         !is.na(junction), junction != "")

if (USE_CANONICAL_CDR3) {
  hu <- hu %>% filter(str_detect(junction, "^C.*[FW]$"))
}

# Require BOTH chains per libid (≥1 TRA AND ≥1 TRB)
hu <- hu %>%
  mutate(chainType = gene %||% if_else(str_detect(v_gene,"TRA"),"TRA",
                                       if_else(str_detect(v_gene,"TRB"),"TRB",NA_character_))) %>%
  filter(chainType %in% c("TRA","TRB"))

ids_both <- hu %>%
  group_by(libid) %>%
  summarize(has_TRA = any(chainType=="TRA"),
            has_TRB = any(chainType=="TRB"),
            .groups = "drop") %>%
  filter(has_TRA & has_TRB) %>%
  pull(libid)

hu2 <- hu %>% filter(libid %in% ids_both)

message("After filters: rows=", nrow(hu2),
        ", libids=", n_distinct(hu2$libid),
        ", with both chains=", length(unique(ids_both)),
        ", TRA rows=", sum(hu2$chainType=="TRA"),
        ", TRB rows=", sum(hu2$chainType=="TRB"))

## ---------- Build global TRA–TRB pairs (for species/metaclone baseline) ----------
tra_min_global <- hu2 %>%
  filter(chainType == "TRA") %>%
  transmute(libid,
            TRA_junc    = clean_cdr3(junction),
            TRA_v       = v_gene,
            TRA_j       = j_gene,
            TRA_species = antigen_species,
            TRA_gene    = antigen_gene)

trb_min_global <- hu2 %>%
  filter(chainType == "TRB") %>%
  transmute(libid,
            TRB_junc    = clean_cdr3(junction),
            TRB_v       = v_gene,
            TRB_j       = j_gene,
            TRB_species = antigen_species,
            TRB_gene    = antigen_gene)

pairs_global <- inner_join(tra_min_global, trb_min_global, by = "libid")

# Metaclones (global, across libids)
alpha_cent_global <- pairs_global %>% distinct(TRA_junc, TRB_junc) %>%
  count(TRA_junc, name="nTRB") %>% filter(nTRB > 1)

beta_cent_global  <- pairs_global %>% distinct(TRB_junc, TRA_junc) %>%
  count(TRB_junc, name="nTRA") %>% filter(nTRA > 1)

to_plot_global <- bind_rows(
  pairs_global %>%
    semi_join(alpha_cent_global, by="TRA_junc") %>%
    transmute(libid, species2 = dplyr::coalesce(TRA_species, "other"), type = "1TRA-multiTRB"),
  pairs_global %>%
    semi_join(beta_cent_global, by="TRB_junc") %>%
    transmute(libid, species2 = dplyr::coalesce(TRB_species, "other"), type = "multiTRA-1TRB")
) %>% distinct()

bucket_keep <- c("CMV","EBV","InfluenzaA","SARS-CoV-2","HCV","Mtb","YFV","HomoSapiens")
to_plot_global <- to_plot_global %>%
  mutate(species2 = if_else(species2 %in% bucket_keep, species2, "other"))

message("Species buckets present (post-filter, global):")
print(to_plot_global %>% count(species2, sort = TRUE))

plot_species <- if (identical(SPECIES_SET, "ALL")) {
  sort(unique(to_plot_global$species2))
} else {
  xs <- intersect(SPECIES_SET, unique(to_plot_global$species2))
  if (!length(xs)) {
    warning(
      paste0(
        "Requested species not found after filtering. Falling back to available buckets: ",
        paste(sort(unique(to_plot_global$species2)), collapse = ", ")
      ),
      call. = FALSE
    )
    xs <- sort(unique(to_plot_global$species2))
  }
  xs
}
message("Species to plot: ", paste(plot_species, collapse = ", "))

## ---------- Plotting + per-species clean CSVs ----------
TRAcol <- viridis(1, begin = 0.12)
TRBcol <- viridis(1, begin = 0.88)
nodePalette <- c(TRA = TRAcol, TRB = TRBcol)

make_species_graph_and_csvs <- function(spec_label){
  message(sprintf("→ %s", spec_label))

  keep_ids <- to_plot_global %>% filter(species2 == spec_label) %>% pull(libid) %>% unique()
  if (!length(keep_ids)) { message("   (no libids)"); return(invisible(NULL)) }

  # Input df (downsampled libids if requested)
  df <- hu2 %>%
    filter(libid %in% keep_ids) %>%
    transmute(
      libid     = as.character(libid),
      chainType = chainType,
      v_gene    = clean_alleles(v_gene),
      j_gene    = clean_alleles(j_gene),
      junction  = clean_cdr3(junction)
    )

  if (!is.na(DOWNSAMPLE_N)) {
    set.seed(SEED)
    ids <- sample(unique(df$libid), size = min(DOWNSAMPLE_N, length(unique(df$libid))), replace = FALSE)
    df <- df %>% filter(libid %in% ids)
  } else {
    ids <- unique(df$libid)
  }

  # Within-sample pairs & top metaclones (partners = plain CDR3s)
  tra_s <- df %>% filter(chainType=="TRA") %>%
    transmute(libid, TRA_junc = junction, TRA_v = v_gene, TRA_j = j_gene)
  trb_s <- df %>% filter(chainType=="TRB") %>%
    transmute(libid, TRB_junc = junction, TRB_v = v_gene, TRB_j = j_gene)

  pairs_s <- inner_join(tra_s, trb_s, by="libid") %>%
    distinct(TRA_junc, TRB_junc, .keep_all = TRUE)

  alpha_counts <- pairs_s %>%
    count(TRA_junc, name="nTRB") %>%
    arrange(desc(nTRB), TRA_junc)

  alpha_partners <- pairs_s %>%
    group_by(TRA_junc) %>%
    summarize(TRB_partners = paste(sort(unique(TRB_junc)), collapse="; "), .groups="drop")

  alpha_top <- alpha_counts %>%
    left_join(alpha_partners, by="TRA_junc") %>%
    slice_head(n = TOP_N_METACLONES)

  beta_counts <- pairs_s %>%
    count(TRB_junc, name="nTRA") %>%
    arrange(desc(nTRA), TRB_junc)

  beta_partners <- pairs_s %>%
    group_by(TRB_junc) %>%
    summarize(TRA_partners = paste(sort(unique(TRA_junc)), collapse="; "), .groups="drop")

  beta_top <- beta_counts %>%
    left_join(beta_partners, by="TRB_junc") %>%
    slice_head(n = TOP_N_METACLONES)

  # Summary row
  summary_row <- tibble(
    species      = spec_label,
    date         = as.character(Sys.Date()),
    n_rows       = nrow(df),
    n_libids     = n_distinct(df$libid),
    n_TRA_rows   = sum(df$chainType=="TRA"),
    n_TRB_rows   = sum(df$chainType=="TRB"),
    n_TRA_junc   = n_distinct(df$junction[df$chainType=="TRA"]),
    n_TRB_junc   = n_distinct(df$junction[df$chainType=="TRB"])
  )

  # ----- GUARANTEED SAVE PATHS (per-species folder) -----
  stamp <- format(Sys.Date(), "%y%m%d")
  spec_dir <- file.path(path.expand(OUTPUT_DIR), spec_label)
  dir.create(spec_dir, showWarnings = FALSE, recursive = TRUE)

  out_summary <- file.path(spec_dir, sprintf("%s_downsample_summary_%s.csv", spec_label, stamp))
  out_alpha   <- file.path(spec_dir, sprintf("%s_top_alphaCentric_%s.csv",   spec_label, stamp))
  out_beta    <- file.path(spec_dir, sprintf("%s_top_betaCentric_%s.csv",    spec_label, stamp))
  out_input   <- file.path(spec_dir, sprintf("%s_TCRgraph_input_%s.csv",     spec_label, stamp))
  out_pdf     <- file.path(spec_dir, sprintf("%s_TCRgraph_%s.pdf",           spec_label, stamp))

  if (SAVE_SUMMARIES) {
    write_csv_safe(summary_row, out_summary)
    write_csv_safe(alpha_top,   out_alpha)   # TRA_junc, nTRB, TRB_partners
    write_csv_safe(beta_top,    out_beta)    # TRB_junc, nTRA, TRA_partners
    if (SAVE_INPUT_TABLE) write_csv_safe(df, out_input)
  }

  # Build & save TCRgraph
  tg <- tcrGraph::makeTcrGraph(df, link = "junction")
  gi <- tcrGraphToIgraph(tg)
  V(gi)$`Number of Cells` <- V(gi)$value
  V(gi)$Chain <- V(gi)$group

  p <- ggraph(gi, layout = "fr") +
    geom_edge_link(edge_width = 1, edge_colour = "grey70") +
    geom_node_point(aes(fill = Chain, size = `Number of Cells`), shape = 21, colour = "black") +
    scale_fill_manual(values = nodePalette) +
    ggtitle(sprintf("%s TCRs (meta-clone enriched)", spec_label)) +
    theme_void(base_size = 12) +
    theme(legend.position = "right", plot.title = element_text(face = "bold"),
	 panel.border = element_rect(colour = "black", fill = NA, linewidth = 1))

  ggsave(out_pdf, p, width = 10, height = 9, units = "in", dpi = 300)
  message("   saved: ", out_pdf)

  invisible(list(summary = summary_row, alpha = alpha_top, beta = beta_top))
}

invisible(lapply(plot_species, make_species_graph_and_csvs))
message("Done.")
