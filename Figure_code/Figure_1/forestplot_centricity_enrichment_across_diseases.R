## =========================
## VDJdb → CSV summaries + Forest Plot (Odds Ratios)
## (≥2 chains, must have α & β; per-species CSVs; global Forest plot)
## =========================

rm(list = ls())

## If plyr is loaded, it can mask dplyr verbs. Detach defensively.
if ("package:plyr" %in% search()) {
  try(detach("package:plyr", unload = TRUE), silent = TRUE)
}

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(readr)
  library(ggplot2)
  library(forcats)
})

msg <- function(...) message(sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), paste0(...)))
`%||%` <- function(x, y) if (is.null(x) || (length(x)==1 && is.na(x))) y else x

## ---------- USER TOGGLES ----------
VDJDB_PATH         <- "~/Desktop/Shubhams_paper_code/data/vdjdb.txt"
OUTPUT_DIR         <- "~/Desktop/Shubhams_paper_code/Figure_PDFs"
SPECIES_SET        <- c("CMV","EBV","HCV","HomoSapiens", "InfluenzaA","SARS-CoV-2","other")    # or  c("CMV","EBV","InfluenzaA","SARS-CoV-2","HomoSapiens","HCV","Mtb","YFV") or "ALL"
DOWNSAMPLE_N       <- 350              # set NA to disable downsampling (for per-species input CSVs)
SEED               <- 42
USE_CANONICAL_CDR3 <- TRUE             # TRUE enforces ^C.*[FW]$
PREFER_CDR3FIX     <- FALSE            # prefer cdr3fix over cdr3 when available
SAVE_INPUT_TABLE   <- TRUE             # also save the downsampled input used for summaries
SAVE_SUMMARIES     <- TRUE             # write summary + top metaclone CSVs
TOP_N_METACLONES   <- 25               # top N α/β-centric junctions per species
## ----------------------------------

dir.create(path.expand(OUTPUT_DIR), showWarnings = FALSE, recursive = TRUE)

## ---------- IO helpers ----------
read_vdjdb <- function(path){
  tryCatch(suppressWarnings(readr::read_tsv(path, show_col_types = FALSE)),
           error = function(e) read.delim(path, stringsAsFactors = FALSE, sep = "\t"))
}
write_csv_safe <- function(df, path){
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  readr::write_csv(df, path, na = "")
  msg(paste("saved:", path))
}

## ---------- Ultra-robust column mapper ----------
std_cols <- function(df){
  df <- tibble::as_tibble(df)
  raw_names <- names(df)
  norm <- function(x) gsub("[^a-z0-9]", "", tolower(x))
  norm_names <- norm(raw_names)
  norm2orig <- setNames(raw_names, norm_names)

  pick <- function(candidates_norm, optional = FALSE, label = NULL){
    for (cand in candidates_norm) {
      hit <- norm2orig[[cand]]
      if (!is.null(hit)) return(hit)
    }
    if (optional) return(NULL)
    stop(
      sprintf("Missing required column %s. Available: %s",
              label %||% paste(candidates_norm, collapse = "/"),
              paste(raw_names, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  c_libid   <- pick(c("complexid","libid","pairid"),                                label = "libid/complex.id/pair.id")
  c_v       <- pick(c("vsegm","vgene","vseg","v","vsegment"),                        label = "V gene")
  c_j       <- pick(c("jsegm","jgene","jseg","j","jsegment"),                        label = "J gene")
  c_cdr3    <- if (PREFER_CDR3FIX)
                 pick(c("cdr3fix","cdr3","aacdr3","cdr3aa","cd3"),                   label = "CDR3")
               else
                 pick(c("cdr3","cdr3fix","aacdr3","cdr3aa","cd3"),                   label = "CDR3")
  c_species <- pick(c("species","speciestaxon","organism"), optional = TRUE)
  # prefer antigen.species; fall back to antigen.epitope if needed
  c_antsp   <- pick(c("antigenspecies","antigenepitope","studygroup","antigenepitopename"), optional = TRUE)
  c_antgene <- pick(c("antigengene"), optional = TRUE)
  c_gene    <- pick(c("gene","tcrchain","chain","genename"), optional = TRUE)

  tibble::tibble(
    libid           = df[[c_libid]],
    v_gene          = df[[c_v]],
    j_gene          = df[[c_j]],
    junction        = df[[c_cdr3]],
    antigen_species = if (!is.null(c_antsp))   df[[c_antsp]]   else NA_character_,
    antigen_gene    = if (!is.null(c_antgene)) df[[c_antgene]] else NA_character_,
    species         = if (!is.null(c_species)) df[[c_species]] else NA_character_,
    gene            = if (!is.null(c_gene))    df[[c_gene]]    else NA_character_
  )
}

## ---------- Utilities ----------
clean_cdr3 <- function(x) {
  x <- as.character(x)
  looks_json <- !is.na(x) & grepl("^\\s*\\{", x)
  out <- x
  if (any(looks_json)) {
    if (!requireNamespace("jsonlite", quietly = TRUE)) {
      stop("Please install 'jsonlite' (install.packages(\"jsonlite\")) to clean JSON-like CDR3 fields.")
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
safe_log2 <- function(x){
  x <- as.numeric(x); x[is.nan(x) | is.infinite(x)] <- NA_real_; log2(x)
}

## ---------- Load & standardize ----------
db_raw <- read_vdjdb(VDJDB_PATH)
msg("Raw columns: ", paste(names(db_raw), collapse=", "))

hu <- std_cols(db_raw)
stopifnot("libid" %in% names(hu), "junction" %in% names(hu))
hu <- dplyr::mutate(hu, libid = as.character(libid))
msg("Mapped columns: ", paste(names(hu), collapse=", "))

# Keep human; normalize label
hu <- dplyr::filter(hu, species %in% c("HomoSapiens","Homo Sapiens","HSapiens","homo sapiens","homosapiens"))
hu$species <- "HomoSapiens"

# Basic filters + canonical junctions (optional)
hu <- dplyr::filter(hu, !is.na(libid), libid != "0", !is.na(junction), junction != "")
if (USE_CANONICAL_CDR3) hu <- dplyr::filter(hu, stringr::str_detect(junction, "^C.*[FW]$"))

# Infer chain if 'gene' missing; keep TRA/TRB
if (!"gene" %in% names(hu) || all(is.na(hu$gene))) {
  hu <- dplyr::mutate(hu, gene = dplyr::case_when(
    stringr::str_detect(v_gene, "TRA") ~ "TRA",
    stringr::str_detect(v_gene, "TRB") ~ "TRB",
    TRUE ~ NA_character_
  ))
}
hu <- dplyr::mutate(hu, chainType = gene)
hu <- dplyr::filter(hu, chainType %in% c("TRA","TRB"))

## ---------- Require BOTH chains per libid (base R; no summarise/pull) ----------
by_lib <- split(hu$chainType, hu$libid)
ids_both <- names(Filter(function(v) any(v == "TRA") && any(v == "TRB"), by_lib))
ids_both <- as.character(ids_both)

hu2 <- dplyr::filter(hu, libid %in% ids_both)

msg(sprintf("After filters: rows=%d, libids=%d, with both chains=%d, TRA rows=%d, TRB rows=%d",
            nrow(hu2), dplyr::n_distinct(hu2$libid), length(unique(ids_both)),
            sum(hu2$chainType=="TRA"), sum(hu2$chainType=="TRB")))

## ---------- Build global TRA–TRB pairs ----------
tra_min_global <- hu2 %>%
  dplyr::filter(chainType == "TRA") %>%
  dplyr::transmute(libid,
                   TRA_junc    = clean_cdr3(junction),
                   TRA_v       = v_gene,
                   TRA_j       = j_gene,
                   TRA_species = antigen_species,
                   TRA_gene    = antigen_gene)

trb_min_global <- hu2 %>%
  dplyr::filter(chainType == "TRB") %>%
  dplyr::transmute(libid,
                   TRB_junc    = clean_cdr3(junction),
                   TRB_v       = v_gene,
                   TRB_j       = j_gene,
                   TRB_species = antigen_species,
                   TRB_gene    = antigen_gene)

pairs_global <- dplyr::inner_join(tra_min_global, trb_min_global, by = "libid")

# Metaclones (global, across libids)
alpha_cent_global <- pairs_global %>%
  dplyr::distinct(TRA_junc, TRB_junc) %>%
  dplyr::group_by(TRA_junc) %>%
  dplyr::summarise(nTRB = dplyr::n(), .groups = "drop") %>%
  dplyr::filter(nTRB > 1)

beta_cent_global  <- pairs_global %>%
  dplyr::distinct(TRB_junc, TRA_junc) %>%
  dplyr::group_by(TRB_junc) %>%
  dplyr::summarise(nTRA = dplyr::n(), .groups = "drop") %>%
  dplyr::filter(nTRA > 1)

to_plot_global <- dplyr::bind_rows(
  pairs_global %>%
    dplyr::semi_join(alpha_cent_global, by="TRA_junc") %>%
    dplyr::transmute(libid, species2 = dplyr::coalesce(TRA_species, "other"), type = "1TRA-multiTRB"),
  pairs_global %>%
    dplyr::semi_join(beta_cent_global, by="TRB_junc") %>%
    dplyr::transmute(libid, species2 = dplyr::coalesce(TRB_species, "other"), type = "multiTRA-1TRB")
) %>% dplyr::distinct()

# Keep named buckets, map others to "other"
bucket_keep <- SPECIES_SET
to_plot_global <- dplyr::mutate(to_plot_global,
                                species2 = ifelse(species2 %in% bucket_keep, species2, "other"),
                                libid    = as.character(libid))

# --- sanity guards ---
stopifnot(is.data.frame(to_plot_global))
stopifnot("libid" %in% names(to_plot_global))
if (!nrow(to_plot_global)) stop("No metaclone libids after filtering. Relax filters or check inputs.")

msg("Species buckets present (post-filter, global):")
species_counts <- to_plot_global %>%
  dplyr::group_by(species2) %>% dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
  dplyr::arrange(dplyr::desc(n))
print(species_counts)

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
msg("Species to use for CSVs: ", paste(plot_species, collapse = ", "))

## ---------- Per-species clean CSVs (no TCRgraph) ----------
make_species_csvs <- function(spec_label){
  msg(sprintf("→ %s", spec_label))

  keep_ids <- unique(to_plot_global$libid[to_plot_global$species2 == spec_label])
  if (!length(keep_ids)) { msg("   (no libids)"); return(invisible(NULL)) }

  df <- hu2 %>%
    dplyr::filter(libid %in% keep_ids) %>%
    dplyr::transmute(
      libid     = as.character(libid),
      chainType = chainType,
      v_gene    = clean_alleles(v_gene),
      j_gene    = clean_alleles(j_gene),
      junction  = clean_cdr3(junction)
    )

  if (!is.na(DOWNSAMPLE_N)) {
    set.seed(SEED)
    ids <- sample(unique(df$libid), size = min(DOWNSAMPLE_N, length(unique(df$libid))), replace = FALSE)
    df <- dplyr::filter(df, libid %in% ids)
  }

  tra_s <- dplyr::filter(df, chainType=="TRA") %>%
    dplyr::transmute(libid, TRA_junc = junction, TRA_v = v_gene, TRA_j = j_gene)
  trb_s <- dplyr::filter(df, chainType=="TRB") %>%
    dplyr::transmute(libid, TRB_junc = junction, TRB_v = v_gene, TRB_j = j_gene)

  pairs_s <- dplyr::inner_join(tra_s, trb_s, by="libid") %>%
    dplyr::distinct(TRA_junc, TRB_junc, .keep_all = TRUE)

  alpha_counts <- pairs_s %>%
    dplyr::group_by(TRA_junc) %>% dplyr::summarise(nTRB = dplyr::n(), .groups = "drop") %>%
    dplyr::arrange(dplyr::desc(nTRB), TRA_junc)
  alpha_partners <- pairs_s %>%
    dplyr::group_by(TRA_junc) %>%
    dplyr::summarise(TRB_partners = paste(sort(unique(TRB_junc)), collapse="; "), .groups = "drop")
  alpha_top <- dplyr::left_join(alpha_counts, alpha_partners, by="TRA_junc")
  if (nrow(alpha_top) > TOP_N_METACLONES) alpha_top <- utils::head(alpha_top, TOP_N_METACLONES)

  beta_counts <- pairs_s %>%
    dplyr::group_by(TRB_junc) %>% dplyr::summarise(nTRA = dplyr::n(), .groups = "drop") %>%
    dplyr::arrange(dplyr::desc(nTRA), TRB_junc)
  beta_partners <- pairs_s %>%
    dplyr::group_by(TRB_junc) %>%
    dplyr::summarise(TRA_partners = paste(sort(unique(TRA_junc)), collapse="; "), .groups = "drop")
  beta_top <- dplyr::left_join(beta_counts, beta_partners, by="TRB_junc")
  if (nrow(beta_top) > TOP_N_METACLONES) beta_top <- utils::head(beta_top, TOP_N_METACLONES)

  summary_row <- tibble::tibble(
    species      = spec_label,
    date         = as.character(Sys.Date()),
    n_rows       = nrow(df),
    n_libids     = dplyr::n_distinct(df$libid),
    n_TRA_rows   = sum(df$chainType=="TRA"),
    n_TRB_rows   = sum(df$chainType=="TRB"),
    n_TRA_junc   = dplyr::n_distinct(df$junction[df$chainType=="TRA"]),
    n_TRB_junc   = dplyr::n_distinct(df$junction[df$chainType=="TRB"])
  )

  ## copy-safe date stamp (YYYYMMDD; no % tokens)
  stamp <- gsub("-", "", as.character(Sys.Date()))

  spec_dir <- file.path(path.expand(OUTPUT_DIR), spec_label)
  dir.create(spec_dir, showWarnings = FALSE, recursive = TRUE)

  out_summary <- file.path(spec_dir, sprintf("%s_downsample_summary_%s.csv", spec_label, stamp))
  out_alpha   <- file.path(spec_dir, sprintf("%s_top_alphaCentric_%s.csv",   spec_label, stamp))
  out_beta    <- file.path(spec_dir, sprintf("%s_top_betaCentric_%s.csv",    spec_label, stamp))
  out_input   <- file.path(spec_dir, sprintf("%s_input_%s.csv",              spec_label, stamp))

  if (SAVE_SUMMARIES) {
    write_csv_safe(summary_row, out_summary)
    write_csv_safe(alpha_top,   out_alpha)
    write_csv_safe(beta_top,    out_beta)
    if (SAVE_INPUT_TABLE) write_csv_safe(df, out_input)
  }

  invisible(list(summary = summary_row, alpha = alpha_top, beta = beta_top))
}

invisible(lapply(plot_species, make_species_csvs))

## ---------- Replacement for TCRgraph: Odds Ratios + Forest Plot ----------
to_plot_global$type <- factor(to_plot_global$type,
                              levels = c("1TRA-multiTRB","multiTRA-1TRB"))

tab <- with(to_plot_global, table(species2, type))
if (!all(c("1TRA-multiTRB","multiTRA-1TRB") %in% colnames(tab))) {
  msg("One of the metaclone types is missing; skipping Forest plot. Saving counts only.")
  counts_csv <- file.path(path.expand(OUTPUT_DIR), "species_type_counts.csv")
  readr::write_csv(as.data.frame.matrix(tab), counts_csv)
} else {
  s1 <- sum(tab[, "1TRA-multiTRB"])
  s2 <- sum(tab[, "multiTRA-1TRB"])
  species_levels <- rownames(tab)

  or_df <- do.call(rbind, lapply(species_levels, function(sp){
    a1 <- tab[sp, "1TRA-multiTRB"]
    a2 <- tab[sp, "multiTRA-1TRB"]
    # 2x2: [this species vs others] x [α-centric vs β-centric]
    m <- matrix(c(a1, s1 - a1,
                  a2, s2 - a2),
                nrow = 2, byrow = FALSE,
                dimnames = list(c("alpha","beta"), c("this","other")))
    ft <- suppressWarnings(stats::fisher.test(m))
    data.frame(label = sp,
               p.value = unname(ft$p.value),
               OR      = unname(ft$estimate),
               lower   = unname(ft$conf.int[1]),
               upper   = unname(ft$conf.int[2]),
               stringsAsFactors = FALSE)
  }))

  levs <- SPECIES_SET
  or_df$label <- factor(or_df$label, levels = levs)

  plot_df <- or_df %>%
    dplyr::mutate(mean  = safe_log2(OR),
                  lower = safe_log2(lower),
                  upper = safe_log2(upper)) %>%
    dplyr::arrange(factor(label, levels = levs)) %>%
    dplyr::mutate(label = forcats::fct_rev(label))

  forest_theme <- ggplot2::theme_bw(32) +
    ggplot2::theme(panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   legend.position = "none",
                   panel.border = element_rect(colour = "black", fill = NA, size = 1))

  p_forest <- ggplot2::ggplot(plot_df, aes(x = label, y = mean, ymin = lower, ymax = upper)) +
      ggplot2::geom_pointrange(size = 1.2) +
      ggplot2::geom_hline(yintercept = 0, linetype = 2) +
      ggplot2::coord_flip() +
      ggplot2::labs(x = NULL, y = "\nlog2(odds ratio) +/- 95% CI") +
      forest_theme

  ## copy-safe date stamp (YYYYMMDD; no % tokens)
  ## ---------- Forest plot save (copy-safe, no sprintf) ----------
  # Date stamp as YYYYMMDD without using % formatting tokens
  stamp <- gsub('-', '', as.character(Sys.Date()))
  od <- path.expand(OUTPUT_DIR)

  forest_pdf <- file.path(od, paste0('Fig_1C_Forest_OR_', stamp, '.pdf'))
  forest_png <- file.path(od, paste0('Fig_1C_Forest_OR_', stamp, '.png'))
  or_csv     <- file.path(od, paste0('Forest_OR_table_', stamp, '.csv'))

  ggplot2::ggsave(forest_pdf, p_forest, width = 10, height = 8, units = 'in', dpi = 300)
  ggplot2::ggsave(forest_png, p_forest, width = 10, height = 8, units = 'in', dpi = 300)
  readr::write_csv(or_df %>% dplyr::arrange(factor(label, levels = levs)), or_csv)

  msg(paste('saved:', forest_pdf))
  msg(paste('saved:', forest_png))
  msg(paste('saved:', or_csv))
  
}


msg("Done.")
