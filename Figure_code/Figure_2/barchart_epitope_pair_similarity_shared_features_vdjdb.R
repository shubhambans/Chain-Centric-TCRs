# ---- Minimal VDJdb loader + multi-epitope plotting (end-to-end, CI cols) ----
rm(list = ls())

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
})

# If plyr is loaded, it masks count(); detach to be safe
if ("package:plyr" %in% search()) detach("package:plyr", unload = TRUE)

# =========================
# GLOBAL FONT / SIZE CONTROLS  (edit these only)
# =========================
FS_BASE        <- 26   # base theme text size
FS_AXIS_TEXT   <- 22
FS_AXIS_TITLE  <- 26
FS_STRIP_TEXT  <- 26
FS_BAR_LABEL   <- 10   # not used here (no bar labels), kept for convenience
FS_PLOT_W      <- 16   # inches
FS_PLOT_H      <- 10   # inches

# =========================
# TOGGLES
# =========================
CHAIN_TOGGLE      <- "TRA"          # "TRA" or "TRB"
METACLONE_TOGGLE  <- "alpha"        # "none","alpha","beta","alpha_only","beta_only"
FACET_PLOT        <- TRUE           # TRUE = facet per contrast; FALSE = single panel
VDJDB_PATH        <- "~/Desktop/Shubhams_paper_code/data/hu_with_metaclone_20250826_095644.csv"
ANCHOR            <- "SARS-CoV-2"   # anchor epitope
N_OTHERS          <- 4              # how many others to test vs anchor

# Restrict to the TOP10 antigen.species?
RESTRICT_TO_TOP10 <- TRUE
TOP10 <- c("CMV","InfluenzaA","EBV","SARS-CoV-2","HomoSapiens")

# =========================
# OUTPUTS  (absolute path)
# =========================
OUT_DIR <- file.path(getwd(), "Figure_Images")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# =========================
# THEME (single place to control typography)
# =========================
set_plot_theme <- function() {
  theme_set(
    theme_bw(base_size = FS_BASE) +
      theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.key       = element_blank(),
        axis.text        = element_text(size = FS_AXIS_TEXT),
        axis.title       = element_text(size = FS_AXIS_TITLE),
        strip.text       = element_text(size = FS_STRIP_TEXT, face = "bold")
      )
  )
}

# Prefer cairo for PDF text rendering
pick_pdf_device <- function(){
  # If Cairo isn't installed, ggsave will still work but may substitute fonts.
  if (requireNamespace("Cairo", quietly = TRUE)) {
    return(Cairo::CairoPDF)
  }
  # fallback
  return(cairo_pdf)
}

# =========================
# Helpers: case-insensitive column pick + smart CSV/TSV reader
# =========================
pick_ci <- function(df, candidates) {
  nms <- names(df); low <- tolower(nms)
  for (cand in candidates) {
    j <- which(low == tolower(cand))
    if (length(j)) return(nms[j[1]])
  }
  NA_character_
}

require_ci <- function(df, candidates, label = "column") {
  hit <- pick_ci(df, candidates)
  if (is.na(hit)) stop(sprintf(
    "Missing required %s. Tried: %s\nAvailable cols: %s",
    label, paste(candidates, collapse = ", "), paste(names(df), collapse = ", ")
  ))
  hit
}

smart_read <- function(path, col_types_chars = FALSE) {
  if (!is.character(path)) stop("smart_read expects a file path (character).")
  if (grepl("\\.csv$", path, ignore.case = TRUE)) {
    if (isTRUE(col_types_chars)) {
      readr::read_csv(path, show_col_types = FALSE,
                      col_types = readr::cols(.default = readr::col_character()))
    } else {
      readr::read_csv(path, show_col_types = FALSE)
    }
  } else {
    x <- if (isTRUE(col_types_chars)) {
      readr::read_tsv(path, comment = "#", show_col_types = FALSE,
                      col_types = readr::cols(.default = readr::col_character()))
    } else {
      readr::read_tsv(path, comment = "#", show_col_types = FALSE)
    }
    if (ncol(x) <= 1) {
      x <- if (isTRUE(col_types_chars)) {
        readr::read_csv(path, show_col_types = FALSE,
                        col_types = readr::cols(.default = readr::col_character()))
      } else {
        readr::read_csv(path, show_col_types = FALSE)
      }
    }
    x
  }
}

# =========================
# Loader: filter by chain; keep Homo sapiens; valid CDR3; return libid/chain/CDR3/V/J/epitope
# =========================
make_vdj_for_plots <- function(data_or_path,
                               chain = c("TRA","TRB"),
                               export_vdjdb = TRUE,
                               vdjdb_name = "vdjdb") {
  chain <- toupper(match.arg(chain))

  df_raw <- if (is.character(data_or_path)) {
    smart_read(data_or_path)  # auto CSV/TSV
  } else if (inherits(data_or_path, "data.frame")) {
    data_or_path
  } else stop("data_or_path must be a file path (character) or a data.frame/tibble.")

  # case-insensitive picks
  c_species <- require_ci(df_raw, c("species","Species","organism","Organism"), "species")
  c_libid   <- pick_ci(df_raw,    c("libid","complex.id"))
  c_cdr3    <- require_ci(df_raw, c("junction","junction_aa","cdr3","cdr3aa","cdr3_aa",
                                    "aa.cdr3","aaCDR3","CDR3","JUNCTION","JUNCTION_AA"), "CDR3/junction")
  c_v       <- require_ci(df_raw, c("v_gene","v.segm","V.GENE.and.allele","v.gene.and.allele"), "V gene")
  c_j       <- require_ci(df_raw, c("j_gene","j.segm","J.GENE.and.allele","j.gene.and.allele"), "J gene")
  c_epi     <- pick_ci(df_raw,    c("epitope","study_group","antigen.species","antigen_species"))
  c_chain   <- pick_ci(df_raw,    c("chain","Gene","gene","tcr.chain","tcr_chain"))

  if (export_vdjdb) assign(vdjdb_name, df_raw, envir = parent.frame())

  out <- df_raw %>%
    mutate(species_norm = stringr::str_replace_all(.data[[c_species]], "\\s+", "")) %>%
    filter(species_norm == "HomoSapiens") %>%
    { if (!is.na(c_chain)) dplyr::filter(., toupper(.data[[c_chain]]) == chain) else . } %>%
    transmute(
      libid    = if (!is.na(c_libid)) as.character(.data[[c_libid]]) else NA_character_,
      chain    = if (!is.na(c_chain)) toupper(.data[[c_chain]]) else chain,
      junction = as.character(.data[[c_cdr3]]),
      v_gene   = as.character(.data[[c_v]]),
      j_gene   = as.character(.data[[c_j]]),
      epitope  = if (!is.na(c_epi)) as.character(.data[[c_epi]]) else NA_character_
    ) %>%
    filter(
      !is.na(junction), junction != "",
      !is.na(v_gene),   v_gene   != "",
      !is.na(j_gene),   j_gene   != "",
      stringr::str_detect(junction, "^[Cc].*[FfWw]$")
    ) %>%
    distinct()

  out
}

# =========================
# Multi-epitope plotter (df already filtered)
# =========================
pairwise_germline_plot_multi <- function(
  df,
  epitopes = c("CMV","EBV","InfluenzaA"),
  anchor_epitope = "SARS-CoV-2",
  cdr3 = "junction", v = "v_gene", j = "j_gene", epi = "epitope",
  facet = TRUE
) {
  ep_all <- unique(c(if (!is.null(anchor_epitope)) anchor_epitope, epitopes))
  ep_all <- ep_all[!is.na(ep_all) & nzchar(ep_all)]
  stopifnot(length(ep_all) >= 2, length(ep_all) <= 7)

  class_levels <- c(
    "Same V, Same J",
    "Same V, Different J",
    "Different V, Same J",
    "Different V, Different J"
  )

  clean_df <- df %>%
    transmute(CDR3 = .data[[cdr3]], V = .data[[v]], J = .data[[j]], Ep = stringr::str_trim(.data[[epi]])) %>%
    separate_rows(Ep, sep = ",\\s*") %>%
    filter(Ep %in% ep_all,
           !is.na(CDR3) & CDR3 != "", !is.na(V) & V != "", !is.na(J) & J != "",
           stringr::str_detect(CDR3, "^[Cc].*[FfWw]$")) %>%
    distinct(CDR3, V, J, Ep) %>%
    as.data.frame()

  pairs_to_do <- if (is.null(anchor_epitope)) {
    t(utils::combn(unique(ep_all), 2))
  } else {
    others <- setdiff(unique(ep_all), anchor_epitope)
    if (!length(others)) stop("Only the anchor provided; add at least one other epitope.")
    cbind(anchor_epitope, others)
  }

  .add_ctx <- function(df0, e1, e2, key) {
    n <- nrow(df0)
    df0$epitope1   <- rep_len(e1,  n)
    df0$epitope2   <- rep_len(e2,  n)
    df0$pair_label <- rep_len(key, n)
    df0
  }

  do_one <- function(e1, e2) {
    clean_pair <- df %>%
      transmute(CDR3 = .data[[cdr3]], V = .data[[v]], J = .data[[j]], Ep = stringr::str_trim(.data[[epi]])) %>%
      separate_rows(Ep, sep = ",\\s*") %>%
      filter(Ep %in% c(e1, e2),
             !is.na(CDR3) & CDR3 != "", !is.na(V) & V != "", !is.na(J) & J != "",
             stringr::str_detect(CDR3, "^[Cc].*[FfWw]$")) %>%
      distinct(CDR3, V, J, Ep)

    shared_pair <- clean_pair %>%
      group_by(CDR3) %>% filter(dplyr::n_distinct(Ep) >= 2) %>% ungroup()

    if (nrow(shared_pair) == 0) {
      pairs <- tibble::tibble(CDR3=character(), V1=character(), J1=character(),
                              V2=character(), J2=character(), class=character())
      summary_pairs <- tibble::tibble(class = factor(class_levels, levels = class_levels),
                                      n = 0L, percent = 0)
    } else {
      pairs <- inner_join(
        shared_pair %>% filter(Ep == e1) %>% select(CDR3, V1 = V, J1 = J),
        shared_pair %>% filter(Ep == e2) %>% select(CDR3, V2 = V, J2 = J),
        by = "CDR3"
      ) %>%
        mutate(
          class = case_when(
            V1 == V2 & J1 == J2 ~ "Same V, Same J",
            V1 == V2 & J1 != J2 ~ "Same V, Different J",
            V1 != V2 & J1 == J2 ~ "Different V, Same J",
            TRUE                ~ "Different V, Different J"
          )
        ) %>%
        distinct(CDR3, V1, J1, V2, J2, .keep_all = TRUE)

      summary_pairs <- pairs %>%
        dplyr::count(class, name = "n") %>%
        tidyr::complete(class = class_levels, fill = list(n = 0)) %>%
        mutate(percent = if (sum(n) > 0) n / sum(n) else 0,
               class = factor(class, levels = class_levels))
    }

    key <- paste(e1, "vs", e2)
    shared_pair   <- .add_ctx(as.data.frame(shared_pair),   e1, e2, key)
    pairs         <- .add_ctx(as.data.frame(pairs),         e1, e2, key)
    summary_pairs <- .add_ctx(as.data.frame(summary_pairs), e1, e2, key)

    list(shared=shared_pair, pairs=pairs, summary=summary_pairs)
  }

  res_list <- apply(pairs_to_do, 1, function(r) do_one(r[1], r[2]))
  shared_all  <- dplyr::bind_rows(lapply(res_list, `[[`, "shared"))
  pairs_all   <- dplyr::bind_rows(lapply(res_list, `[[`, "pairs"))
  summary_all <- dplyr::bind_rows(lapply(res_list, `[[`, "summary"))

  pal <- c(
    "Same V, Same J"            = "#1f77b4",
    "Same V, Different J"       = "#ff7f0e",
    "Different V, Same J"       = "#2ca02c",
    "Different V, Different J"  = "#d62728"
  )

  if (facet) {
    p <- ggplot(summary_all, aes(class, percent, fill = class)) +
      geom_col(width = 0.6, show.legend = FALSE) +
      scale_x_discrete(limits = rev(class_levels)) +
      scale_fill_manual(values = pal, drop = FALSE) +
      scale_y_continuous(limits = c(0, 1),
                         breaks = c(0, 0.25, 0.5, 0.75, 1)) +
      labs(title = "", y = "\nFraction shared junctions", x = NULL) +
      coord_flip() +
      facet_wrap(~ pair_label, ncol = 2)
  } else {
    summary_total <- pairs_all %>%
      dplyr::count(class, name = "n") %>%
      tidyr::complete(class = class_levels, fill = list(n = 0)) %>%
      mutate(percent = if (sum(n) > 0) n / sum(n) else 0,
             class = factor(class, levels = class_levels))

    p <- ggplot(summary_total, aes(class, percent, fill = class)) +
      geom_col(width = 0.6, show.legend = FALSE) +
      scale_x_discrete(limits = rev(class_levels)) +
      scale_fill_manual(values = pal, drop = FALSE) +
      scale_y_continuous(limits = c(0, 1),
                         breaks = c(0, 0.25, 0.5, 0.75, 1)) +
      labs(title = "", y = "\nFraction shared junctions", x = NULL) +
      coord_flip()
  }

  invisible(list(plot = p,
                 clean = clean_df,
                 shared = shared_all,
                 pairs = pairs_all,
                 summary = summary_all))
}

# =========================
# VDJdb → metaclone classifier (case-insensitive cols)
# =========================
vdjdb_metaclones_only <- function(data_or_path, read_as_character = TRUE) {
  df_raw <- if (is.character(data_or_path)) {
    smart_read(data_or_path, col_types_chars = read_as_character)
  } else if (inherits(data_or_path, "data.frame")) {
    data_or_path
  } else stop("data_or_path must be a file path or a data.frame")

  df <- df_raw
  ci_rename <- function(df, target, candidates) {
    hit <- pick_ci(df, candidates)
    if (!is.na(hit) && !(target %in% names(df))) names(df)[names(df) == hit] <- target
    df
  }

  df <- ci_rename(df, "junction", c("junction","junction_aa","cdr3","cdr3aa","cdr3_aa","aa.cdr3","aaCDR3","CDR3","JUNCTION","JUNCTION_AA"))
  df <- ci_rename(df, "libid",    c("libid","complex.id"))
  df <- ci_rename(df, "v_gene",   c("v_gene","v.segm","V.GENE.and.allele","v.gene.and.allele"))
  df <- ci_rename(df, "j_gene",   c("j_gene","j.segm","J.GENE.and.allele","j.gene.and.allele"))
  df <- ci_rename(df, "chain",    c("chain","Gene","gene","tcr.chain","tcr_chain"))
  df <- ci_rename(df, "epitope",  c("epitope","study_group","antigen.species","antigen_species"))
  df <- ci_rename(df, "species",  c("species","Species","organism","Organism"))

  if (!("species" %in% names(df))) stop("Expected a 'species' column in VDJdb export")

  df <- df %>%
    mutate(species_norm = stringr::str_replace_all(.data$species, "\\s+", "")) %>%
    filter(species_norm == "HomoSapiens") %>%
    filter(!is.na(libid), libid != 0, libid != "0")

  if (!"chain" %in% names(df)) {
    if (!("v_gene" %in% names(df)) && !("j_gene" %in% names(df)))
      stop("Need either 'chain' or V/J gene columns to infer chain.")
    df <- df %>%
      mutate(chain = dplyr::case_when(
        stringr::str_detect(coalesce(v_gene, ""), "^TRAV|^TRAJ") |
          stringr::str_detect(coalesce(j_gene, ""), "^TRAV|^TRAJ") ~ "TRA",
        stringr::str_detect(coalesce(v_gene, ""), "^TRBV|^TRBJ") |
          stringr::str_detect(coalesce(j_gene, ""), "^TRBV|^TRBJ") ~ "TRB",
        TRUE ~ NA_character_
      ))
  }

  df <- df %>%
    mutate(chain = toupper(chain)) %>%
    filter(chain %in% c("TRA","TRB"))

  df_clean <- df %>%
    filter(!is.na(junction), junction != "") %>%
    filter(stringr::str_detect(junction, "^[Cc].*[FfWw]$")) %>%
    distinct(libid, chain, junction, v_gene, j_gene, .keep_all = TRUE)

  both_ids <- df_clean %>%
    distinct(libid, chain) %>% count(libid) %>% filter(n >= 2) %>% pull(libid)
  df_clean <- df_clean %>% filter(libid %in% both_ids)

  tra_all <- df_clean %>% filter(chain == "TRA") %>%
    transmute(libid, junction.tra = junction)
  trb_all <- df_clean %>% filter(chain == "TRB") %>%
    transmute(libid, junction.trb = junction)
  paired_all <- dplyr::inner_join(tra_all, trb_all, by = "libid")

  alpha_ids <- paired_all %>% group_by(junction.tra) %>%
    filter(dplyr::n_distinct(junction.trb) > 1) %>% ungroup() %>% distinct(libid) %>% pull(libid)
  beta_ids  <- paired_all %>% group_by(junction.trb) %>%
    filter(dplyr::n_distinct(junction.tra) > 1) %>% ungroup() %>% distinct(libid) %>% pull(libid)

  metaclones <- tibble::tibble(
    libid = c(alpha_ids, beta_ids),
    type  = c(rep("alphaCentric", length(alpha_ids)), rep("betaCentric", length(beta_ids)))
  ) %>% distinct(libid, type)

  list(vdjdb_metaclone_rows = df_raw, metaclones = metaclones)
}

# =========================
# Apply metaclone selection to plotting df (by libid)
# =========================
apply_metaclone_selection <- function(df_plot, metaclones, selection = "none") {
  selection <- match.arg(selection, c("none","alpha","beta","alpha_only","beta_only"))
  if (selection == "none" || !"libid" %in% names(df_plot)) return(df_plot)

  alpha_ids <- metaclones %>% filter(type == "alphaCentric") %>% pull(libid)
  beta_ids  <- metaclones %>% filter(type == "betaCentric")  %>% pull(libid)
  keep <- switch(selection,
                 alpha      = alpha_ids,
                 beta       = beta_ids,
                 alpha_only = setdiff(alpha_ids, intersect(alpha_ids, beta_ids)),
                 beta_only  = setdiff(beta_ids,  intersect(alpha_ids, beta_ids)))
  df_plot %>% semi_join(tibble(libid = keep), by = "libid")
}

# =========================
# END-TO-END
# =========================
set_plot_theme()

# 1) Compute metaclones from your CSV/TSV
res_meta <- vdjdb_metaclones_only(VDJDB_PATH, read_as_character = TRUE)

# 2) Load chain-specific df (libid + chain + junction/v/j/epitope)
df_chain <- make_vdj_for_plots(res_meta$vdjdb_metaclone_rows, chain = CHAIN_TOGGLE)

# 3) Optional: restrict to TOP10 antigen.species
if (RESTRICT_TO_TOP10) {
  df_chain <- df_chain %>% filter(epitope %in% TOP10)
}

# 4) Apply metaclone filter toggle
df_chain <- apply_metaclone_selection(df_chain, res_meta$metaclones, METACLONE_TOGGLE)

# 5) Choose “others” epitopes present vs anchor
anchor <- ANCHOR
others <- df_chain %>%
  filter(!is.na(epitope)) %>%
  mutate(epi = as.character(epitope)) %>%
  filter(epi != anchor) %>%
  count(epi, sort = TRUE) %>%
  pull(epi) %>%
  { if (RESTRICT_TO_TOP10) intersect(., TOP10) else . } %>%
  head(N_OTHERS)

if (!length(others)) stop("No 'other' epitopes found after filters; adjust TOP10/CHAIN/METACLONE toggles.")

# 6) Plot
res <- pairwise_germline_plot_multi(
  df_chain,
  epitopes       = others,
  anchor_epitope = anchor,
  facet          = FACET_PLOT
)

print(res$plot)

# 7) Save plot to OUT_DIR (THIS is the corrected save)
out_pdf <- file.path(OUT_DIR, sprintf("Fig2D_S10_%s_%s_%s.pdf",
                                      CHAIN_TOGGLE, METACLONE_TOGGLE,
                                      gsub("[^A-Za-z0-9-]+","_", ANCHOR)))

cat("OUT_DIR:", normalizePath(OUT_DIR, mustWork = FALSE), "\n")
cat("Saving to:", out_pdf, "\n")
stopifnot(dir.exists(OUT_DIR))

pdf_dev <- pick_pdf_device()

tryCatch({
  ggsave(filename = out_pdf, plot = res$plot,
         width = FS_PLOT_W, height = FS_PLOT_H, units = "in",
         device = pdf_dev, useDingbats = FALSE)
  cat("Saved OK:", out_pdf, "\n")
}, error = function(e) {
  message("ggsave FAILED: ", conditionMessage(e))
})

print(list.files(OUT_DIR, pattern = "Fig2D_S10", full.names = TRUE))
