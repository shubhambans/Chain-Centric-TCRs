## =========================
## End-to-end: Epitope Overlap Heatmap + Significance
## + explicit OUT_DIR handling (all outputs written there)
## + global min-set normalization across TRA & TRB (optional)
## + explicit p-value report for TRA: InfluenzaA vs SARS-CoV-2
## + export tidy p-value tables (post-normalization)
## =========================

rm(list = ls())

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(readr)
  library(reshape2); library(tibble); library(scales); library(stringr); library(rlang)
})

## ---- Messaging / utils ---------------------------------------------------
msg <- function(...) message(sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), paste0(...)))
`%||%` <- function(x, y) if (is.null(x)) y else x
sanitize_chr <- function(x) {
  if (is.list(x)) vapply(x, function(el) paste(as.character(unlist(el, use.names = FALSE)), collapse = ", "), character(1))
  else as.character(x)
}

## ---- CONFIG --------------------------------------------------------------
vdj_path    <- "~/Desktop/Shubhams_paper_code/data/hu_with_metaclone_20250826_095644.csv"

# >>> OUTPUT DIRECTORY (NEW) <<<
OUT_DIR <- file.path(getwd(), "Figure_Images")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# All outputs will go here:
OUT_pdf_HEATMAP <- file.path(OUT_DIR, "epitope_overlap_TRA_TRB.pdf")
OUT_CSV_ALL     <- file.path(OUT_DIR, "overlap_pvalues_downsampled.csv")
OUT_CSV_SIG     <- file.path(OUT_DIR, "overlap_pvalues_downsampled_sig.csv")
OUT_TXT_PREFIX  <- file.path(OUT_DIR, "Overlap") # used to build list filenames

CHAIN_SET   <- c("TRA","TRB")
KEY_MODE    <- "cdr3"                 # "cdr3","cdr3_v","cdr3_j","cdr3_v_j"
METACLONE   <- "alpha"                # "all","none","alpha","beta","both"
CDR3_REGEX  <- NULL
REQUIRE_CD4 <- FALSE

FILL_MODE   <- "jaccard"              # "count" or "jaccard"
SHOW_COUNTS <- TRUE
SHOW_STARS  <- TRUE

# Significance method
PVAL_METHOD <- "perm"                 # "hypergeom" (one-sided), "fisher_two_sided", "perm"
N_PERM      <- 5000
set.seed(12345)

# Multiple-testing correction (applied AFTER merging chains)
APPLY_CORRECTION        <- TRUE
CORRECTION_METHOD       <- "BH"       # "BH", "holm", "bonferroni"
ADJUST_SCOPE            <- "global"   # "global" or "per_chain"
EXCLUDE_DIAGONAL_IN_ADJ <- TRUE

# Star logic
STAR_LEVELS <- c(`1e-4`="****", `1e-3`="***", `1e-2`="**", `5e-2`="*")
USE_FLOORS  <- TRUE
MIN_INTER   <- 2
MIN_JAC     <- 0.01

# Epitopes to include/order
TOP_EPI <- c("CMV","InfluenzaA","EBV","SARS-CoV-2","HomoSapiens")

# Reporting target
REPORT_CHAIN <- "TRA"
REPORT_EPI_A <- "InfluenzaA"
REPORT_EPI_B <- "SARS-CoV-2"

# ---- Normalization --------------------------------------------------------
NORMALIZE_TO_MIN_SET <- TRUE          # downsample each epitope’s key set to min size
NORMALIZE_SCOPE <- "global"           # "per_chain" or "global"

## ---- Canon helpers -------------------------------------------------------
norm_token <- function(x){
  x <- sanitize_chr(x)
  x <- gsub("[\u2010-\u2015\u2212\uFE63\uFF0D]", "-", x, perl=TRUE)
  x <- tolower(trimws(x))
  gsub("[^a-z0-9-]+","",x)
}
canon_epitope <- function(x){
  tok <- norm_token(x)
  map_one <- function(t){
    if (t %in% c("cmv","hcmv","humancytomegalovirus","cytomegalovirus")) return("CMV")
    if (t %in% c("influenzaa","influenzaavirus","iav","influenza-a"))     return("InfluenzaA")
    if (t %in% c("influenzab","influenza-b"))                             return("InfluenzaB")
    if (t %in% c("ebv","epsteinbarr","epsteinbarrvirus","hhv4"))          return("EBV")
    if (t %in% c("sars-cov-2","sarscov2","covid19","covid-19","2019-ncov","ncov2019","betacoronavirus")) return("SARS-CoV-2")
    if (t %in% c("homosapiens","homo-sapiens","human"))                   return("HomoSapiens")
    if (t %in% c("hiv","hiv1","hiv-1"))                                   return("HIV")
    if (t %in% c("hcv","hepatitisc","hepatitiscvirus"))                   return("HCV")
    if (t %in% c("mtb","mycobacteriumtuberculosis"))                      return("Mtb")
    if (t %in% c("yfv","yellowfever","yellowfevervirus"))                 return("YFV")
    NA_character_
  }
  vapply(tok, map_one, character(1))
}
.read_any_table <- function(path, guess_max=2e5){
  df <- tryCatch(readr::read_csv(path, show_col_types=FALSE, progress=FALSE, guess_max=guess_max),
                 error=function(e) NULL)
  if (is.null(df) || ncol(df)==1) df <- readr::read_tsv(path, show_col_types=FALSE, progress=FALSE, guess_max=guess_max)
  df
}

## ---- Standardize ---------------------------------------------------------
standardize_vdjdb <- function(data_or_path){
  df_raw <- if (is.character(data_or_path)) .read_any_table(data_or_path)
  else if (inherits(data_or_path, "data.frame")) data_or_path
  else stop("Provide a CSV/TSV path or a data.frame")
  df_raw <- tibble::as_tibble(df_raw); names(df_raw) <- tolower(trimws(names(df_raw)))

  pick1 <- function(cands, optional=FALSE, msg=NULL){
    m <- intersect(cands, names(df_raw))
    if (length(m)) m[1] else if (optional) NULL else stop(msg %||% paste("Missing:", paste(cands, collapse=", ")))
  }

  c_species <- pick1(c("species","organism","species.taxon"), TRUE)
  c_cdr3    <- pick1(c("junction","cdr3","aa.cdr3","aacdr3"), msg="Need junction/CDR3 column")
  c_v       <- pick1(c("v_gene","v.segm","v.seg","v"))
  c_j       <- pick1(c("j_gene","j.segm","j.seg","j"))
  c_epi     <- pick1(c("epitope"), TRUE)
  c_antsp   <- pick1(c("antigen.species","antigen.epitope","study_group","antigen.epitope.name"), TRUE)
  c_chain   <- pick1(c("chain","gene","tcr.chain","gene.chain","gene_name"), TRUE)
  c_libid   <- pick1(c("complex.id","libid","complex_id","complexid","pair.id"), TRUE)
  c_cell    <- pick1(c("cell_type","cell.type","celltype"), TRUE)

  ep_col <- if (!is.null(c_epi)) c_epi else if (!is.null(c_antsp)) c_antsp else stop("Need `epitope` or `antigen.species`")

  tibble::tibble(
    species     = sanitize_chr(if (!is.null(c_species)) df_raw[[c_species]] else NA),
    junction    = sanitize_chr(df_raw[[c_cdr3]]),
    v_gene      = sanitize_chr(df_raw[[c_v]]),
    j_gene      = sanitize_chr(df_raw[[c_j]]),
    epitope_raw = sanitize_chr(df_raw[[ep_col]]),
    chain_in    = sanitize_chr(if (!is.null(c_chain)) df_raw[[c_chain]] else NA),
    libid       = sanitize_chr(if (!is.null(c_libid)) df_raw[[c_libid]] else NA),
    cell_type   = sanitize_chr(if (!is.null(c_cell)) df_raw[[c_cell]] else NA)
  )
}

## ---- Metaclone filtering -------------------------------------------------
subset_by_metaclone <- function(df_std, cdr3_regex=NULL, require_cd4_only=FALSE,
                                metaclone=c("all","none","alpha","beta","both"),
                                require_v=FALSE, require_j=FALSE){
  metaclone <- match.arg(metaclone)

  if (!all(is.na(df_std$species))) {
    sp <- gsub("\\s+","", df_std$species)
    df_std <- df_std[sp %in% c("HomoSapiens","Homosapiens","HomoSapiens"), , drop=FALSE]
  }
  if (require_cd4_only && !all(is.na(df_std$cell_type))) {
    df_std <- df_std[grepl("CD4", df_std$cell_type, ignore.case=TRUE), , drop=FALSE]
  }

  chain_inf <- ifelse(!is.na(df_std$chain_in) & nzchar(df_std$chain_in), toupper(df_std$chain_in), NA_character_)
  chain_inf <- ifelse(!is.na(chain_inf), chain_inf,
                      ifelse(grepl("^trav|^traj", df_std$v_gene, TRUE) | grepl("^trav|^traj", df_std$j_gene, TRUE), "TRA",
                             ifelse(grepl("^trbv|^trbj", df_std$v_gene, TRUE) | grepl("^trbj|^trbv", df_std$j_gene, TRUE), "TRB", NA_character_)))
  df_std$chain <- chain_inf

  ok <- !is.na(df_std$junction) & nzchar(df_std$junction) &
    df_std$chain %in% c("TRA","TRB") &
    (if (require_v) !is.na(df_std$v_gene) & nzchar(df_std$v_gene) else TRUE) &
    (if (require_j) !is.na(df_std$j_gene) & nzchar(df_std$j_gene) else TRUE)
  df <- df_std[ok, , drop=FALSE]

  if (!is.null(cdr3_regex)) df <- df[grepl(cdr3_regex, df$junction), , drop=FALSE]
  if (metaclone %in% c("all","none")) return(tibble::as_tibble(df))

  if (all(is.na(df$libid) | !nzchar(df$libid))) {
    warning("Metaclone filter requested but no 'libid' present; returning unfiltered.")
    return(tibble::as_tibble(df))
  }

  df_clean <- df %>% dplyr::distinct(libid, chain, junction, v_gene, j_gene, epitope_raw, .keep_all=TRUE)

  # robust “has both chains”
  lib_has_TRA <- df_clean %>% filter(chain=="TRA") %>% distinct(libid)
  lib_has_TRB <- df_clean %>% filter(chain=="TRB") %>% distinct(libid)
  both_ids <- intersect(lib_has_TRA$libid, lib_has_TRB$libid)

  df_pair  <- df_clean %>% dplyr::filter(libid %in% both_ids)

  tra <- df_pair %>% dplyr::filter(chain=="TRA") %>% dplyr::transmute(libid, junction.tra=junction)
  trb <- df_pair %>% dplyr::filter(chain=="TRB") %>% dplyr::transmute(libid, junction.trb=junction)
  paired <- dplyr::inner_join(tra, trb, by="libid")

  alpha_meta_ids <- paired %>% dplyr::group_by(junction.tra) %>%
    dplyr::filter(dplyr::n_distinct(junction.trb)>1) %>%
    dplyr::ungroup() %>% dplyr::distinct(libid) %>% dplyr::pull(libid)

  beta_meta_ids  <- paired %>% dplyr::group_by(junction.trb) %>%
    dplyr::filter(dplyr::n_distinct(junction.tra)>1) %>%
    dplyr::ungroup() %>% dplyr::distinct(libid) %>% dplyr::pull(libid)

  keep_ids <- switch(metaclone,
                     alpha = alpha_meta_ids,
                     beta  = beta_meta_ids,
                     both  = union(alpha_meta_ids, beta_meta_ids))
  tibble::as_tibble(df_pair %>% dplyr::filter(libid %in% keep_ids))
}

## ---- Per-chain tidy frame -----------------------------------------------
make_chain_df <- function(df_std, chain=c("TRA","TRB"),
                          key_mode=c("cdr3","cdr3_v","cdr3_j","cdr3_v_j"),
                          cdr3_regex=NULL, require_cd4_only=FALSE,
                          metaclone=c("all","none","alpha","beta","both")){
  chain <- match.arg(chain); key_mode <- match.arg(key_mode); metaclone <- match.arg(metaclone)
  need_v <- key_mode %in% c("cdr3_v","cdr3_v_j"); need_j <- key_mode %in% c("cdr3_j","cdr3_v_j")

  df_filt <- subset_by_metaclone(df_std, cdr3_regex, require_cd4_only, metaclone, need_v, need_j)

  df_out <- df_filt %>%
    dplyr::mutate(epitope_raw = sanitize_chr(epitope_raw),
                  ep_list     = strsplit(epitope_raw, "[,;]\\s*")) %>%
    tidyr::unnest_longer(ep_list, values_to="ep_one") %>%
    dplyr::mutate(ep_canon = canon_epitope(ep_one),
                  v_gene   = sanitize_chr(v_gene),
                  j_gene   = sanitize_chr(j_gene),
                  junction = sanitize_chr(junction)) %>%
    dplyr::filter(chain==!!chain, !is.na(ep_canon), ep_canon %in% TOP_EPI)

  tibble::as_tibble(df_out %>% dplyr::select(chain, v_gene, j_gene, junction, ep_canon, libid))
}

explode_keys_epitopes <- function(df_chain, key_mode=c("cdr3","cdr3_v","cdr3_j","cdr3_v_j")){
  key_mode <- match.arg(key_mode)
  df_chain <- tibble::as_tibble(df_chain) %>%
    dplyr::mutate(junction=sanitize_chr(junction),
                  v_gene=sanitize_chr(v_gene),
                  j_gene=sanitize_chr(j_gene),
                  ep_canon=sanitize_chr(ep_canon))
  KEY <- switch(key_mode,
                cdr3     = df_chain$junction,
                cdr3_v   = paste(df_chain$v_gene, df_chain$junction, sep="_"),
                cdr3_j   = paste(df_chain$j_gene, df_chain$junction, sep="_"),
                cdr3_v_j = paste(df_chain$v_gene, df_chain$j_gene, df_chain$junction, sep="_")
  )
  tibble::tibble(KEY=KEY, EPI=df_chain$ep_canon) %>% dplyr::distinct(KEY, EPI)
}

## ---- Build incidence matrix (with optional normalization) ----------------
build_bin_matrix <- function(KE, epis, normalize = FALSE, fixed_k = NULL) {
  epis_present <- intersect(epis, unique(KE$EPI))
  if (!length(epis_present)) {
    return(list(bin = matrix(FALSE, nrow = 0, ncol = 0), epi_ids = character(0)))
  }
  KE <- KE %>% dplyr::filter(EPI %in% epis_present)

  bin_tbl <- table(KE$KEY, KE$EPI)
  # ensure full epi column set
  for (ep in epis_present) {
    if (!ep %in% colnames(bin_tbl)) {
      bin_tbl <- cbind(bin_tbl, setNames(list(integer(nrow(bin_tbl))), ep))
    }
  }
  bin_tbl <- bin_tbl[, epis_present, drop = FALSE]
  bin <- bin_tbl > 0

  if (!normalize) return(list(bin = bin, epi_ids = colnames(bin)))

  col_sums <- colSums(bin)
  pos_cols <- col_sums[col_sums > 0]
  if (!length(pos_cols)) return(list(bin = bin, epi_ids = colnames(bin)))
  min_k <- if (!is.null(fixed_k)) fixed_k else min(pos_cols)

  keep_rows <- rep(FALSE, nrow(bin))
  sel_idx_list <- vector("list", ncol(bin))
  for (j in seq_len(ncol(bin))) {
    idx <- which(bin[, j])
    if (length(idx) > 0) {
      if (length(idx) > min_k) idx <- sample(idx, min_k)
      sel_idx_list[[j]] <- idx
      keep_rows[idx] <- TRUE
    } else {
      sel_idx_list[[j]] <- integer(0)
    }
  }
  if (!any(keep_rows)) return(list(bin = bin[FALSE, , drop = FALSE], epi_ids = colnames(bin)))
  bin2 <- bin[keep_rows, , drop = FALSE]
  for (j in seq_len(ncol(bin2))) {
    sel <- sel_idx_list[[j]]
    if (length(sel)) {
      new_rows <- match(sel, which(keep_rows))
      bin2[, j] <- FALSE
      bin2[new_rows, j] <- TRUE
    } else {
      bin2[, j] <- FALSE
    }
  }
  list(bin = bin2, epi_ids = colnames(bin2), min_k = min_k)
}

## ---- Global min across TRA & TRB -----------------------------------------
compute_global_min_k <- function(dfs, key_mode, epis) {
  mins <- c()
  for (chain in names(dfs)) {
    df_chain <- dfs[[chain]]
    if (is.null(df_chain) || !nrow(df_chain)) next
    KE <- explode_keys_epitopes(df_chain, key_mode) %>% dplyr::filter(EPI %in% epis)
    if (!nrow(KE)) next
    bin_tbl <- table(KE$KEY, KE$EPI)
    if (!ncol(bin_tbl)) next
    col_sums <- colSums(bin_tbl > 0)
    pos <- col_sums[col_sums > 0]
    if (length(pos)) mins <- c(mins, min(pos))
  }
  if (length(mins)) min(mins) else NA_integer_
}

## ---- P-values ------------------------------------------------------------
overlap_p_hypergeom <- function(a, A, B, N) stats::phyper(q=a-1, m=A, n=N-A, k=B, lower.tail=FALSE)
fisher_two_sided <- function(a, A, B, N){
  b <- A - a; c <- B - a; d <- N - (a + b + c)
  if (min(a,b,c,d) < 0) return(NA_real_)
  stats::fisher.test(matrix(c(a,b,c,d), 2, 2), alternative="two.sided")$p.value
}
permute_overlap_p <- function(a, A, B, N, bin_mat, col_i, col_j, nperm=1000){
  if (A == 0 || B == 0) return(NA_real_)
  obs <- a
  perm_vals <- replicate(nperm, {
    perm_idx <- sample(N)
    sum(bin_mat[perm_idx, col_i] & bin_mat[, col_j])
  })
  (sum(perm_vals >= obs) + 1) / (nperm + 1)
}
p_to_stars <- function(x, thresholds=STAR_LEVELS){
  out <- rep("", length(x)); thr <- as.numeric(names(thresholds)); lab <- unname(unlist(thresholds))
  for (i in seq_along(thr)) out[x < thr[i]] <- lab[i]
  out
}

## ---- Overlap computation per chain (raw p only) --------------------------
compute_epitope_overlap <- function(df_chain, chain_label,
                                    key_mode=c("cdr3","cdr3_v","cdr3_j","cdr3_v_j"),
                                    fixed_k=NULL){
  key_mode <- match.arg(key_mode)
  KE <- explode_keys_epitopes(df_chain, key_mode)
  if (!nrow(KE)) return(tibble(E1=character(),E2=character(),Inter=integer(),Union=integer(),Jac=numeric(),
                               A_size=integer(),B_size=integer(),N_universe=integer(),p=numeric(),Chain=character()))
  epis_present <- intersect(TOP_EPI, unique(KE$EPI))
  if (!length(epis_present)) return(tibble(E1=character(),E2=character(),Inter=integer(),Union=integer(),Jac=numeric(),
                                           A_size=integer(),B_size=integer(),N_universe=integer(),p=numeric(),Chain=character()))
  KE <- KE %>% dplyr::filter(EPI %in% epis_present)

  bm <- build_bin_matrix(
    KE, epis_present,
    normalize = NORMALIZE_TO_MIN_SET,
    fixed_k   = if (isTRUE(NORMALIZE_TO_MIN_SET) && NORMALIZE_SCOPE=="global") fixed_k else NULL
  )
  bin <- bm$bin
  if (!length(bin)) return(tibble(E1=character(),E2=character(),Inter=integer(),Union=integer(),Jac=numeric(),
                                  A_size=integer(),B_size=integer(),N_universe=integer(),p=numeric(),Chain=character()))
  epi_ids <- colnames(bin); n <- length(epi_ids); N_univ <- nrow(bin)

  Inter <- Union <- Jac <- matrix(0, n, n, dimnames=list(epi_ids, epi_ids))
  A_sz <- B_sz <- Pval <- matrix(NA_real_, n, n, dimnames=list(epi_ids, epi_ids))

  for (i in seq_len(n)) for (j in seq_len(i)) {
    a_vec <- bin[, i]; b_vec <- bin[, j]
    a <- sum(a_vec & b_vec); A <- sum(a_vec); B <- sum(b_vec); u <- sum(a_vec | b_vec)
    jac <- if (u==0) NA_real_ else a/u
    p <- NA_real_
    if (A > 0 && B > 0) {
      if (PVAL_METHOD == "hypergeom")             p <- overlap_p_hypergeom(a, A, B, N_univ)
      else if (PVAL_METHOD == "fisher_two_sided") p <- fisher_two_sided(a, A, B, N_univ)
      else if (PVAL_METHOD == "perm")             p <- permute_overlap_p(a, A, B, N_univ, bin, i, j, N_PERM)
      else stop("Unknown PVAL_METHOD")
    }
    Inter[i,j] <- a; Union[i,j] <- u; Jac[i,j] <- jac; A_sz[i,j] <- A; B_sz[i,j] <- B; Pval[i,j] <- p
  }

  reshape2::melt(Inter, varnames=c("E1","E2"), value.name="Inter") %>%
    dplyr::mutate(E1i=match(E1, epi_ids), E2i=match(E2, epi_ids),
                  Union=as.vector(Union[cbind(E1i,E2i)]),
                  Jac  =as.vector(Jac[cbind(E1i,E2i)]),
                  A_size=as.vector(A_sz[cbind(E1i,E2i)]),
                  B_size=as.vector(B_sz[cbind(E1i,E2i)]),
                  N_universe=N_univ,
                  p=as.vector(Pval[cbind(E1i,E2i)])) %>%
    dplyr::filter(E1i >= E2i) %>%
    dplyr::mutate(Chain = sanitize_chr(chain_label)) %>%
    dplyr::select(E1,E2,Inter,Union,Jac,A_size,B_size,N_universe,p,Chain) %>%
    tibble::as_tibble()
}

## ---- Plot ----------------------------------------------------------------
plot_epitope_heatmap_one_triangle <- function(df_in, fill_mode=c("count","jaccard")){
  fill_mode <- match.arg(fill_mode); stopifnot(nrow(df_in)>0)

  present <- sort(unique(c(df_in$E1, df_in$E2)))
  lvl     <- intersect(TOP_EPI, present)
  df_in$E1 <- factor(df_in$E1, levels=lvl)
  df_in$E2 <- factor(df_in$E2, levels=lvl)
  df_in <- df_in %>% dplyr::mutate(Diag=(E1==E2), stars=ifelse(is.na(stars),"",stars))

  idx <- data.frame(E1i=as.integer(df_in$E1), E2i=as.integer(df_in$E2))
  df_in <- df_in[idx$E1i >= idx$E2i, , drop=FALSE]

  if (fill_mode=="count"){
    non_diag_max <- suppressWarnings(max(df_in$Inter[!df_in$Diag], na.rm=TRUE)); if (!is.finite(non_diag_max)) non_diag_max <- 1
    fill_scale <- scale_fill_gradientn(colours=c("white","#fddbc7","#ef8a62","#b2182b"),
                                       values=if (non_diag_max>1) scales::rescale(c(0,1,5,non_diag_max)) else c(0,.33,.66,1),
                                       na.value="white", name="Shared keys")
    fill_var <- "Inter"
  } else {
    non_diag_max <- suppressWarnings(max(df_in$Jac[!df_in$Diag], na.rm=TRUE)); if (!is.finite(non_diag_max)) non_diag_max <- 1
    fill_scale <- scale_fill_gradientn(colours=c("white","red","darkred"),
                                       values=c(0,.5,1), limits=c(0,non_diag_max),
                                       na.value="white", breaks=c(0,non_diag_max/2,non_diag_max),
                                       labels=number(c(0,non_diag_max/2,non_diag_max), accuracy=0.01),
                                       name="Jaccard")
    fill_var <- "Jac"
  }

  p <- ggplot(df_in, aes(E1,E2)) +
    geom_tile(fill="white", color="black") +
    geom_tile(data=df_in %>% dplyr::filter(Diag), fill="grey75", color="black") +
    geom_text(data=df_in %>% dplyr::filter(Diag), label="—", size=10, fontface="italic") +
    geom_tile(data=df_in %>% dplyr::filter(!Diag), aes_string(fill=fill_var), color="black") +
    fill_scale

  # Diagonal label: per-epitope size (post-normalization)
  p <- p +
    geom_tile(data = df_in %>% dplyr::filter(Diag),
              fill = "grey75", color = "black") +
    geom_text(data = df_in %>% dplyr::filter(Diag),
              aes(label = A_size), size = 7)

  if (SHOW_COUNTS) p <- p + geom_text(data=df_in %>% dplyr::filter(!Diag),
                                      aes(label=ifelse(is.na(!!sym(fill_var)),"",as.character(Inter))),
                                      size=7, vjust=1.2)
  if (SHOW_STARS)  p <- p + geom_text(data=df_in %>% dplyr::filter(!Diag),
                                      aes(label=stars), size=10, vjust=-0.2)

  p + coord_fixed() + facet_wrap(~ Chain, ncol=2) +
    theme_minimal(base_size=14) +
    theme(axis.text.x=element_text(angle=70, hjust=1, size=22),
          axis.text.y=element_text(size=22),
          strip.text=element_text(size=28, face="bold"),
          panel.grid=element_blank(),
          axis.title=element_blank()) +
    scale_x_discrete(drop=FALSE) + scale_y_discrete(drop=FALSE)
}

## ---- Optional: write overlap list files ----------------------------------
write_overlap_lists <- function(df_chain, chain, key_mode="cdr3",
                                epiA="InfluenzaA", epiB="SARS-CoV-2",
                                out_dir=OUT_DIR){
  KE <- explode_keys_epitopes(df_chain, key_mode)
  A  <- unique(KE$KEY[KE$EPI == epiA]); B <- unique(KE$KEY[KE$EPI == epiB])
  inter_keys <- sort(intersect(A, B))
  get_cdr3 <- function(keys) if (key_mode=="cdr3") keys else sub(".*_", "", keys)
  inter_cdr3 <- sort(unique(get_cdr3(inter_keys)))
  fn <- file.path(out_dir, sprintf("Overlap_%s_vs_%s_%s.txt", epiA, epiB, chain))
  writeLines(inter_cdr3, con = fn)
  msg(sprintf("%s: wrote %d overlapping CDR3s to %s", chain, length(inter_cdr3), normalizePath(fn)))
}

## =============================== RUN =====================================
msg("Loading VDJdb-like file...")
hu <- standardize_vdjdb(vdj_path)

dfs <- list()
if ("TRA" %in% CHAIN_SET) dfs$TRA <- make_chain_df(hu, "TRA", KEY_MODE, CDR3_REGEX, REQUIRE_CD4, METACLONE)
if ("TRB" %in% CHAIN_SET) dfs$TRB <- make_chain_df(hu, "TRB", KEY_MODE, CDR3_REGEX, REQUIRE_CD4, METACLONE)

# Quick sanity table
summ_keys <- function(df, chain){
  if (is.null(df) || !nrow(df)) return(tibble(EPI=character(), A_or_B=integer(), Chain=character()))
  KE <- explode_keys_epitopes(df, KEY_MODE) %>% dplyr::filter(EPI %in% TOP_EPI)
  KE %>% dplyr::count(EPI, name="A_or_B") %>%
    dplyr::arrange(dplyr::desc(A_or_B)) %>%
    dplyr::mutate(Chain = chain)
}
print(rbind(summ_keys(dfs$TRA,"TRA"), summ_keys(dfs$TRB,"TRB")))

# Compute the global min set size across chains if requested
global_min_k <- NA_integer_
if (isTRUE(NORMALIZE_TO_MIN_SET) && NORMALIZE_SCOPE=="global") {
  global_min_k <- compute_global_min_k(dfs, KEY_MODE, TOP_EPI)
  if (is.na(global_min_k)) {
    msg("Global min set size could not be computed (no positive columns).")
  } else {
    msg(sprintf("Global normalization enabled. Using min_k = %d across TRA & TRB.", global_min_k))
  }
}

check_pair <- function(df_chain, key_mode, tag){
  if (is.null(df_chain) || !nrow(df_chain)) return(invisible(NULL))
  KE <- explode_keys_epitopes(df_chain, key_mode)
  infl <- unique(KE$KEY[KE$EPI=="InfluenzaA"])
  cov2 <- unique(KE$KEY[KE$EPI=="SARS-CoV-2"])
  msg(sprintf("[%s] InfluenzaA=%d, SARS-CoV-2=%d, INTER=%d",
              tag, length(infl), length(cov2), length(intersect(infl, cov2))))
}
msg(sprintf("KEY_MODE=%s  METACLONE=%s  PVAL_METHOD=%s  N_PERM=%d  NORMALIZE=%s (%s)",
            KEY_MODE, METACLONE, PVAL_METHOD, N_PERM,
            as.character(NORMALIZE_TO_MIN_SET), NORMALIZE_SCOPE))
if (!is.null(dfs$TRA)) check_pair(dfs$TRA, KEY_MODE, "TRA")
if (!is.null(dfs$TRB)) check_pair(dfs$TRB, KEY_MODE, "TRB")

# Compute raw p's per chain (using the same global min_k if applicable)
out_list <- list()
if (!is.null(dfs$TRA) && nrow(dfs$TRA)) out_list$TRA <- compute_epitope_overlap(dfs$TRA, "TRA", KEY_MODE,
                                                                                fixed_k = if (NORMALIZE_SCOPE=="global") global_min_k else NULL)
if (!is.null(dfs$TRB) && nrow(dfs$TRB)) out_list$TRB <- compute_epitope_overlap(dfs$TRB, "TRB", KEY_MODE,
                                                                                fixed_k = if (NORMALIZE_SCOPE=="global") global_min_k else NULL)
dfJ_all <- bind_rows(out_list, .id = NULL)

if (!nrow(dfJ_all)) {
  msg("No rows to plot. Try KEY_MODE='cdr3', set CDR3_REGEX <- NULL, or check TOP_EPI mapping.")
  quit(save="no")
}

## ---- Multiple-testing correction (global/per-chain; off-diagonals only) --
dfJ_all$q <- NA_real_
if (APPLY_CORRECTION) {
  if (ADJUST_SCOPE == "global") {
    idx <- with(dfJ_all, !is.na(p) & (if (EXCLUDE_DIAGONAL_IN_ADJ) E1 != E2 else TRUE))
    dfJ_all$q[idx] <- p.adjust(dfJ_all$p[idx], method = CORRECTION_METHOD)
  } else if (ADJUST_SCOPE == "per_chain") {
    dfJ_all <- dfJ_all %>%
      group_by(Chain) %>%
      mutate(q = {
        idx <- !is.na(p) & (if (EXCLUDE_DIAGONAL_IN_ADJ) E1 != E2 else TRUE)
        out <- rep(NA_real_, n())
        out[idx] <- p.adjust(p[idx], method = CORRECTION_METHOD)
        out
      }) %>% ungroup()
  } else stop("Unknown ADJUST_SCOPE: ", ADJUST_SCOPE)
}

# Choose p for stars
dfJ_all <- dfJ_all %>% mutate(p_used = if (APPLY_CORRECTION) q else p)
dfJ_all$stars <- p_to_stars(dfJ_all$p_used)
if (USE_FLOORS) dfJ_all$stars <- ifelse(dfJ_all$Inter >= MIN_INTER & dfJ_all$Jac >= MIN_JAC, dfJ_all$stars, "")

## ---- PLOT + SAVE ---------------------------------------------------------
p <- plot_epitope_heatmap_one_triangle(dfJ_all, FILL_MODE)
print(p)
ggsave(OUT_pdf_HEATMAP, p, width=12, height=8, dpi=300)
msg(paste0("Saved heatmap: ", normalizePath(OUT_pdf_HEATMAP)))

## ---- REPORT: detailed p-value calcs for TRA: InfluenzaA vs SARS-CoV-2 ----
report_pair <- function(df_chain, key_mode,
                        epiA, epiB,
                        nperm = max(10000, N_PERM),
                        fixed_k = if (isTRUE(NORMALIZE_TO_MIN_SET) && NORMALIZE_SCOPE=="global") global_min_k else NULL) {
  KE <- explode_keys_epitopes(df_chain, key_mode) %>% dplyr::filter(EPI %in% c(epiA, epiB))
  if (!nrow(KE)) stop("No keys for selected epitopes in the requested chain.")

  bm <- build_bin_matrix(KE, c(epiA, epiB),
                         normalize = NORMALIZE_TO_MIN_SET,
                         fixed_k   = fixed_k)
  bin <- bm$bin
  if (!length(bin)) stop("Empty bin after filtering/normalization.")

  a_vec <- bin[, 1]; b_vec <- bin[, 2]
  N <- nrow(bin); A <- sum(a_vec); B <- sum(b_vec); a <- sum(a_vec & b_vec)
  U <- sum(a_vec | b_vec); jac <- if (U == 0) NA_real_ else a / U
  Exp <- A * B / N

  p_hyp    <- if (A > 0 && B > 0) overlap_p_hypergeom(a, A, B, N) else NA_real_
  p_fisher <- if (A > 0 && B > 0) fisher_two_sided(a, A, B, N)    else NA_real_

  perm_vals <- replicate(nperm, { sum(sample(a_vec) & b_vec) })
  p_perm <- (sum(perm_vals >= a) + 1) / (nperm + 1)

  tibble::tibble(
    E1 = epiA, E2 = epiB,
    N = N, A = A, B = B, a = a, Union = U,
    Expected = Exp, Jaccard = jac,
    p_hypergeom = p_hyp, p_fisher2s = p_fisher, p_perm = p_perm,
    perm_median = median(perm_vals), perm_p95 = as.numeric(quantile(perm_vals, .95))
  )
}

if (!is.null(dfs[[REPORT_CHAIN]]) && nrow(dfs[[REPORT_CHAIN]])) {
  report_pair_one <- report_pair(dfs[[REPORT_CHAIN]], KEY_MODE, REPORT_EPI_A, REPORT_EPI_B)
  msg(paste0(REPORT_CHAIN, " ", REPORT_EPI_A, " vs ", REPORT_EPI_B, " overlap statistics:"))
  print(report_pair_one)

  pair_row <- dfJ_all %>%
    dplyr::filter(Chain == REPORT_CHAIN,
                  (E1 == REPORT_EPI_A & E2 == REPORT_EPI_B) |
                    (E1 == REPORT_EPI_B & E2 == REPORT_EPI_A)) %>%
    dplyr::select(Chain, E1, E2, Inter, A_size, B_size, N_universe, p, q, p_used, stars)

  msg(paste0("Heatmap-adjusted p (q) for ", REPORT_CHAIN, " ", REPORT_EPI_A, " vs ", REPORT_EPI_B, ":"))
  print(pair_row)
} else {
  msg("Requested REPORT_CHAIN has no data; skipping report.")
}

## ---- EXPORT: tidy p-value tables (post-normalization) --------------------
export_overlap_p_table <- function(df, file,
                                   include_diagonal = FALSE) {
  stopifnot(all(c("E1","E2","Inter","Union","Jac","A_size","B_size",
                  "N_universe","p","q","p_used","stars","Chain") %in% names(df)))

  out <- df %>%
    dplyr::mutate(diag = (E1 == E2)) %>%
    { if (!include_diagonal) dplyr::filter(., !diag) else . } %>%
    dplyr::select(Chain, E1, E2,
                  Inter, Union, Jac,
                  A_size, B_size, N_universe,
                  p_raw = p, q_BH = q, p_used, stars) %>%
    dplyr::arrange(Chain, E1, E2)

  readr::write_csv(out, file)
  msg(sprintf("Wrote %d rows → %s", nrow(out), normalizePath(file)))
  invisible(out)
}

p_table_all <- export_overlap_p_table(dfJ_all, OUT_CSV_ALL, include_diagonal = FALSE)

sig_df <- dfJ_all %>% dplyr::filter(!is.na(p_used) & p_used <= 0.05 & E1 != E2)
p_table_sig <- export_overlap_p_table(sig_df, OUT_CSV_SIG, include_diagonal = FALSE)

## ---- Optional overlap list files -----------------------------------------
if (!is.null(dfs$TRA) && nrow(dfs$TRA)) write_overlap_lists(dfs$TRA, "TRA", KEY_MODE, REPORT_EPI_A, REPORT_EPI_B, OUT_DIR)
if (!is.null(dfs$TRB) && nrow(dfs$TRB)) write_overlap_lists(dfs$TRB, "TRB", KEY_MODE, REPORT_EPI_A, REPORT_EPI_B, OUT_DIR)

msg(paste0("All outputs written to: ", normalizePath(OUT_DIR)))
