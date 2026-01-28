###############################################################################
# END-TO-END: Nested models (cumulative features) line plot
# - Fits nested logistic regression models on TRA only
# - Computes AUC on TRAIN and TEST for each nested model
# - Produces Nature-style black/white line plot (solid=train, dashed=test)
# - Saves to your OUT_DIR
###############################################################################

rm(list = ls())

suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(readr)
  library(ggplot2)
  library(pROC)
  library(Peptides)
})

###############################################################################
# 0) GLOBAL FONT / SIZE CONTROLS (edit here only)
###############################################################################
FS_BASE       <- 22
FS_AXIS_TEXT  <- 18
FS_AXIS_TITLE <- 22
FS_LEGEND     <- 18
FS_POINT      <- 4.2
FS_LINE       <- 1.2
PLOT_W_IN     <- 6
PLOT_H_IN     <- 6

###############################################################################
# 1) USER SETTINGS
###############################################################################
DATA_FILE <- "~/Desktop/Shubhams_paper_code/data/hu_with_metaclone_20250826_095644.csv"
PGEN_FILE <- "~/Desktop/Shubhams_paper_code/data/igor_output_combined.csv"

# If Pgen column is not auto-detected, set it explicitly (e.g. "Pgen" or "pgen_nt")
PGEN_COL_MANUAL <- NULL

# Target definition:
# "alpha_vs_noncentric" : alphaCentric vs (none + betaCentric)
# "alpha_vs_none"       : alphaCentric vs none (drops betaCentric)
TARGET_MODE <- "alpha_vs_noncentric"

# Train/test split
SET_SEED   <- 42
TRAIN_FRAC <- 0.80

# OUTPUTS
OUT_DIR <- file.path(getwd(), "Figure_Images")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

###############################################################################
# 2) HELPERS
###############################################################################
stop_if_missing <- function(df, cols, df_name = "data.frame") {
  miss <- setdiff(cols, names(df))
  if (length(miss) > 0) {
    stop(sprintf("Missing required columns in %s: %s", df_name, paste(miss, collapse = ", ")))
  }
}

strip_allele <- function(x) toupper(sub("\\*.*$", "", trimws(x)))

detect_pgen_col <- function(df){
  preferred <- c("Pgen","pgen","pgen_nt","Pgen_nt","pgen_est","Pgen_est",
                 "log10_Pgen","log10pgen","lnPgen","logPgen","Pgen_value")
  hit <- preferred[preferred %in% names(df)]
  if (length(hit) > 0) return(hit[1])
  cand <- names(df)[sapply(df, is.numeric) & grepl("pgen", names(df), ignore.case = TRUE)]
  if (length(cand) > 0) return(cand[1])
  NA_character_
}

# Flexibility (Vihinen VINM940101)
aa_flex <- c(
  A=0.984, R=1.008, N=1.048, D=1.068, C=0.906,
  Q=1.037, E=1.094, G=1.031, H=0.950, I=0.927,
  L=0.935, K=1.102, M=0.952, F=0.915, P=1.049,
  S=1.046, T=0.997, W=0.904, Y=0.929, V=0.931
)
vihinen_flex <- function(seq){
  aa <- strsplit(seq, "")[[1]]
  vals <- aa_flex[aa]
  if (any(is.na(vals))) return(NA_real_)
  mean(vals)
}

# Train-only scaling then apply to test
scale_train_apply <- function(train_df, test_df, num_cols){
  mu <- sapply(train_df[num_cols], mean, na.rm = TRUE)
  sdv <- sapply(train_df[num_cols], sd, na.rm = TRUE)
  sdv[sdv == 0 | is.na(sdv)] <- 1

  train_s <- train_df
  test_s  <- test_df
  for (cc in num_cols){
    train_s[[cc]] <- (train_s[[cc]] - mu[[cc]]) / sdv[[cc]]
    test_s[[cc]]  <- (test_s[[cc]]  - mu[[cc]]) / sdv[[cc]]
  }
  list(train_s=train_s, test_s=test_s, mu=mu, sdv=sdv)
}

auc_on_df <- function(fit, df){
  pr <- predict(fit, newdata = df, type = "response")
  as.numeric(pROC::auc(pROC::roc(df$y, pr, quiet = TRUE)))
}

###############################################################################
# 3) LOAD + STANDARDIZE HU/VDJdb-DERIVED FILE
###############################################################################
hu <- read.csv(DATA_FILE, check.names = FALSE, stringsAsFactors = FALSE)

stop_if_missing(
  hu,
  cols = c("complex.id", "cdr3", "v.segm", "j.segm", "metaclone"),
  df_name = "HU/VDJdb CSV"
)

hu <- hu %>%
  mutate(
    libid      = as.character(complex.id),
    junction   = as.character(cdr3),
    v_gene_raw = as.character(v.segm),
    j_gene_raw = as.character(j.segm),
    Type       = as.character(metaclone)
  ) %>%
  mutate(
    Type = case_when(
      str_detect(tolower(Type), "alpha") ~ "alphaCentric",
      str_detect(tolower(Type), "beta")  ~ "betaCentric",
      TRUE                               ~ "none"
    )
  ) %>%
  filter(str_detect(junction, "^C"),
         str_detect(junction, "[FW]$")) %>%
  mutate(
    v_gene = strip_allele(v_gene_raw),
    j_gene = strip_allele(j_gene_raw),
    chain = case_when(
      grepl("^TRAV", v_gene) | grepl("^TRAJ", j_gene) ~ "TRA",
      grepl("^TRBV", v_gene) | grepl("^TRBJ", j_gene) ~ "TRB",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(chain %in% c("TRA","TRB"))

# Remove MAIT + iNKT (TRA-only heuristic)
mait_inkt_libs <- hu %>%
  filter(
    chain == "TRA" &
      (
        junction == "CVVSDRGSTLGRLYF" |
          (v_gene == "TRAV1-2" & j_gene %in% c("TRAJ33","TRAJ20","TRAJ12"))
      )
  ) %>% pull(libid)

hu <- hu %>% filter(!(libid %in% mait_inkt_libs))

# Require paired TRA+TRB
paired_ids <- hu %>%
  distinct(libid, chain) %>%
  count(libid, name = "n_chain") %>%
  filter(n_chain >= 2) %>%
  pull(libid)

hu <- hu %>% filter(libid %in% paired_ids)

# Deduplicate
df <- hu %>% distinct(libid, chain, junction, v_gene, j_gene, Type, .keep_all = TRUE)

###############################################################################
# 4) FEATURE ENGINEERING (AA features)
###############################################################################
df <- df %>%
  mutate(
    len  = nchar(junction),
    hyd  = Peptides::hydrophobicity(junction, scale = "Eisenberg"),
    flex = vapply(junction, vihinen_flex, numeric(1))
  ) %>%
  filter(!is.na(hyd), !is.na(flex))

###############################################################################
# 5) KEEP TRA + DEFINE y
###############################################################################
tra <- df %>%
  filter(chain == "TRA") %>%
  mutate(
    y = case_when(
      TARGET_MODE == "alpha_vs_noncentric" ~ ifelse(Type == "alphaCentric", 1L, 0L),
      TARGET_MODE == "alpha_vs_none" ~ ifelse(Type == "alphaCentric", 1L,
                                             ifelse(Type == "none", 0L, NA_integer_)),
      TRUE ~ NA_integer_
    ),
    j_gene = factor(j_gene)
  ) %>%
  filter(!is.na(y))

###############################################################################
# 6) LOAD PGEN + MERGE BY libid
###############################################################################
values <- read.csv(PGEN_FILE, stringsAsFactors = FALSE)
stop_if_missing(values, cols = c("TCR_name"), df_name = "Pgen CSV")

values$libid <- as.character(gsub("^TCR_", "", values$TCR_name))

pgen_col <- if (!is.null(PGEN_COL_MANUAL)) PGEN_COL_MANUAL else detect_pgen_col(values)
if (is.na(pgen_col) || !(pgen_col %in% names(values))) {
  stop("Could not detect a Pgen column in PGEN_FILE. Set PGEN_COL_MANUAL to the correct column name.")
}

pgen_tbl <- values %>%
  transmute(
    libid = as.character(libid),
    Pgen_raw = suppressWarnings(as.numeric(.data[[pgen_col]]))
  ) %>%
  distinct(libid, .keep_all = TRUE)

tra <- tra %>%
  left_join(pgen_tbl, by = "libid")

eps <- 1e-300
tra <- tra %>%
  mutate(log10Pgen = ifelse(is.na(Pgen_raw), NA_real_, log10(Pgen_raw + eps)))

tra_model <- tra %>%
  filter(!is.na(log10Pgen)) %>%
  select(libid, y, len, hyd, flex, log10Pgen, j_gene)

cat("\nRows in tra_model:", nrow(tra_model), "\n")
cat("Positive fraction:", mean(tra_model$y == 1), "\n")

###############################################################################
# 7) TRAIN/TEST SPLIT (STRATIFIED)
###############################################################################
set.seed(SET_SEED)
idx_pos <- which(tra_model$y == 1)
idx_neg <- which(tra_model$y == 0)

train_pos <- sample(idx_pos, size = ceiling(length(idx_pos) * TRAIN_FRAC))
train_neg <- sample(idx_neg, size = ceiling(length(idx_neg) * TRAIN_FRAC))
train_idx <- sort(c(train_pos, train_neg))
test_idx  <- setdiff(seq_len(nrow(tra_model)), train_idx)

train <- tra_model[train_idx, , drop = FALSE]
test  <- tra_model[test_idx,  , drop = FALSE]

# Scale numeric predictors using TRAIN stats only
num_cols <- c("len","hyd","flex","log10Pgen")
scaled <- scale_train_apply(train, test, num_cols)
train_s <- scaled$train_s
test_s  <- scaled$test_s

###############################################################################
# 8) FIT NESTED MODELS + COMPUTE AUC ON TRAIN & TEST
###############################################################################
formulas <- list(
  null      = y ~ 1,
  plus_len   = y ~ len,
  plus_flex  = y ~ len + flex,
  plus_hyd   = y ~ len + flex + hyd,
  plus_pgen  = y ~ len + flex + hyd + log10Pgen
  # If you want to include TRAJ usage as a final step, uncomment:
  # , plus_traj = y ~ len + flex + hyd + log10Pgen + j_gene
)

model_order <- names(formulas)

auc_rows <- lapply(model_order, function(nm){
  fit <- glm(formulas[[nm]], data = train_s, family = binomial())
  tibble(
    model = nm,
    auc_train = auc_on_df(fit, train_s),
    auc_test  = auc_on_df(fit, test_s)
  )
}) %>% bind_rows() %>%
  mutate(
    step = match(model, model_order),
    model = factor(model, levels = model_order)
  )

print(auc_rows)

# Long format for ggplot
auc_long <- auc_rows %>%
  select(step, model, auc_train, auc_test) %>%
  tidyr::pivot_longer(cols = c(auc_train, auc_test),
                      names_to = "dataset", values_to = "auc") %>%
  mutate(
    dataset = recode(dataset, auc_train = "training", auc_test = "test"),
    dataset = factor(dataset, levels = c("training","test"))
  )

###############################################################################
# 9) PLOT (match your attached style)
###############################################################################
p <- ggplot(auc_long, aes(x = model, y = auc, group = dataset,
                          linetype = dataset, shape = dataset)) +
  geom_line(color = "black", linewidth = FS_LINE) +
  geom_point(color = "black", size = FS_POINT) +
  scale_linetype_manual(values = c(training = "solid", test = "dashed")) +
  scale_shape_manual(values = c(training = 16, test = 17)) +
  labs(
    x = "\nModel features, cumulative",
    y = "AUC\n",
    linetype = "Dataset",
    shape = "Dataset"
  ) +
  theme_classic(base_size = FS_BASE) +
  theme(
    axis.line   = element_line(linewidth = 1.1),
    axis.text.x = element_text(angle = 35, hjust = 1, size = FS_AXIS_TEXT),
    axis.text.y = element_text(size = FS_AXIS_TEXT),
    axis.title  = element_text(size = FS_AXIS_TITLE, face = "bold"),
    legend.title = element_text(size = FS_LEGEND),
    legend.text  = element_text(size = FS_LEGEND),
    legend.position = c(0.82, 0.55)
  )

print(p)

###############################################################################
# 10) SAVE (PDF + PNG) TO OUT_DIR
###############################################################################
out_stub <- sprintf("Fig1E_nested_models_lineplot_%s_seed%d", TARGET_MODE, SET_SEED)
out_pdf  <- file.path(OUT_DIR, paste0(out_stub, ".pdf"))
out_png  <- file.path(OUT_DIR, paste0(out_stub, ".png"))

cat("\nOUT_DIR:", normalizePath(OUT_DIR, mustWork = FALSE), "\n")
cat("Saving PDF:", out_pdf, "\n")
cat("Saving PNG:", out_png, "\n")

# Prefer cairo PDF for consistent text rendering
pdf_device <- if (requireNamespace("Cairo", quietly = TRUE)) Cairo::CairoPDF else cairo_pdf

tryCatch({
  ggsave(out_pdf, plot = p, width = PLOT_W_IN, height = PLOT_H_IN, units = "in",
         device = pdf_device, useDingbats = FALSE)
  ggsave(out_png, plot = p, width = PLOT_W_IN, height = PLOT_H_IN, units = "in",
         dpi = 300)
  cat("Saved OK.\n")
}, error = function(e){
  message("SAVE FAILED: ", conditionMessage(e))
})

cat("\nFiles now in OUT_DIR matching 'Fig1E_nested_models':\n")
print(list.files(OUT_DIR, pattern = "Fig1E_nested_models", full.names = TRUE))

###############################################################################
# END
###############################################################################
