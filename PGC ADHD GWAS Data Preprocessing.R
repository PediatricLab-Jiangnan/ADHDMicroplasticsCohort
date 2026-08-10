library(data.table)
library(TwoSampleMR)
library(dplyr)

# --- 1. Environment Setup ----------------------------------------------------
work_dir <- "C:/Users/xieru/Desktop/PEandADHD"  # Use forward slashes for cross-platform compatibility
setwd(work_dir)

# --- 2. Load Raw Data --------------------------------------------------------
raw_file <- "ADHD2022_iPSYCH_deCODE_PGC.meta.gz"

if (!file.exists(raw_file)) {
  stop("Raw data file not found: ", file.path(work_dir, raw_file))
}

cat("Loading data:", raw_file, "\n")
data <- fread(raw_file, header = TRUE, sep = "\t", data.table = FALSE)

cat("Raw data dimensions:", nrow(data), "rows ×", ncol(data), "columns\n")
cat("Key columns check:", 
    all(c("SNP", "OR", "SE", "A1", "A2", "P", "CHR", "BP", "EUR_FRQ") %in% names(data)), 
    "\n\n")

# --- 3. Core Transformation: OR → Beta ---------------------------------------
# Note: Odds Ratios must be natural-log-transformed to log(OR)
data$BETA <- log(data$OR)

# Check for invalid OR values (OR ≤ 0 would cause log() to fail)
if (any(data$OR <= 0, na.rm = TRUE)) {
  warning(sum(data$OR <= 0, na.rm = TRUE), " non-positive OR value(s) detected and set to NA")
  data$BETA[data$OR <= 0] <- NA
}

# --- 4. EAF (EUR_FRQ) Cleaning ------------------------------------------------
cat("\n--- EUR_FRQ Quality Check ---\n")

# Inspect original status
cat("Original EUR_FRQ type:", class(data$EUR_FRQ), "\n")
cat("Missing values:", sum(is.na(data$EUR_FRQ)), "\n")

# Force to character first, then numeric (handles factors/scientific notation issues)
data$EUR_FRQ_clean <- suppressWarnings(as.numeric(as.character(data$EUR_FRQ)))

# Check for conversion failures
failed_conversion <- sum(is.na(data$EUR_FRQ_clean)) - sum(is.na(data$EUR_FRQ))
if (failed_conversion > 0) {
  warning(failed_conversion, " EUR_FRQ value(s) could not be converted to numeric and were set to NA")
  
  # Optional: inspect specific problematic values (for debugging)
  bad_vals <- unique(data$EUR_FRQ[is.na(data$EUR_FRQ_clean) & !is.na(data$EUR_FRQ)])
  if (length(bad_vals) > 0) {
    cat("Problematic values (examples):", paste(head(bad_vals, 5), collapse = ", "), "\n")
  }
}

# EAF range validation (theoretically should be within [0, 1])
out_of_range <- sum(data$EUR_FRQ_clean < 0 | data$EUR_FRQ_clean > 1, na.rm = TRUE)
if (out_of_range > 0) {
  warning(out_of_range, " EUR_FRQ value(s) fall outside [0,1] range — please verify")
}

cat("Cleaned EUR_FRQ range:", 
    min(data$EUR_FRQ_clean, na.rm = TRUE), "-", 
    max(data$EUR_FRQ_clean, na.rm = TRUE), "\n")
cat("Final missing count:", sum(is.na(data$EUR_FRQ_clean)), "\n\n")

# --- 5. Format to TwoSampleMR Outcome Format ---------------------------------
cat("Converting to TwoSampleMR format...\n")

outcome_data <- format_data(
  dat         = data,
  type        = "outcome",
  snp_col     = "SNP",
  beta_col    = "BETA",           # Use the log(OR) transformed column
  se_col      = "SE",
  effect_allele_col  = "A1",
  other_allele_col   = "A2",
  pval_col    = "P",
  chr_col     = "CHR",
  pos_col     = "BP",
  eaf_col     = "EUR_FRQ_clean"   # Use cleaned EAF
)

# --- 6. Pre-export Quality Summary -------------------------------------------
cat("\n--- Output Data Quality Summary ---\n")
cat("Valid SNPs:", nrow(outcome_data), "\n")
cat("Beta range:", round(min(outcome_data$beta.outcome, na.rm = TRUE), 4), 
    "to", round(max(outcome_data$beta.outcome, na.rm = TRUE), 4), "\n")
cat("P-value range:", format(min(outcome_data$pval.outcome, na.rm = TRUE), digits = 2), 
    "to", format(max(outcome_data$pval.outcome, na.rm = TRUE), digits = 2), "\n")
cat("EAF missing rate:", round(sum(is.na(outcome_data$eaf.outcome)) / nrow(outcome_data) * 100, 2), "%\n")

# --- 7. Save Results ---------------------------------------------------------
output_file <- "ADHD_MRdata.txt"
write.table(
  outcome_data,
  file      = output_file,
  sep       = "\t",
  quote     = FALSE,
  row.names = FALSE
)

cat("\n✅ File saved:", file.path(work_dir, output_file), "\n")
