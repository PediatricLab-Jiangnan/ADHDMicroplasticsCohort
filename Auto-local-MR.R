# Load packages
library(VariantAnnotation)
library(gwasvcf)
library(gwasglue)
library(TwoSampleMR)
library(ieugwasr)
library(dplyr)
library(purrr)

# Set working directory (optional, does not affect subfolder reading)
setwd("C:/Users/xieru/Desktop/ADHD")

##############################
####### Process Exposure Data in Batches ######
##############################

# --- Core modification: Read all vcf.gz files from the 'exposure' subfolder ---
exposure_folder <- "./exposure" 

exposure_files <- list.files(
  path = exposure_folder, 
  pattern = "\\.vcf\\.gz$", 
  full.names = TRUE
)

# Check if files were found 
if(length(exposure_files) == 0) {
  stop("No .vcf.gz files found in the specified folder. Please check the path and file format.")
}
print(paste("Found", length(exposure_files), "exposure files:"))
print(exposure_files)

# Define processing function
process_exposure <- function(vcf_file) {
  
  # Read and process single exposure
  exposure_vcf <- readVcf(vcf_file)
  exposure_vcf_p_filter <- query_gwas(vcf = exposure_vcf, pval = 1e-05)
  exposure_data <- gwasvcf_to_TwoSampleMR(exposure_vcf_p_filter)
  
  # Clumping
  exposure_data_clumped <- ld_clump(
    dat = tibble(rsid = exposure_data$SNP,  
                 pval = exposure_data$pval.exposure, 
                 id = exposure_data$exposure),
    clump_kb = 10000,
    clump_r2 = 0.001,
    clump_p = 1,
    bfile = "./EUR",
    plink_bin = "./plink"
  )
  
  exposure_data %>% 
    filter(SNP %in% exposure_data_clumped$rsid) %>% 
    mutate(exposure_name = tools::file_path_sans_ext(basename(vcf_file)))
}

# Process all exposures in batch
exposure_list <- map(exposure_files, process_exposure)
names(exposure_list) <- tools::file_path_sans_ext(basename(exposure_files))

##############################
####### Process Outcome Data #########
##############################

# Read outcome data
outcome_vcf <- readVcf("ASD.vcf.gz")
outcome_data <- gwasvcf_to_TwoSampleMR(outcome_vcf)

##############################
####### MR Analysis Workflow ###########
##############################

# Define analysis function
run_mr_analysis <- function(exposure_data) {
  
  # Data integration
  data_common <- merge(exposure_data, outcome_data, by = "SNP")
  
  # Format outcome data
  outcome_formatted <- format_data(
    outcome_data,
    type = "outcome",
    snps = data_common$SNP,
    snp_col = "SNP",
    beta_col = "beta.exposure",
    se_col = "se.exposure",
    eaf_col = "eaf.exposure",
    effect_allele_col = "effect_allele.exposure",
    other_allele_col = "other_allele.exposure",
    pval_col = "pval.exposure",
    samplesize_col = "samplesize.exposure",
    ncase_col = "ncase.exposure",
    id_col = "id.exposure"
  )
  
  # Data harmonization
  harmonised_data <- harmonise_data(
    exposure_dat = exposure_data,
    outcome_dat = outcome_formatted
  )
  
  # Perform MR analysis
  list(
    mr_results = mr(harmonised_data),
    heterogeneity = mr_heterogeneity(harmonised_data),
    pleiotropy = mr_pleiotropy_test(harmonised_data),
    scatter_plot = mr_scatter_plot(mr(harmonised_data), harmonised_data)[[1]],
    sensitivity_plots = list(
      forest = mr_forest_plot(mr_singlesnp(harmonised_data))[[1]],
      leaveoneout = mr_leaveoneout_plot(mr_leaveoneout(harmonised_data))[[1]],
      funnel = mr_funnel_plot(mr_singlesnp(harmonised_data))[[1]]
    )
  )
}

##############################
####### Execute Batch Analysis #########
##############################

# Analyze all exposures
results <- map(exposure_list, safely(run_mr_analysis))

##############################
####### Output Results #############
##############################

# Create output directory
output_dir <- "MR_Results"
dir.create(output_dir, showWarnings = FALSE)

# Define save function
save_results <- function(result, exposure_name) {
  # Create exposure-specific directory
  exp_dir <- file.path(output_dir, exposure_name)
  dir.create(exp_dir, showWarnings = FALSE)
  
  # Save results
  write.csv(result$mr_results, file.path(exp_dir, "mr_results.csv"))
  write.csv(result$heterogeneity, file.path(exp_dir, "heterogeneity.csv"))
  write.csv(result$pleiotropy, file.path(exp_dir, "pleiotropy.csv"))
  
  # Save plots
  # Modified part
  ggplot2::ggsave(file.path(exp_dir, "scatter_plot.pdf"), result$scatter_plot, 
                  device = cairo_pdf, width=10, height=8, dpi=300)
  ggplot2::ggsave(file.path(exp_dir, "forest_plot.pdf"), result$sensitivity_plots$forest,
                  device = cairo_pdf, width=12, height=15, dpi=300)
  ggplot2::ggsave(file.path(exp_dir, "leaveoneout_plot.pdf"), result$sensitivity_plots$leaveoneout,
                  device = cairo_pdf, width=10, height=12, dpi=300)
  ggplot2::ggsave(file.path(exp_dir, "funnel_plot.pdf"), result$sensitivity_plots$funnel,
                  device = cairo_pdf, width=10, height=8, dpi=300)
}

# Save all results
iwalk(results, ～if(!is.null(.x$result)) save_results(.x$result, .y))
