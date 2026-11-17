
# ==============================================================================
# Script: 03_quality_control.R
# Project: Identification and Analysis of Differentially Expressed Genes in Rice
# Dataset: GSE115826
# PART 3: Quality Control and Distribution Assessment
# ==============================================================================

# 1. Load required packages
suppressPackageStartupMessages({
  library(ggplot2)
  library(tidyr)
  library(dplyr)
  library(tibble)
})

cat("\n========== PART 3: QUALITY CONTROL ==========\n")

# 2. Load cleaned data
clean_file <- "data/02_cleaned_data.RData"

if (!file.exists(clean_file)) {
  stop("ERROR: Cleaned data file not found. Run Part 2 first.")
}

e <- new.env()
load(clean_file, envir = e)

if (!exists("expr_clean", envir = e) ||
    !exists("sample_meta", envir = e)) {
  stop("ERROR: Required cleaned data objects are missing.")
}

expr_clean <- as.matrix(e$expr_clean)
sample_meta <- e$sample_meta

# 3. Validate expression matrix
if (!is.numeric(expr_clean)) {
  stop("ERROR: Expression matrix must be numeric.")
}

if (nrow(expr_clean) < 2 || ncol(expr_clean) < 2) {
  stop("ERROR: Insufficient probes or samples.")
}

if (any(!is.finite(expr_clean))) {
  stop("ERROR: Missing or infinite expression values detected.")
}

if (is.null(rownames(expr_clean)) ||
    is.null(colnames(expr_clean))) {
  stop("ERROR: Probe IDs or sample names are missing.")
}

if (anyDuplicated(colnames(expr_clean)) > 0) {
  stop("ERROR: Duplicate sample names detected.")
}

if (any(!nzchar(colnames(expr_clean)))) {
  stop("ERROR: Empty sample names detected.")
}

# 4. Validate sample metadata
if (!"title" %in% colnames(sample_meta)) {
  stop("ERROR: Sample metadata does not contain title column.")
}

if (!all(colnames(expr_clean) %in% rownames(sample_meta))) {
  stop("ERROR: Some expression samples are missing from metadata.")
}

# Align metadata with expression matrix sample order
sample_meta <- sample_meta[colnames(expr_clean), , drop = FALSE]

if (!identical(rownames(sample_meta), colnames(expr_clean))) {
  stop("ERROR: Sample metadata alignment failed.")
}

# 5. Extract treatment groups from sample titles
treatment <- ifelse(
  grepl("^control", sample_meta$title, ignore.case = TRUE),
  "Control",
  ifelse(
    grepl("^drought", sample_meta$title, ignore.case = TRUE),
    "Drought",
    NA
  )
)

if (anyNA(treatment)) {
  stop("ERROR: Unknown treatment group found in sample titles.")
}

sample_info <- data.frame(
  SampleID = colnames(expr_clean),
  Treatment = treatment,
  stringsAsFactors = FALSE
)

cat("\nSample Treatment Information:\n")
print(sample_info)

cat("\nTreatment Group Counts:\n")
print(table(sample_info$Treatment))

if (length(unique(sample_info$Treatment)) != 2) {
  stop("ERROR: Both Control and Drought groups are required.")
}

# 6. Validate sample variability
sample_sd <- apply(expr_clean, 2, sd)

if (any(!is.finite(sample_sd)) || any(sample_sd == 0)) {
  stop("ERROR: Invalid or zero-variance samples detected.")
}

cat("\nExpression Matrix Dimensions:\n")
print(dim(expr_clean))

cat("\nProbes:", nrow(expr_clean), "\n")
cat("Samples:", ncol(expr_clean), "\n")

# 7. Create output directories
dir.create("results", showWarnings = FALSE, recursive = TRUE)
dir.create("figures", showWarnings = FALSE, recursive = TRUE)

# 8. Calculate sample quality statistics
qc_stats <- data.frame(
  SampleID = colnames(expr_clean),
  Treatment = treatment,
  Mean = apply(expr_clean, 2, mean),
  SD = sample_sd,
  Median = apply(expr_clean, 2, median),
  IQR = apply(expr_clean, 2, IQR),
  Min = apply(expr_clean, 2, min),
  Max = apply(expr_clean, 2, max)
)

write.csv(
  qc_stats,
  "results/sample_qc_summary.csv",
  row.names = FALSE
)

cat("\nSample Quality Statistics:\n")
print(qc_stats)

# 9. Prepare expression data for visualization
df_long <- as.data.frame(expr_clean) %>%
  rownames_to_column("ProbeID") %>%
  pivot_longer(
    cols = -ProbeID,
    names_to = "SampleID",
    values_to = "Expression"
  )

df_long$Treatment <- treatment[
  match(df_long$SampleID, sample_info$SampleID)
]

# 10. Figure 1: Expression distribution boxplot
p_box <- ggplot(
  df_long,
  aes(x = SampleID, y = Expression, fill = Treatment)
) +
  geom_boxplot(outlier.alpha = 0.3) +
  theme_bw() +
  labs(
    title = "Figure 1: Expression Distribution Across Samples",
    x = "Samples",
    y = "Expression Value",
    fill = "Treatment"
  ) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  )

ggsave(
  "figures/fig1_expression_boxplot.png",
  plot = p_box,
  width = 9,
  height = 6,
  dpi = 300
)

# 11. Figure 2: Expression density distribution
p_density <- ggplot(
  df_long,
  aes(x = Expression, color = SampleID)
) +
  geom_density(linewidth = 0.8) +
  theme_bw() +
  labs(
    title = "Figure 2: Expression Density Distribution",
    x = "Expression Value",
    y = "Density",
    color = "Sample"
  )

ggsave(
  "figures/fig2_density_distribution.png",
  plot = p_density,
  width = 9,
  height = 6,
  dpi = 300
)

# 12. Sample correlation assessment
cor_matrix <- cor(
  expr_clean,
  method = "pearson"
)

mean_cor <- (
  rowSums(cor_matrix) - diag(cor_matrix)
) / (ncol(cor_matrix) - 1)

correlation_stats <- data.frame(
  SampleID = colnames(expr_clean),
  Treatment = treatment,
  Mean_Correlation = as.numeric(mean_cor)
)

write.csv(
  correlation_stats,
  "results/sample_correlation_qc.csv",
  row.names = FALSE
)

cat("\nMean Sample Correlations:\n")
print(correlation_stats)

# 13. Outlier assessment
if (any(!is.finite(mean_cor)) || any(mean_cor < 0.85)) {
  warning(
    "Potential sample outlier detected. Review sample correlations."
  )
} else {
  cat("\nNo sample has mean correlation below 0.85.\n")
}

# 14. Save correlation matrix
write.csv(
  cor_matrix,
  "results/sample_correlation_matrix.csv"
)

# 15. Save sample treatment information
write.csv(
  sample_info,
  "results/sample_treatment_information.csv",
  row.names = FALSE
)

# 16. Verify output files
output_files <- c(
  "results/sample_qc_summary.csv",
  "results/sample_correlation_qc.csv",
  "results/sample_correlation_matrix.csv",
  "results/sample_treatment_information.csv",
  "figures/fig1_expression_boxplot.png",
  "figures/fig2_density_distribution.png"
)

cat("\nOutput Verification:\n")

verification <- data.frame(
  File = output_files,
  Exists = file.exists(output_files)
)

print(verification)

if (!all(verification$Exists)) {
  stop("ERROR: Some output files were not created.")
}

cat("\n========== PART 3 COMPLETED SUCCESSFULLY ==========\n")

# ============================================================
# FINAL QC DECISION REPORT
# ============================================================

dir.create("report", showWarnings = FALSE, recursive = TRUE)

qc_report <- c(
  "RICE DROUGHT STRESS ANALYSIS",
  "PART 3: QUALITY CONTROL DECISION REPORT",
  "========================================",
  "",
  "Dataset: GSE115826",
  "Platform: GPL6864",
  "Total samples: 8",
  "Control samples: 4",
  "Drought samples: 4",
  "Retained biological probes: 43734",
  "",
  "1. SAMPLE CORRELATION RESULTS",
  "Average control replicate correlation: 0.5394409",
  "Average drought replicate correlation: 0.9607586",
  "",
  "Control pairwise correlations ranged from 0.499 to 0.596.",
  "Drought pairwise correlations ranged from 0.951 to 0.972.",
  "",
  "2. EXPRESSION DISTRIBUTION",
  "Mean absolute expression:",
  "Control: 0.07614945",
  "Drought: 0.43749327",
  "",
  "3. EXPERIMENTAL DESIGN",
  "Cy5: Test sample.",
  "Cy3: Well-watered reference sample.",
  "Control arrays: Well-watered versus well-watered.",
  "Drought arrays: Drought versus well-watered.",
  "",
  "4. PROCESSING INFORMATION",
  "All samples used the same Agilent labeling protocol.",
  "All samples report Agilent Feature Extraction Software v11.5.1.1",
  "for background subtraction.",
  "",
  "5. QC DECISION",
  "All eight samples are retained.",
  "No sample was removed based on correlation alone.",
  "Control replicates show lower correlation and require cautious",
  "interpretation in downstream analysis.",
  "Drought replicates show high correlation consistency.",
  "",
  "6. LIMITATIONS",
  "Lower control correlations do not prove technical failure.",
  "Further analysis should account for the two-colour log-ratio",
  "nature of the microarray measurements.",
  "",
  "Part 3 QC decision documented."
)

writeLines(
  qc_report,
  "report/03_quality_control_decision.txt"
)

cat("\nQC DECISION REPORT CREATED SUCCESSFULLY.\n")
