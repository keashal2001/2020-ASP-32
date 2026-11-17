# ==============================================================================
# Script: 02_data_cleaning.R
# Project: Identification and Analysis of Differentially Expressed Genes in Rice
#          (Oryza sativa) Under Drought Stress Using Bioinformatics Approaches
# Part 2: Data Preparation and Cleaning
# ==============================================================================

# 1. Load required package
suppressPackageStartupMessages({
  library(Biobase)
})

# 2. Set project directories
data_dir <- "data"
raw_file <- file.path(data_dir, "01_raw_eset.RData")
output_file <- file.path(data_dir, "02_cleaned_data.RData")

# 3. Check input file
if (!file.exists(raw_file)) {
  stop("Missing 01_raw_eset.RData. Check your working directory.")
}

# 4. Load raw ExpressionSet
e <- new.env()
load(raw_file, envir = e)

if (!exists("eset", envir = e)) {
  stop("The ExpressionSet object 'eset' was not found.")
}

eset <- e$eset

# 5. Extract expression matrix and sample metadata
expr_raw <- Biobase::exprs(eset)
sample_meta <- Biobase::pData(eset)

# 6. Load platform annotation
platform_file <- file.path("data", "raw", "GPL6864_platform.txt")

if (!file.exists(platform_file)) {
  stop("Missing GPL6864_platform.txt.")
}

platform <- read.delim(
  platform_file,
  skip = 14,
  header = TRUE,
  sep = "\t",
  quote = "",
  check.names = FALSE,
  stringsAsFactors = FALSE,
  fill = TRUE
)

# 7. Validate probe IDs and annotation columns
if (is.null(rownames(expr_raw))) {
  stop("Expression matrix has no probe IDs.")
}

if (anyDuplicated(rownames(expr_raw)) > 0) {
  stop("Duplicate probe IDs detected.")
}

required_columns <- c("ID", "ControlType")
missing_columns <- setdiff(required_columns, colnames(platform))

if (length(missing_columns) > 0) {
  stop(
    "Platform annotation is missing required column(s): ",
    paste(missing_columns, collapse = ", ")
  )
}

# 8. Match platform annotation to expression probes
annot_index <- match(rownames(expr_raw), as.character(platform$ID))

if (anyNA(annot_index)) {
  stop("Some expression probe IDs do not match platform annotation.")
}

annot_sub <- platform[annot_index, , drop = FALSE]

# 9. Keep biological probes only (ControlType == 0)
control_type <- suppressWarnings(
  as.numeric(as.character(annot_sub$ControlType))
)

biological_mask <- !is.na(control_type) & control_type == 0

expr_clean <- expr_raw[biological_mask, , drop = FALSE]
annot_clean <- annot_sub[biological_mask, , drop = FALSE]

# 10. Validate cleaned data and sample alignment
if (nrow(expr_clean) == 0) {
  stop("No biological probes remain after filtering.")
}

if (anyNA(expr_clean) || any(!is.finite(expr_clean))) {
  stop("Missing or infinite values detected in cleaned expression data.")
}

if (!identical(rownames(expr_clean), as.character(annot_clean$ID))) {
  stop("Expression and annotation probe IDs are not aligned.")
}

if (!identical(colnames(expr_clean), rownames(sample_meta))) {
  stop("Expression columns and sample metadata rows are not aligned.")
}

if (!"title" %in% colnames(sample_meta)) {
  stop("Sample metadata does not contain a 'title' column.")
}

# 11. Report cleaning results
cat("\n===== PART 2: DATA CLEANING SUMMARY =====\n")
cat("Original probes:", nrow(expr_raw), "\n")
cat("Biological probes retained:", nrow(expr_clean), "\n")
cat("Removed probes:", nrow(expr_raw) - nrow(expr_clean), "\n")
cat("Samples:", ncol(expr_clean), "\n")
cat("Missing values:", sum(is.na(expr_clean)), "\n")
cat("Infinite values:", sum(!is.finite(expr_clean)), "\n")
cat(
  "Control samples:",
  sum(grepl("^control", sample_meta$title, ignore.case = TRUE)),
  "\n"
)
cat(
  "Drought samples:",
  sum(grepl("^drought", sample_meta$title, ignore.case = TRUE)),
  "\n"
)

# 12. Save cleaned data
save(
  expr_clean,
  annot_clean,
  sample_meta,
  file = output_file
)

cat("\nCleaned data saved to:", output_file, "\n")
cat("PART 2 SCRIPT COMPLETED SUCCESSFULLY.\n")
