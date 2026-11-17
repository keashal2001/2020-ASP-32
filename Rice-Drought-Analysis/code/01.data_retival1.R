
# ==============================================================================
# Script: 01_data_retrieval.R
# Project: Identification and Analysis of Differentially Expressed Genes in Rice
#          (Oryza sativa) Under Drought Stress Using Bioinformatics Approaches
# Step: PART 1 - GEO Data Retrieval and Initial Preparation
# Dataset: GSE115826
# Platform: GPL6864
# ==============================================================================

# 1. INITIALIZATION
# ==============================================================================

message("==============================================================")
message("PART 1: GEO DATA RETRIEVAL")
message("==============================================================")

# Create project directories
if (!dir.exists("data")) {
  dir.create("data", recursive = TRUE)
}

if (!dir.exists("results")) {
  dir.create("results", recursive = TRUE)
}

# Required packages
required_packages <- c(
  "GEOquery",
  "Biobase"
)

# Install missing packages
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (pkg == "GEOquery") {
      if (!requireNamespace("BiocManager", quietly = TRUE)) {
        install.packages("BiocManager")
      }
      BiocManager::install("GEOquery", ask = FALSE, update = FALSE)
    } else {
      if (!requireNamespace("BiocManager", quietly = TRUE)) {
        install.packages("BiocManager")
      }
      BiocManager::install("Biobase", ask = FALSE, update = FALSE)
    }
  }
}

suppressPackageStartupMessages({
  library(GEOquery)
  library(Biobase)
})


# 2. DATASET INFORMATION
# ==============================================================================

gse_id <- "GSE115826"
gpl_id <- "GPL6864"

cat("\nDataset accession:", gse_id, "\n")
cat("Platform accession:", gpl_id, "\n")
cat("Organism: Oryza sativa (Rice)\n")


# 3. DOWNLOAD GEO SERIES MATRIX
# ==============================================================================

message("\nSTEP 1.1: Downloading GEO expression dataset...")

gse_list <- GEOquery::getGEO(
  GEO = gse_id,
  GSEMatrix = TRUE,
  getGPL = FALSE
)

if (length(gse_list) == 0) {
  stop("No GEO expression dataset was returned.")
}

# Select the correct platform
platform_names <- vapply(
  gse_list,
  function(x) Biobase::annotation(x),
  character(1)
)

cat("\nAvailable platforms:\n")
print(platform_names)

platform_index <- which(platform_names == gpl_id)

if (length(platform_index) != 1) {
  stop(
    "GPL6864 was not uniquely identified in the downloaded dataset. ",
    "Check the available platform names."
  )
}

eset <- gse_list[[platform_index]]

if (!inherits(eset, "ExpressionSet")) {
  stop("Downloaded object is not an ExpressionSet.")
}

cat("\nGEO ExpressionSet downloaded successfully.\n")

cat("\nExpressionSet summary:\n")
print(eset)


# 4. EXTRACT EXPRESSION MATRIX
# ==============================================================================

message("\nSTEP 1.2: Extracting expression matrix...")

expr_raw <- Biobase::exprs(eset)

expr_raw <- as.matrix(expr_raw)

storage.mode(expr_raw) <- "numeric"

if (is.null(rownames(expr_raw)) ||
    is.null(colnames(expr_raw))) {
  stop("Expression matrix is missing probe or sample identifiers.")
}

cat("\nExpression matrix dimensions:\n")
print(dim(expr_raw))

cat("\nFirst six probes:\n")
print(head(rownames(expr_raw)))

cat("\nFirst six samples:\n")
print(head(colnames(expr_raw)))

cat("\nExpression value summary:\n")
print(summary(as.vector(expr_raw)))


# 5. EXTRACT SAMPLE METADATA
# ==============================================================================

message("\nSTEP 1.3: Extracting sample metadata...")

sample_meta <- Biobase::pData(eset)

if (nrow(sample_meta) != ncol(expr_raw)) {
  stop("Sample metadata does not match expression matrix columns.")
}

# Ensure sample order is identical
sample_meta <- sample_meta[
  match(colnames(expr_raw), rownames(sample_meta)),
  ,
  drop = FALSE
]

stopifnot(
  identical(colnames(expr_raw), rownames(sample_meta))
)

cat("\nSample metadata dimensions:\n")
print(dim(sample_meta))

cat("\nSample metadata columns:\n")
print(colnames(sample_meta))

cat("\nSample titles:\n")
print(head(sample_meta$title))


# 6. DOWNLOAD GPL6864 PLATFORM ANNOTATION
# ==============================================================================

message("\nSTEP 1.4: Downloading GPL6864 platform annotation...")

# Download complete platform table
gpl <- GEOquery::getGEO(
  GEO = gpl_id,
  destdir = "data"
)

annot_table <- Biobase::Table(gpl)

if (!is.data.frame(annot_table) || nrow(annot_table) == 0) {
  stop("GPL6864 annotation table is empty.")
}

# Check required columns
if (!"ID" %in% colnames(annot_table)) {
  stop("GPL6864 annotation does not contain an ID column.")
}

if (!"ControlType" %in% colnames(annot_table)) {
  stop("GPL6864 annotation does not contain ControlType.")
}

annot_table$ID <- as.character(annot_table$ID)

cat("\nPlatform annotation dimensions:\n")
print(dim(annot_table))

cat("\nPlatform annotation columns:\n")
print(colnames(annot_table))

cat("\nFirst ten annotation records:\n")
print(head(annot_table, 10))


# 7. CHECK PLATFORM AND EXPRESSION PROBE MATCHING
# ==============================================================================

message("\nSTEP 1.5: Checking probe identifiers...")

if (anyDuplicated(rownames(expr_raw))) {
  warning("Duplicate expression probe identifiers detected.")
}

if (anyDuplicated(annot_table$ID)) {
  warning("Duplicate platform annotation IDs detected.")
}

matching_probes <- intersect(
  rownames(expr_raw),
  annot_table$ID
)

cat("\nProbe matching report:\n")

cat("Expression probes:",
    nrow(expr_raw), "\n")

cat("Platform probes:",
    nrow(annot_table), "\n")

cat("Common probes:",
    length(matching_probes), "\n")

cat("Unmatched expression probes:",
    sum(!rownames(expr_raw) %in% annot_table$ID), "\n")

if (length(matching_probes) == 0) {
  stop("No matching probes between expression matrix and GPL6864.")
}


# 8. INITIAL DATA QUALITY CHECKS
# ==============================================================================

message("\nSTEP 1.6: Initial data quality assessment...")

cat("\nMissing values:\n")
print(sum(is.na(expr_raw)))

cat("\nNaN values:\n")
print(sum(is.nan(expr_raw)))

cat("\nInfinite values:\n")
print(sum(is.infinite(expr_raw)))

cat("\nDuplicate expression probes:\n")
print(sum(duplicated(rownames(expr_raw))))


# 9. SAVE RAW INTERMEDIATE DATASET
# ==============================================================================

message("\nSTEP 1.7: Saving retrieved dataset...")

output_file <- file.path(
  "data",
  "01_raw_eset.RData"
)

save(
  eset,
  expr_raw,
  sample_meta,
  annot_table,
  platform_path = file.path("data", "GPL6864_platform.txt"),
  file = output_file
)

# Save platform annotation separately as a tab-delimited file
platform_path <- file.path(
  "data",
  "GPL6864_platform.txt"
)

write.table(
  annot_table,
  file = platform_path,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE,
  na = ""
)

if (!file.exists(output_file)) {
  stop("Failed to save GEO intermediate dataset.")
}

cat("\nSaved intermediate file:\n")
cat(output_file, "\n")

cat("\nSaved platform annotation:\n")
cat(platform_path, "\n")


# 10. FINAL VERIFICATION
# ==============================================================================

message("\nSTEP 1.8: Final verification...")

# Reload and verify saved data
check_env <- new.env()

load(
  output_file,
  envir = check_env
)

stopifnot(
  exists("eset", envir = check_env),
  exists("expr_raw", envir = check_env),
  exists("sample_meta", envir = check_env),
  exists("annot_table", envir = check_env),
  exists("platform_path", envir = check_env)
)

cat("\nVerified expression matrix:\n")
print(dim(check_env$expr_raw))

cat("\nVerified sample metadata:\n")
print(dim(check_env$sample_meta))

cat("\nVerified platform annotation:\n")
print(dim(check_env$annot_table))

cat("\nVerified platform file:\n")
print(check_env$platform_path)

message("\n==============================================================")
message("PART 1: DATA RETRIEVAL COMPLETED SUCCESSFULLY")
message("==============================================================")

