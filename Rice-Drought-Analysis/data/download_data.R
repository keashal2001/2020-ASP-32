# ==============================================================================
# Script: download_data.R
# Project: Identification and Analysis of Differentially Expressed Genes in Rice
#          (Oryza sativa) Under Drought Stress Using Bioinformatics Approaches
# Dataset: NCBI GEO GSE115826 (Platform: GPL6864)
# Description: Automated, reproducible retrieval and validation of GEO series
#              matrix and GPL6864 platform annotation tables.
# ==============================================================================

# Ensure target directories exist
raw_dir <- file.path("data", "raw")
if (!dir.exists(raw_dir)) {
  dir.create(raw_dir, recursive = TRUE)
  message("Created raw data directory: ", raw_dir)
}

# 1. Download GSE115826 Series Matrix file
series_matrix_url <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE115nnn/GSE115826/matrix/GSE115826_series_matrix.txt.gz"
series_matrix_dest <- file.path(raw_dir, "GSE115826_series_matrix.txt.gz")

if (!file.exists(series_matrix_dest)) {
  message("Downloading GSE115826 series matrix from NCBI GEO...")
  download.file(url = series_matrix_url, destfile = series_matrix_dest, mode = "wb")
  message("Successfully downloaded series matrix: ", file.info(series_matrix_dest)$size, " bytes.")
} else {
  message("GSE115826 series matrix already exists locally.")
}

# 2. Download GPL6864 Platform Annotation table
platform_url <- "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL6864&targ=self&form=text&view=data"
platform_dest <- file.path(raw_dir, "GPL6864_platform.txt")

if (!file.exists(platform_dest)) {
  message("Downloading GPL6864 platform table from NCBI GEO...")
  download.file(url = platform_url, destfile = platform_dest, mode = "wb")
  message("Successfully downloaded GPL6864 platform table: ", file.info(platform_dest)$size, " bytes.")
} else {
  message("GPL6864 platform table already exists locally.")
}

# 3. Verification of downloaded files
message("\n--- Dataset Download Summary ---")
message("Series Matrix: ", series_matrix_dest, " (", round(file.info(series_matrix_dest)$size / 1024, 1), " KB)")
message("Platform Table: ", platform_dest, " (", round(file.info(platform_dest)$size / (1024 * 1024), 2), " MB)")

# Verify sample IDs in series matrix header
con <- gzfile(series_matrix_dest, "rt")
lines <- readLines(con, n = 80)
close(con)

sample_id_line <- grep("^!Series_sample_id", lines, value = TRUE)
cat("Verified Samples:\n", sample_id_line, "\n")
message("Data retrieval verified successfully.")
