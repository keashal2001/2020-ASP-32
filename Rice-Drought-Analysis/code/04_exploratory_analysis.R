
# ==============================================================================
# Script: 04_exploratory_analysis.R
# Project: Rice Drought Stress Analysis - GSE115826
# PART 4: Exploratory Data Analysis (EDA)
# Analysis: PCA and Sample Correlation
# ==============================================================================

options(stringsAsFactors = FALSE)

cat("\n========== PART 4 STARTED ==========\n")

# 1. Check project directory
if (!dir.exists("data") && dir.exists("../data")) {
  setwd("..")
}

if (!dir.exists("data")) {
  stop("ERROR: Project root not found. Set working directory to project root.")
}

# 2. Load required package
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages(
    "ggplot2",
    repos = "https://cloud.r-project.org"
  )
}

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("ERROR: ggplot2 installation failed.")
}

library(ggplot2)

# 3. Create output folders
dir.create("results", showWarnings = FALSE, recursive = TRUE)
dir.create("figures", showWarnings = FALSE, recursive = TRUE)

# 4. Load cleaned data safely
clean_file <- "data/02_cleaned_data.RData"

if (!file.exists(clean_file)) {
  stop("ERROR: data/02_cleaned_data.RData not found.")
}

data_env <- new.env()

loaded_objects <- load(
  clean_file,
  envir = data_env
)

if (!all(c("expr_clean", "sample_meta") %in% loaded_objects)) {
  stop("ERROR: expr_clean or sample_meta missing from cleaned data.")
}

expr_clean <- as.matrix(data_env$expr_clean)
storage.mode(expr_clean) <- "numeric"

sample_meta <- as.data.frame(data_env$sample_meta)

# Validate expression matrix
if (is.null(rownames(expr_clean)) ||
    is.null(colnames(expr_clean))) {
  stop("ERROR: Expression matrix needs probe and sample names.")
}

if (any(!is.finite(expr_clean))) {
  stop("ERROR: Expression matrix contains missing or infinite values.")
}

if (anyDuplicated(rownames(expr_clean)) ||
    anyDuplicated(colnames(expr_clean))) {
  stop("ERROR: Duplicate probe or sample names detected.")
}

# Validate sample metadata
if (is.null(rownames(sample_meta))) {
  stop("ERROR: Sample metadata has no sample ID row names.")
}

if (!all(colnames(expr_clean) %in% rownames(sample_meta))) {
  stop("ERROR: Expression sample IDs do not match metadata.")
}

# Align metadata with expression matrix
sample_meta <- sample_meta[
  colnames(expr_clean),
  ,
  drop = FALSE
]

cat("\nDATASET INFORMATION\n")
cat("Number of probes:", nrow(expr_clean), "\n")
cat("Number of samples:", ncol(expr_clean), "\n")

# 5. Identify treatment groups using verified sample titles

if (!"title" %in% colnames(sample_meta)) {
  stop("ERROR: title column missing from sample metadata.")
}

sample_titles <- as.character(sample_meta$title)

treatment <- ifelse(
  grepl("^control", sample_titles, ignore.case = TRUE),
  "Control",
  ifelse(
    grepl("^drought", sample_titles, ignore.case = TRUE),
    "Drought",
    NA
  )
)

if (anyNA(treatment)) {
  print(data.frame(
    SampleID = rownames(sample_meta),
    Title = sample_titles,
    Treatment = treatment
  ))
  
  stop("ERROR: Could not identify all treatment groups.")
}

sample_meta$Treatment <- factor(
  treatment,
  levels = c("Control", "Drought")
)

sample_info <- data.frame(
  SampleID = colnames(expr_clean),
  Treatment = sample_meta$Treatment
)

cat("\nTREATMENT GROUPS\n")
print(table(sample_info$Treatment))

if (any(table(sample_info$Treatment) != 4)) {
  stop("ERROR: Expected four Control and four Drought samples.")
}

# 6. Principal Component Analysis (PCA)

cat("\nCalculating PCA...\n")

# Calculate variance for each probe
probe_variance <- apply(expr_clean, 1, var)

# Remove probes with zero or invalid variance
valid_probes <- is.finite(probe_variance) &
  probe_variance > 0

expr_pca <- expr_clean[
  valid_probes,
  ,
  drop = FALSE
]

if (nrow(expr_pca) < 2 || ncol(expr_pca) < 3) {
  stop("ERROR: Not enough variable probes or samples for PCA.")
}

cat("Variable probes used for PCA:", nrow(expr_pca), "\n")

# PCA on samples
# Do not scale individual probes because these are processed
# two-colour microarray log-ratio measurements.

pca_res <- prcomp(
  t(expr_pca),
  center = TRUE,
  scale. = FALSE
)

if (ncol(pca_res$x) < 2) {
  stop("ERROR: PCA produced fewer than two components.")
}

# Calculate percentage of variance explained
var_explained <- 100 *
  pca_res$sdev^2 / sum(pca_res$sdev^2)

pca_df <- data.frame(
  SampleID = rownames(pca_res$x),
  PC1 = pca_res$x[, 1],
  PC2 = pca_res$x[, 2],
  Treatment = sample_info$Treatment[
    match(rownames(pca_res$x), sample_info$SampleID)
  ]
)

if (anyNA(pca_df$Treatment)) {
  stop("ERROR: PCA samples could not be matched to treatments.")
}

cat("\nPCA RESULTS\n")
cat("PC1 variance:", round(var_explained[1], 2), "%\n")
cat("PC2 variance:", round(var_explained[2], 2), "%\n")

print(pca_df)

# Save PCA scores
write.csv(
  pca_df,
  "results/pca_sample_scores.csv",
  row.names = FALSE
)

# Save variance explained
variance_df <- data.frame(
  Principal_Component = paste0(
    "PC",
    seq_along(var_explained)
  ),
  Variance_Explained_Percent = round(
    var_explained,
    4
  )
)

write.csv(
  variance_df,
  "results/pca_variance_explained.csv",
  row.names = FALSE
)

# 7. Create PCA plot

treatment_colors <- c(
  "Control" = "#1e3a8a",
  "Drought" = "#c2410c"
)

p_pca <- ggplot(
  pca_df,
  aes(
    x = PC1,
    y = PC2,
    color = Treatment,
    shape = Treatment
  )
) +
  geom_point(size = 4) +
  geom_text(
    aes(label = SampleID),
    vjust = -1,
    size = 3,
    show.legend = FALSE
  ) +
  scale_color_manual(values = treatment_colors) +
  scale_shape_manual(
    values = c(
      "Control" = 16,
      "Drought" = 17
    )
  ) +
  theme_bw() +
  labs(
    title = "PCA of Rice Shoot Expression Profiles",
    subtitle = "GSE115826: Control vs Drought",
    x = paste0(
      "PC1 (",
      round(var_explained[1], 2),
      "%)"
    ),
    y = paste0(
      "PC2 (",
      round(var_explained[2], 2),
      "%)"
    ),
    color = "Treatment",
    shape = "Treatment"
  ) +
  theme(
    axis.text = element_text(color = "black"),
    plot.title = element_text(face = "bold")
  )

ggsave(
  "figures/fig3_pca_plot.png",
  plot = p_pca,
  width = 8,
  height = 6,
  dpi = 300
)

cat("PCA plot saved successfully.\n")

# 8. Sample Pearson correlation

cat("\nCalculating sample correlations...\n")

cor_matrix <- cor(
  expr_clean,
  method = "pearson"
)

if (any(!is.finite(cor_matrix))) {
  stop("ERROR: Invalid values in correlation matrix.")
}

# Use a separate filename to preserve Part 3 results
write.csv(
  cor_matrix,
  "results/eda_sample_correlation_matrix.csv"
)

# Convert matrix to data frame for plotting
cor_data <- as.data.frame(
  as.table(cor_matrix)
)

colnames(cor_data) <- c(
  "Sample1",
  "Sample2",
  "Correlation"
)

p_cor <- ggplot(
  cor_data,
  aes(
    x = Sample1,
    y = Sample2,
    fill = Correlation
  )
) +
  geom_tile(color = "white") +
  geom_text(
    aes(label = sprintf("%.2f", Correlation)),
    size = 3
  ) +
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0,
    limits = c(-1, 1)
  ) +
  theme_bw() +
  labs(
    title = "Sample Pearson Correlation",
    subtitle = "GSE115826 Rice Shoot Microarray",
    x = "Samples",
    y = "Samples",
    fill = "Correlation"
  ) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    ),
    axis.text.y = element_text(size = 8),
    plot.title = element_text(face = "bold")
  )

ggsave(
  "figures/fig4_sample_correlation_heatmap.png",
  plot = p_cor,
  width = 9,
  height = 8,
  dpi = 300
)

cat("Sample correlation heatmap saved successfully.\n")

# 9. Verify output files

output_files <- c(
  "results/pca_sample_scores.csv",
  "results/pca_variance_explained.csv",
  "results/eda_sample_correlation_matrix.csv",
  "figures/fig3_pca_plot.png",
  "figures/fig4_sample_correlation_heatmap.png"
)

verification <- data.frame(
  File = output_files,
  Exists = file.exists(output_files),
  Size_Bytes = ifelse(
    file.exists(output_files),
    file.info(output_files)$size,
    NA
  )
)

cat("\n========== OUTPUT VERIFICATION ==========\n")
print(verification)

if (!all(verification$Exists) ||
    any(verification$Size_Bytes <= 0, na.rm = TRUE)) {
  stop("ERROR: One or more output files are missing or empty.")
}

# 10. Final summary

cat("\n========== PART 4 COMPLETED SUCCESSFULLY ==========\n")

cat("Probes in cleaned matrix:", nrow(expr_clean), "\n")
cat("Variable probes used for PCA:", nrow(expr_pca), "\n")
cat("Samples analysed:", ncol(expr_clean), "\n")
cat("Control samples:", sum(treatment == "Control"), "\n")
cat("Drought samples:", sum(treatment == "Drought"), "\n")
cat("All five output files verified.\n")

