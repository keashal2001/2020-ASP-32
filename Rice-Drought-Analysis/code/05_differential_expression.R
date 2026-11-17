
# ==============================================================================
# Script: 05_differential_expression.R
# Project: Rice Drought Stress Analysis - GSE115826
# PART 5: Differential Expression Analysis using limma
# ==============================================================================

options(stringsAsFactors = FALSE)

cat("\n========== PART 5 STARTED ==========\n")

# 1. Check project directory
if (!dir.exists("data") && dir.exists("../data")) {
  setwd("..")
}

if (!dir.exists("data")) {
  stop("ERROR: Project root not found. Set working directory to project root.")
}

# 2. Check limma package
if (!requireNamespace("limma", quietly = TRUE)) {
  stop("ERROR: limma package is missing. Install it using BiocManager.")
}

suppressPackageStartupMessages(library(limma))

# 3. Load cleaned data safely
clean_file <- "data/02_cleaned_data.RData"

if (!file.exists(clean_file)) {
  stop("ERROR: data/02_cleaned_data.RData not found.")
}

data_env <- new.env()
loaded_objects <- load(clean_file, envir = data_env)

required_objects <- c("expr_clean", "sample_meta", "annot_clean")

missing_objects <- setdiff(required_objects, loaded_objects)

if (length(missing_objects) > 0) {
  stop(
    "ERROR: Missing object(s) in cleaned data: ",
    paste(missing_objects, collapse = ", ")
  )
}

expr_clean <- as.matrix(data_env$expr_clean)
storage.mode(expr_clean) <- "numeric"

sample_meta <- as.data.frame(data_env$sample_meta)
annot_clean <- as.data.frame(data_env$annot_clean)

# 4. Validate expression data
if (any(!is.finite(expr_clean))) {
  stop("ERROR: Expression matrix contains missing or infinite values.")
}

if (is.null(rownames(expr_clean)) ||
    is.null(colnames(expr_clean))) {
  stop("ERROR: Expression matrix must have probe and sample names.")
}

if (anyDuplicated(rownames(expr_clean)) ||
    anyDuplicated(colnames(expr_clean))) {
  stop("ERROR: Duplicate probe or sample IDs detected.")
}

if (is.null(rownames(sample_meta))) {
  stop("ERROR: Sample metadata has no sample ID row names.")
}

if (!all(colnames(expr_clean) %in% rownames(sample_meta))) {
  stop("ERROR: Expression sample IDs do not match sample metadata.")
}

sample_meta <- sample_meta[
  colnames(expr_clean),
  ,
  drop = FALSE
]

cat("\nDATASET INFORMATION\n")
cat("Biological probes:", nrow(expr_clean), "\n")
cat("Samples:", ncol(expr_clean), "\n")

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

group <- factor(
  treatment,
  levels = c("Control", "Drought")
)

cat("\nSAMPLE GROUPS\n")
print(table(group))

if (any(table(group) != 4)) {
  stop("ERROR: Expected four Control and four Drought samples.")
}

# 6. Create design matrix
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)
rownames(design) <- colnames(expr_clean)

cat("\nDESIGN MATRIX\n")
print(design)

# 7. Fit linear models
cat("\nFitting limma linear models...\n")

fit <- lmFit(expr_clean, design)

# 8. Create Drought vs Control contrast
contrast_matrix <- makeContrasts(
  Drought_vs_Control = Drought - Control,
  levels = design
)

cat("\nCONTRAST MATRIX\n")
print(contrast_matrix)

fit_contrasts <- contrasts.fit(fit, contrast_matrix)

# 9. Empirical Bayes moderation
fit_ebayes <- eBayes(fit_contrasts)

# 10. Extract all probe-level results
all_results <- topTable(
  fit_ebayes,
  coef = "Drought_vs_Control",
  number = Inf,
  adjust.method = "BH",
  sort.by = "P"
)

all_results$ProbeID <- rownames(all_results)

# The input values are log10 ratios, so limma logFC is log10 fold change.
all_results$FoldChange_Drought_vs_Control <- 10^all_results$logFC

# 11. Validate platform annotation
required_annotation_columns <- c(
  "ID", "GB_ACC", "ORF", "ProbeName",
  "GeneName", "Description", "Accessions"
)

missing_annotation_columns <- setdiff(
  required_annotation_columns,
  colnames(annot_clean)
)

if (length(missing_annotation_columns) > 0) {
  stop(
    "ERROR: Missing annotation column(s): ",
    paste(missing_annotation_columns, collapse = ", ")
  )
}

annot_sub <- annot_clean[, required_annotation_columns]
annot_sub$ID <- as.character(annot_sub$ID)

if (anyDuplicated(annot_sub$ID)) {
  warning(
    "Duplicate annotation IDs found. Using the first annotation ",
    "entry for each probe to prevent duplicated result rows."
  )
}

# Match annotation without expanding the number of probe rows
annotation_index <- match(
  all_results$ProbeID,
  annot_sub$ID
)

annotation_data <- annot_sub[
  annotation_index,
  setdiff(required_annotation_columns, "ID"),
  drop = FALSE
]

all_annotated <- cbind(
  all_results,
  annotation_data
)

if (nrow(all_annotated) != nrow(expr_clean)) {
  stop("ERROR: Result row count does not match the number of probes.")
}

# Sort by adjusted p-value and then absolute fold change
all_annotated <- all_annotated[
  order(
    all_annotated$adj.P.Val,
    -abs(all_annotated$logFC)
  ),
]

# 12. Apply significance thresholds
pval_cutoff <- 0.05

# log10(2) corresponds to a 2-fold change
lfc_cutoff <- log10(2)

sig_up <- subset(
  all_annotated,
  adj.P.Val < pval_cutoff & logFC >= lfc_cutoff
)

sig_down <- subset(
  all_annotated,
  adj.P.Val < pval_cutoff & logFC <= -lfc_cutoff
)

sig_all <- subset(
  all_annotated,
  adj.P.Val < pval_cutoff & abs(logFC) >= lfc_cutoff
)

# 13. Print differential expression summary
cat("\n========== DIFFERENTIAL EXPRESSION SUMMARY ==========\n")

cat("Total biological probes tested:", nrow(all_annotated), "\n")
cat("FDR threshold:", pval_cutoff, "\n")
cat("Absolute fold-change threshold: 2-fold\n")
cat("Equivalent absolute log10FC threshold:", round(lfc_cutoff, 4), "\n")

cat("Significant probes:", nrow(sig_all), "\n")
cat("Upregulated probes:", nrow(sig_up), "\n")
cat("Downregulated probes:", nrow(sig_down), "\n")

cat("=====================================================\n")

# 14. Display top upregulated probes
cat("\nTOP 10 UPREGULATED PROBES\n")

if (nrow(sig_up) > 0) {
  print(head(
    sig_up[, c(
      "ProbeID", "logFC", "FoldChange_Drought_vs_Control",
      "AveExpr", "P.Value", "adj.P.Val", "Description"
    )],
    10
  ))
} else {
  cat("No probes meet the upregulation thresholds.\n")
}

# 15. Display top downregulated probes
cat("\nTOP 10 DOWNREGULATED PROBES\n")

if (nrow(sig_down) > 0) {
  print(head(
    sig_down[, c(
      "ProbeID", "logFC", "FoldChange_Drought_vs_Control",
      "AveExpr", "P.Value", "adj.P.Val", "Description"
    )],
    10
  ))
} else {
  cat("No probes meet the downregulation thresholds.\n")
}

# 16. Save results
dir.create("results", showWarnings = FALSE, recursive = TRUE)

write.csv(
  all_annotated,
  "results/all_genes_limma_results.csv",
  row.names = FALSE
)

write.csv(
  sig_all,
  "results/significant_degs_summary.csv",
  row.names = FALSE
)

write.csv(
  sig_up,
  "results/upregulated_probes.csv",
  row.names = FALSE
)

write.csv(
  sig_down,
  "results/downregulated_probes.csv",
  row.names = FALSE
)

# Save statistical model and results
save(
  expr_clean,
  sample_meta,
  fit_ebayes,
  all_annotated,
  sig_all,
  sig_up,
  sig_down,
  file = "data/05_limma_results.RData"
)

# 17. Verify output files
output_files <- c(
  "results/all_genes_limma_results.csv",
  "results/significant_degs_summary.csv",
  "results/upregulated_probes.csv",
  "results/downregulated_probes.csv",
  "data/05_limma_results.RData"
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

cat("\nOUTPUT FILE VERIFICATION\n")
print(verification)

if (!all(verification$Exists) ||
    any(verification$Size_Bytes <= 0, na.rm = TRUE)) {
  stop("ERROR: One or more output files are missing or empty.")
}

cat("\n========== PART 5 COMPLETED SUCCESSFULLY ==========\n")
cat("All five output files verified.\n")
