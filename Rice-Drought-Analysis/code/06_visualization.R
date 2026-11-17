
# ==============================================================================
# Script: 06_visualization.R
# Project: Rice Drought Stress Analysis - GSE115826
# PART 6: Scientific Visualization
# Figures 5, 6 and 7
# ==============================================================================

options(stringsAsFactors = FALSE)

cat("\n========== PART 6 STARTED ==========\n")

# 1. Check project directory
if (!dir.exists("data") && dir.exists("../data")) {
  setwd("..")
}

if (!dir.exists("data")) {
  stop("ERROR: Project root not found.")
}

# 2. Check required packages
required_packages <- c(
  "ggplot2",
  "pheatmap",
  "RColorBrewer",
  "dplyr",
  "tidyr",
  "tibble"
)

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(
      "ERROR: Missing package: ", pkg,
      ". Install it before running this script."
    )
  }
}

suppressPackageStartupMessages({
  library(ggplot2)
  library(pheatmap)
  library(RColorBrewer)
  library(dplyr)
  library(tidyr)
  library(tibble)
})

# 3. Create output folder
dir.create("figures", showWarnings = FALSE, recursive = TRUE)

# 4. Load Part 5 results safely
res_file <- "data/05_limma_results.RData"

if (!file.exists(res_file)) {
  stop("ERROR: data/05_limma_results.RData not found.")
}

data_env <- new.env()
loaded_objects <- load(res_file, envir = data_env)

required_objects <- c(
  "expr_clean",
  "sample_meta",
  "all_annotated",
  "sig_all",
  "sig_up",
  "sig_down"
)

missing_objects <- setdiff(required_objects, loaded_objects)

if (length(missing_objects) > 0) {
  stop(
    "ERROR: Missing objects: ",
    paste(missing_objects, collapse = ", ")
  )
}

expr_clean <- as.matrix(data_env$expr_clean)
sample_meta <- as.data.frame(data_env$sample_meta)
all_annotated <- as.data.frame(data_env$all_annotated)
sig_all <- as.data.frame(data_env$sig_all)
sig_up <- as.data.frame(data_env$sig_up)
sig_down <- as.data.frame(data_env$sig_down)

storage.mode(expr_clean) <- "numeric"

# 5. Validate data
if (any(!is.finite(expr_clean))) {
  stop("ERROR: Expression matrix contains non-finite values.")
}

if (is.null(rownames(expr_clean)) ||
    is.null(colnames(expr_clean))) {
  stop("ERROR: Expression matrix requires probe and sample names.")
}

if (!all(colnames(expr_clean) %in% rownames(sample_meta))) {
  stop("ERROR: Sample IDs do not match metadata.")
}

sample_meta <- sample_meta[
  colnames(expr_clean),
  ,
  drop = FALSE
]

# Identify treatment from verified sample titles
if (!"title" %in% colnames(sample_meta)) {
  stop("ERROR: Sample title column not found.")
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
  stop("ERROR: Could not identify all treatment groups.")
}

sample_info <- data.frame(
  SampleID = colnames(expr_clean),
  Treatment = factor(
    treatment,
    levels = c("Control", "Drought")
  )
)

cat("\nDATA VALIDATION\n")
cat("Probes:", nrow(expr_clean), "\n")
cat("Samples:", ncol(expr_clean), "\n")
print(table(sample_info$Treatment))

# Confirm Part 5 result columns
required_result_columns <- c(
  "ProbeID", "logFC", "adj.P.Val", "Description"
)

if (!all(required_result_columns %in% colnames(all_annotated))) {
  stop("ERROR: Required Part 5 result columns are missing.")
}

# Thresholds matching Part 5
pval_cutoff <- 0.05
lfc_cutoff <- log10(2)

# ==============================================================================
# FIGURE 5: VOLCANO PLOT
# ==============================================================================

cat("\nGenerating Figure 5: Volcano Plot...\n")

volcano_df <- all_annotated %>%
  mutate(
    neg_log10_padj = -log10(pmax(adj.P.Val, 1e-300)),
    Status = case_when(
      adj.P.Val < pval_cutoff & logFC >= lfc_cutoff ~ "Upregulated",
      adj.P.Val < pval_cutoff & logFC <= -lfc_cutoff ~ "Downregulated",
      TRUE ~ "Not Significant"
    )
  )

volcano_df$Status <- factor(
  volcano_df$Status,
  levels = c(
    "Not Significant",
    "Upregulated",
    "Downregulated"
  )
)

# Select the most statistically significant probes for labels
label_genes <- bind_rows(
  volcano_df %>%
    filter(Status == "Upregulated") %>%
    arrange(adj.P.Val) %>%
    head(6),
  volcano_df %>%
    filter(Status == "Downregulated") %>%
    arrange(adj.P.Val) %>%
    head(6)
)

label_genes$ShortName <- ifelse(
  is.na(label_genes$Description) |
    label_genes$Description == "",
  label_genes$ProbeID,
  substr(label_genes$Description, 1, 25)
)

volcano_colors <- c(
  "Not Significant" = "#94a3b8",
  "Upregulated" = "#dc2626",
  "Downregulated" = "#2563eb"
)

p_volcano <- ggplot(
  volcano_df,
  aes(
    x = logFC,
    y = neg_log10_padj,
    color = Status
  )
) +
  geom_point(alpha = 0.55, size = 1.2) +
  geom_vline(
    xintercept = c(-lfc_cutoff, lfc_cutoff),
    linetype = "dashed",
    color = "grey35"
  ) +
  geom_hline(
    yintercept = -log10(pval_cutoff),
    linetype = "dashed",
    color = "grey35"
  ) +
  geom_point(
    data = label_genes,
    color = "black",
    shape = 1,
    size = 2.5
  ) +
  geom_text(
    data = label_genes,
    aes(label = ShortName),
    color = "black",
    size = 2.8,
    vjust = -0.8,
    check_overlap = TRUE
  ) +
  scale_color_manual(values = volcano_colors) +
  theme_bw(base_size = 12) +
  labs(
    title = "Figure 5: Volcano Plot of Rice Drought-Responsive Probes",
    subtitle = "GSE115826 | FDR < 0.05 and at least 2-fold change",
    x = expression(log[10] ~ "Fold Change (Drought / Control)"),
    y = expression(-log[10] ~ "(Adjusted P-value)"),
    color = "Expression Status"
  ) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    legend.position = "top",
    panel.grid.minor = element_blank()
  )

ggsave(
  "figures/fig5_volcano_plot.png",
  plot = p_volcano,
  width = 9,
  height = 7,
  dpi = 300
)

cat("Figure 5 saved successfully.\n")

# ==============================================================================
# FIGURE 6: TOP 50 PROBE HEATMAP
# ==============================================================================

cat("\nGenerating Figure 6: Top 50 Probe Heatmap...\n")

# Select top 50 significant probes by adjusted p-value
top50_probes <- head(
  as.character(sig_all$ProbeID),
  50
)

# Validate probe IDs
top50_probes <- intersect(
  top50_probes,
  rownames(expr_clean)
)

if (length(top50_probes) < 2) {
  stop("ERROR: Not enough significant probes for heatmap.")
}

top50_mat <- expr_clean[
  top50_probes,
  ,
  drop = FALSE
]

# Remove probes with zero variance
row_variance <- apply(top50_mat, 1, var)

valid_heatmap_probes <- is.finite(row_variance) &
  row_variance > 0

top50_mat <- top50_mat[
  valid_heatmap_probes,
  ,
  drop = FALSE
]

if (nrow(top50_mat) < 2) {
  stop("ERROR: Not enough variable probes for heatmap.")
}

# Row-wise Z-score normalization
top50_z <- t(scale(t(top50_mat)))

if (any(!is.finite(top50_z))) {
  stop("ERROR: Non-finite values detected in heatmap matrix.")
}

# Sample annotation
annotation_col <- data.frame(
  Treatment = sample_info$Treatment
)

rownames(annotation_col) <- sample_info$SampleID

annotation_colors <- list(
  Treatment = c(
    "Control" = "#1e3a8a",
    "Drought" = "#c2410c"
  )
)

# Add probe descriptions
annotation_index <- match(
  rownames(top50_z),
  as.character(all_annotated$ProbeID)
)

probe_desc <- all_annotated$Description[annotation_index]

probe_desc[is.na(probe_desc) | probe_desc == ""] <- "Unknown"

rownames(top50_z) <- paste0(
  rownames(top50_z),
  ": ",
  substr(probe_desc, 1, 35)
)

png(
  "figures/fig6_deg_heatmap.png",
  width = 2700,
  height = 3000,
  res = 300
)

pheatmap(
  top50_z,
  main = paste(
    "Figure 6: Top Drought-Responsive Rice Probes",
    "Row Z-scores across 8 samples"
  ),
  color = colorRampPalette(
    rev(brewer.pal(11, "RdBu"))
  )(100),
  annotation_col = annotation_col,
  annotation_colors = annotation_colors,
  show_colnames = TRUE,
  fontsize_row = 7,
  fontsize_col = 9,
  clustering_distance_rows = "correlation",
  clustering_distance_cols = "euclidean",
  clustering_method = "complete",
  border_color = NA
)

dev.off()

cat("Figure 6 saved successfully.\n")
cat("Variable probes displayed:", nrow(top50_z), "\n")

# ==============================================================================
# FIGURE 7: SELECTED PROBE EXPRESSION PROFILES
# ==============================================================================

cat("\nGenerating Figure 7: Selected Probe Profiles...\n")

selected_probes <- c(
  "39577",
  "44040",
  "14199",
  "41300",
  "44259",
  "35235"
)

gene_names <- c(
  "39577" = "Isocitrate Lyase",
  "44040" = "Heat Shock Protein 70",
  "14199" = "16.9 kDa Heat Shock Protein",
  "41300" = "Male Sterility MS5",
  "44259" = "Laccase",
  "35235" = "OsNramp1"
)

# Confirm all selected probes exist
missing_probes <- setdiff(
  selected_probes,
  rownames(expr_clean)
)

if (length(missing_probes) > 0) {
  stop(
    "ERROR: Selected probes missing: ",
    paste(missing_probes, collapse = ", ")
  )
}

profile_mat <- expr_clean[
  selected_probes,
  ,
  drop = FALSE
]

profile_df <- as.data.frame(profile_mat) %>%
  rownames_to_column(var = "ProbeID") %>%
  pivot_longer(
    cols = -ProbeID,
    names_to = "SampleID",
    values_to = "Log10Ratio"
  ) %>%
  left_join(
    sample_info,
    by = "SampleID"
  ) %>%
  mutate(
    GeneLabel = factor(
      gene_names[ProbeID],
      levels = gene_names[selected_probes]
    )
  )

if (anyNA(profile_df$Treatment)) {
  stop("ERROR: Treatment information missing in profile data.")
}

# Calculate mean and standard error
profile_summary <- profile_df %>%
  group_by(GeneLabel, Treatment) %>%
  summarise(
    MeanExpr = mean(Log10Ratio),
    SEM = sd(Log10Ratio) / sqrt(n()),
    .groups = "drop"
  )

p_profiles <- ggplot(
  profile_summary,
  aes(
    x = Treatment,
    y = MeanExpr,
    fill = Treatment
  )
) +
  geom_col(
    width = 0.55,
    alpha = 0.85,
    color = "black"
  ) +
  geom_errorbar(
    aes(
      ymin = MeanExpr - SEM,
      ymax = MeanExpr + SEM
    ),
    width = 0.2
  ) +
  geom_jitter(
    data = profile_df,
    aes(
      x = Treatment,
      y = Log10Ratio,
      color = Treatment
    ),
    width = 0.12,
    size = 1.8,
    alpha = 0.75,
    inherit.aes = FALSE
  ) +
  facet_wrap(
    ~ GeneLabel,
    scales = "free_y",
    ncol = 3
  ) +
  scale_fill_manual(
    values = c(
      "Control" = "#1e3a8a",
      "Drought" = "#c2410c"
    )
  ) +
  scale_color_manual(
    values = c(
      "Control" = "#1e3a8a",
      "Drought" = "#c2410c"
    )
  ) +
  theme_bw(base_size = 11) +
  labs(
    title = "Figure 7: Expression Profiles of Selected Rice Probes",
    subtitle = "Mean and SEM across four biological replicates per condition",
    x = "Treatment",
    y = expression(log[10] ~ "(Cy5 / Cy3 ratio)")
  ) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    strip.text = element_text(face = "bold"),
    legend.position = "none",
    panel.grid.minor = element_blank()
  )

ggsave(
  "figures/fig7_top_genes_profiles.png",
  plot = p_profiles,
  width = 10,
  height = 7,
  dpi = 300
)

cat("Figure 7 saved successfully.\n")

# ==============================================================================
# FINAL OUTPUT VERIFICATION
# ==============================================================================

output_files <- c(
  "figures/fig5_volcano_plot.png",
  "figures/fig6_deg_heatmap.png",
  "figures/fig7_top_genes_profiles.png"
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
  stop("ERROR: One or more figure files are missing or empty.")
}

cat("\n========== PART 6 COMPLETED SUCCESSFULLY ==========\n")
cat("Three figures generated and verified.\n")

