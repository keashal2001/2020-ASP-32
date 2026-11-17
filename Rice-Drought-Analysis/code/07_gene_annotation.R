
# ==============================================================================
# Script: 07_gene_annotation.R
# Part 7: Gene Annotation and Functional Profiling
# Dataset: Rice drought stress - GSE115826
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
})

cat("\n========== PART 7 STARTED ==========\n")

# 1. Check project directory and output folders
if (!dir.exists("data") && dir.exists("../data")) {
  setwd("..")
}

if (!dir.exists("data")) {
  stop("ERROR: Project root not found.")
}

dir.create("results", showWarnings = FALSE, recursive = TRUE)

# 2. Load Part 5 results safely
res_file <- "data/05_limma_results.RData"

if (!file.exists(res_file)) {
  stop("ERROR: Part 5 results file not found.")
}

data_env <- new.env()
loaded_objects <- load(res_file, envir = data_env)

if (!"all_annotated" %in% loaded_objects) {
  stop("ERROR: all_annotated object not found.")
}

all_annotated <- as.data.frame(data_env$all_annotated)

# 3. Validate required columns
required_cols <- c(
  "ProbeID", "logFC", "P.Value", "adj.P.Val",
  "GB_ACC", "ProbeName", "GeneName",
  "Description", "Accessions"
)

missing_cols <- setdiff(required_cols, names(all_annotated))

if (length(missing_cols) > 0) {
  stop(
    "Missing annotation columns: ",
    paste(missing_cols, collapse = ", ")
  )
}

if (any(!is.finite(all_annotated$logFC)) ||
    any(!is.finite(all_annotated$adj.P.Val))) {
  stop("ERROR: Non-finite differential expression results detected.")
}

cat("\nDATA INFORMATION\n")
cat("Total probes:", nrow(all_annotated), "\n")

# 4. Extract RAP-DB and GenBank identifiers
cat("\nExtracting gene identifiers...\n")

accessions <- as.character(all_annotated$Accessions)

accession_parts <- strsplit(
  ifelse(is.na(accessions), "", accessions),
  "\\|",
  fixed = FALSE
)

rap_ids <- vapply(
  accession_parts,
  function(x) {
    if (length(x) >= 1 &&
        grepl("^Os[0-9]{2}g[0-9]+$", x[1])) {
      x[1]
    } else {
      NA_character_
    }
  },
  character(1)
)

gb_from_accessions <- vapply(
  accession_parts,
  function(x) {
    if (length(x) >= 2 &&
        nzchar(x[2]) &&
        x[2] != "NA") {
      x[2]
    } else {
      NA_character_
    }
  },
  character(1)
)

all_annotated$RAP_LocusID <- rap_ids

existing_gb <- as.character(all_annotated$GB_ACC)

use_existing_gb <- !is.na(existing_gb) &
  nzchar(existing_gb) &
  existing_gb != "NA"

all_annotated$GenBank_Acc <- ifelse(
  use_existing_gb,
  existing_gb,
  gb_from_accessions
)

# Do not estimate MSU IDs from RAP-DB coordinates.
# No unverified cross-reference mapping is introduced.

cat("Probes with RAP-DB IDs:",
    sum(!is.na(all_annotated$RAP_LocusID)), "\n")

cat("Probes with GenBank accessions:",
    sum(!is.na(all_annotated$GenBank_Acc)), "\n")

# 5. Assign preliminary functional categories
assign_pathway <- function(desc) {
  
  if (is.na(desc) || !nzchar(trimws(desc))) {
    return("Unclassified / Unknown")
  }
  
  desc <- tolower(desc)
  
  if (grepl("heat shock|chaperon|hsp|dna[jk]", desc)) {
    return("Chaperones and Heat Shock Response")
    
  } else if (grepl(
    "embryogenesis|lea|dehydrin|osmoprotect|aquaporin|proline|trehalose|mannitol|betaine",
    desc
  )) {
    return("Osmotic Adjustment and LEA Proteins")
    
  } else if (grepl(
    "transcription factor|zinc finger|bzip|nac|myb|wrky|ap2|erebp|homeobox|mads",
    desc
  )) {
    return("Transcription Factors and Gene Regulation")
    
  } else if (grepl(
    "peroxidase|superoxide|catalase|glutathione|ascorbate|redox|antioxidant",
    desc
  )) {
    return("ROS and Redox Scavenging")
    
  } else if (grepl(
    "laccase|pectin|cellulose|xyloglucan|cell wall|extensin|fasciclin",
    desc
  )) {
    return("Cell Wall Remodeling and Lignification")
    
  } else if (grepl(
    "chlorophyll|photosystem|rubisco|thylakoid|light-harvesting|plastid|calvin",
    desc
  )) {
    return("Photosynthesis and Chloroplast Function")
    
  } else if (grepl(
    "kinase|phosphatase|calcium|calmodulin|receptor",
    desc
  )) {
    return("Signal Transduction and Kinases")
    
  } else if (grepl(
    "transporter|channel|carrier|nramp|abc transporter|cation|anion|atpase",
    desc
  )) {
    return("Membrane Transport and Ion Homeostasis")
    
  } else if (grepl(
    "methylase|chromat|histone|dna methyl",
    desc
  )) {
    return("Epigenetics and Chromatin Remodeling")
    
  } else if (grepl(
    "isocitrate|glycolysis|tca|metabolism|synthase|reductase|lyase|dehydrogenase",
    desc
  )) {
    return("Primary and Secondary Metabolism")
    
  } else {
    return("Other Biological Functions")
  }
}

all_annotated$FunctionalCategory <- vapply(
  all_annotated$Description,
  assign_pathway,
  character(1)
)

# 6. Identify significant probes
# Match Part 5: FDR < 0.05 and at least 2-fold change.
fdr_cutoff <- 0.05
logfc_cutoff <- log10(2)

sig_annotated <- all_annotated %>%
  filter(
    adj.P.Val < fdr_cutoff,
    abs(logFC) >= logfc_cutoff
  )

up_annotated <- sig_annotated %>%
  filter(logFC >= logfc_cutoff) %>%
  arrange(adj.P.Val)

down_annotated <- sig_annotated %>%
  filter(logFC <= -logfc_cutoff) %>%
  arrange(adj.P.Val)

cat("\nDIFFERENTIAL EXPRESSION SUMMARY\n")
cat("Significant probes:", nrow(sig_annotated), "\n")
cat("Upregulated probes:", nrow(up_annotated), "\n")
cat("Downregulated probes:", nrow(down_annotated), "\n")

# 7. Count unique RAP-DB loci separately from probes
unique_rap <- function(x) {
  length(unique(na.omit(x$RAP_LocusID)))
}

cat("\nUNIQUE RAP-DB LOCUS COUNTS\n")
cat("All significant probes:",
    nrow(sig_annotated), "\n")
cat("Unique RAP-DB loci:",
    unique_rap(sig_annotated), "\n")
cat("Upregulated unique RAP-DB loci:",
    unique_rap(up_annotated), "\n")
cat("Downregulated unique RAP-DB loci:",
    unique_rap(down_annotated), "\n")

# 8. Export annotated probe tables
target_cols <- c(
  "ProbeID",
  "RAP_LocusID",
  "GenBank_Acc",
  "GB_ACC",
  "ORF",
  "ProbeName",
  "GeneName",
  "Description",
  "logFC",
  "FoldChange_Drought_vs_Control",
  "P.Value",
  "adj.P.Val",
  "FunctionalCategory"
)

target_cols <- intersect(target_cols, names(all_annotated))

write.csv(
  sig_annotated[, target_cols, drop = FALSE],
  "results/significant_probes_annotated.csv",
  row.names = FALSE,
  na = ""
)

write.csv(
  up_annotated[, target_cols, drop = FALSE],
  "results/upregulated_probes_annotated.csv",
  row.names = FALSE,
  na = ""
)

write.csv(
  down_annotated[, target_cols, drop = FALSE],
  "results/downregulated_probes_annotated.csv",
  row.names = FALSE,
  na = ""
)

# 9. Functional category summary
pathway_summary <- sig_annotated %>%
  group_by(FunctionalCategory) %>%
  summarise(
    Total_DEG_Probes = n(),
    Upregulated_Probes = sum(logFC >= logfc_cutoff),
    Downregulated_Probes = sum(logFC <= -logfc_cutoff),
    Unique_RAP_Loci = n_distinct(
      RAP_LocusID[!is.na(RAP_LocusID)]
    ),
    Top_Upregulated_Description =
      Description[which.max(logFC)][1],
    Maximum_log10FC = round(max(logFC), 4),
    Top_Downregulated_Description =
      Description[which.min(logFC)][1],
    Minimum_log10FC = round(min(logFC), 4),
    .groups = "drop"
  ) %>%
  arrange(desc(Total_DEG_Probes))

write.csv(
  pathway_summary,
  "results/functional_category_summary.csv",
  row.names = FALSE,
  na = ""
)

cat("\nFUNCTIONAL CATEGORY DISTRIBUTION\n")
print(pathway_summary)

# 10. Annotation coverage summary
annotation_summary <- data.frame(
  Metric = c(
    "Total biological probes",
    "Significant DEG probes",
    "Significant probes with RAP-DB IDs",
    "Significant probes with GenBank accessions",
    "Unique RAP-DB loci among significant probes",
    "Significant probes without RAP-DB IDs"
  ),
  Count = c(
    nrow(all_annotated),
    nrow(sig_annotated),
    sum(!is.na(sig_annotated$RAP_LocusID)),
    sum(!is.na(sig_annotated$GenBank_Acc)),
    unique_rap(sig_annotated),
    sum(is.na(sig_annotated$RAP_LocusID))
  )
)

write.csv(
  annotation_summary,
  "results/annotation_coverage_summary.csv",
  row.names = FALSE
)

cat("\nANNOTATION COVERAGE\n")
print(annotation_summary)

# 11. Save final annotated workspace
save(
  all_annotated,
  sig_annotated,
  up_annotated,
  down_annotated,
  pathway_summary,
  annotation_summary,
  file = "data/07_annotated_results.RData"
)

# 12. Verify output files
output_files <- c(
  "results/significant_probes_annotated.csv",
  "results/upregulated_probes_annotated.csv",
  "results/downregulated_probes_annotated.csv",
  "results/functional_category_summary.csv",
  "results/annotation_coverage_summary.csv",
  "data/07_annotated_results.RData"
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

cat("\nOUTPUT VERIFICATION\n")
print(verification)

if (!all(verification$Exists) ||
    any(verification$Size_Bytes <= 0, na.rm = TRUE)) {
  stop("ERROR: One or more output files are missing or empty.")
}

cat("\n========== PART 7 COMPLETED SUCCESSFULLY ==========\n")
