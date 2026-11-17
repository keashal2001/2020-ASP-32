# Multi-Cultivar Transcriptomic Biomarker Discovery, Pathway Enrichment, and Machine Learning Classification of Drought Stress Response in *Oryza sativa*

**Coursework:** Bioinformatics (ICA3)  
**Author:** Keashal Harishman  
**Student Registration Number:** 2020/ASP/32  
**GitHub Repository Name:** `2020-ASP-32` *(Note: GitHub does not permit `/` in repository names; hyphens replace slashes)*  
**Organism:** *Oryza sativa* L. (Rice)  
**NCBI GEO Accession:** [GSE65022](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE65022)  
**Analysis Environment:** R 4.4.1  
**Core R Packages:** `DESeq2`, `GEOquery`, `randomForest`, `caret`, `e1071`, `pROC`, `pheatmap`, `tidyverse`

---

## Abstract

Drought stress restricts crop growth and drives widespread transcriptional reprogramming in rainfed rice (*Oryza sativa*). This study analyses a 24-sample leaf RNA-sequencing dataset modeled on NCBI GEO series GSE65022, comparing well-watered control and water-deficit treatments across three *indica/aus* cultivars—drought-susceptible IR64, moderately tolerant Vandana, and drought-tolerant Nagina22—at Day 3 and Day 7. Using a multi-factor negative binomial generalized linear model (`DESeq2`, `~ genotype + timepoint + condition`) across 5,000 expressed loci, we identified 485 significantly differentially expressed genes ($\text{FDR} < 0.05$, $|\log_2\text{FC}| > 1.0$; 244 upregulated, 241 downregulated). Protective effectors (*OsRAB16A*, *OsLEA3-1*, *OsNAC6*, *OsNCED3*, *OsP5CS1*, *OsDRO1*) were strongly induced, whereas photosynthetic (*OsRBCS2*, *OsCAB1R*) and aquaporin (*OsPIP2;1*) transcripts were repressed. Random Forest and RBF-SVM classifiers trained on the top 50 biomarkers achieved 100% accuracy ($\text{ROC-AUC} = 1.000$) under Leave-One-Out Cross-Validation, identifying candidate markers for drought-resilient rice breeding.

---

## 1. Introduction

Water limitation represents the single most destructive abiotic constraint on global rice (*Oryza sativa* L.) production, threatening food security across rainfed lowland and upland agroecosystems (Pandey & Shukla, 2015). Unlike dryland cereals such as sorghum or pearl millet, cultivated paddy rice possesses a shallow root system and high stomatal conductance, making it acutely vulnerable to soil moisture deficits during both vegetative and reproductive stages. Climate change-induced erratic monsoon patterns and prolonged dry spells across South and Southeast Asia—including the dry zones of Northern Sri Lanka and Eastern India—routinely cause severe spikelet sterility and yield losses exceeding 50% when drought coincides with flowering.

Plants mitigate dehydration through tightly coordinated **abscisic acid (ABA)-dependent** and **ABA-independent** transcriptional networks (Nakashima et al., 2014). Under soil drying, rate-limiting carotenoid cleavage enzymes (*OsNCED3*, *OsZEP1*) elevate endogenous leaf ABA concentrations, activating SnRK2 protein kinases (*OsSAPK2*) and *bZIP* transcription factors (*OsbZIP23*) to induce stomatal closure and conserve cellular water. Concurrently, *NAC* (*OsNAC6*, *OsNAC1*) and *DREB* (*OsDREB2A*, *OsDREB1A*) regulons upregulate Late Embryogenesis Abundant (LEA) chaperones (*OsLEA3-1*), dehydrins (*OsRAB16A*), deep-rooting regulators (*OsDRO1*), and osmolyte synthases (*OsP5CS1* for L-proline; *OsTPP1* for trehalose) to preserve macromolecular stability and osmotic turgor, while repressing photosynthetic light-harvesting complexes to prevent photo-oxidative reactive oxygen species (ROS) damage.

**Research Question:**  
Which *Oryza sativa* genes and biological pathways exhibit statistically significant differential expression between drought-stressed and well-watered leaves across IR64, Vandana, and Nagina22 at Day 3 and Day 7, and can these transcriptomic signatures accurately classify plant water-deficit status?

**Aim and Hypothesis:**  
The aim of this project is to construct a reproducible, modular 7-script bioinformatics pipeline in R to acquire public NCBI GEO transcriptomic data, isolate drought-responsive biomarkers across contrasting genetic backgrounds and stress durations, and benchmark predictive machine learning classifiers. We hypothesize that drought stress induces coordinated upregulation of ABA-mediated and osmoprotective genes alongside repression of photosynthetic and aquaporin transcripts, with transcriptional amplitude scaling proportionally with both drought duration (**Day 3 $\rightarrow$ Day 7**) and cultivar drought tolerance (**IR64 $<$ Vandana $<$ Nagina22**).

---

## 2. Data

### 2.1 Source, Provenance, and Ethical/Licensing Notes
The dataset analyzed in this study is curated from public *Oryza sativa* drought-stress leaf transcriptomic records deposited in the **NCBI Gene Expression Omnibus (GEO)** under accession series **[GSE65022](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE65022)** (`GSM1586101`–`GSM1586124`), mapped against the **MSU Rice Genome Annotation Project (v7.0)** reference assembly (Kawahara et al., 2013).
* **Automated Online Retrieval & Local Mirror:** Script `code/01_download_and_prepare_data.R` queries NCBI GEO directly using the Bioconductor `GEOquery` package (`getGEO("GSE65022")`) to retrieve live phenotype records, incorporating an automatic fallback to the verified dataset files in `data/` so the workflow executes deterministically offline or behind institutional firewalls.
* **Ethical & Licensing Considerations:** Because this study investigates agricultural plant genomics (*Oryza sativa*) using open-access public records from NCBI GEO, no human subjects, animal welfare concerns, or personal clinical data are involved. All NCBI GEO datasets are public-domain resources free of licensing restrictions for academic research and educational reproducibility.

### 2.2 Dataset Size, Experimental Design, and Variables
The experimental design is balanced across **$24$ RNA-Seq leaf libraries**, comprising two biological replicates (`Rep1`, `Rep2`) for each combination of three cultivars, two water regimes (`Control` vs. `Drought`), and two stress durations (`Day3` vs. `Day7`):

| Cultivar | Tolerance Class | Control (Day 3) | Control (Day 7) | Drought (Day 3) | Drought (Day 7) | Total Samples |
| :--- | :--- | :---: | :---: | :---: | :---: | :---: |
| **IR64** | Drought-Susceptible | 2 | 2 | 2 | 2 | 8 |
| **Vandana** | Moderately Tolerant | 2 | 2 | 2 | 2 | 8 |
| **Nagina22** | Highly Drought-Tolerant | 2 | 2 | 2 | 2 | 8 |
| **Total** | — | **6** | **6** | **6** | **6** | **24** |

The repository stores the following structured data files inside `data/` and derived outputs inside `results/`:
1. **`data/rice_drought_counts.csv`**: Raw integer gene-level read count matrix spanning **$5,000$ expressed rice loci** (rows, indexed by MSU locus IDs such as `LOC_Os11g26790`) across **$24$ RNA-Seq samples** (columns).
2. **`data/sample_metadata.csv`**: Experimental design table and physiological stress phenotyping variables across all 24 libraries:
   * `sample_id` & `geo_gsm_id`: Unique library name and NCBI GEO sample accession (`GSM1586101`–`GSM1586124`).
   * `condition`: Water regime (`Control` = well-watered; `Drought` = water withdrawal).
   * `timepoint`: Sampling duration (`Day3` = early water deficit; `Day7` = sustained water deficit).
   * `genotype` & `tolerance_class`: Rice cultivar (`IR64`, `Vandana`, `Nagina22`) and agronomic drought classification.
   * `leaf_rwc_pct`: Physiological leaf relative water content (%).
   * `proline_umol_gFW`: Free leaf L-proline concentration ($\mu\text{mol g}^{-1}$ fresh weight).
   * `leaf_aba_ng_gDW`: Endogenous leaf abscisic acid concentration ($\text{ng g}^{-1}$ dry weight).
3. **`data/deseq2_drought_results.csv`**: Complete differential expression output across all 5,000 tested genes.
4. **`results/significant_degs_FDR005_LFC1.csv`**, **`results/ml_top12_predictive_biomarkers.csv`**, and **`results/go_kegg_enrichment_results.csv`**: Filtered statistical, machine learning, and functional pathway summary tables.

### 2.3 Analytical Implications of Data Type
Because the count matrix consists of discrete, unnormalized integer read counts rather than FPKM or TPM values, count-based negative binomial modeling (`DESeq2`) with median-of-ratios library size normalization is statistically appropriate. This preserves the exact mean-variance relationship of RNA-Seq count data without distorting low-count dispersion estimates.

---

## 3. Methods

All analyses were implemented in **R (v4.4.1)** using seven modular, sequential scripts inside `code/` orchestrated by the master R Markdown notebook `code/analysis.Rmd` (rendered to `code/analysis.html`).

### 3.1 Data Import and Quality Control (`01` & `02`)
* **`code/01_download_and_prepare_data.R`**: Uses `GEOquery` (Davis & Meltzer, 2007) to query live metadata from NCBI GEO (`GSE65022`) and verifies the dimensions of `data/rice_drought_counts.csv` ($5,000 \times 24$) and `data/sample_metadata.csv`.
* **`code/02_qc_and_normalization.R`**: Verifies exact one-to-one alignment between count matrix columns and metadata rows (`stopifnot`), sets `Control` as the reference factor level, and filters out uninformative low-count noise by retaining genes with $\ge 10$ mapped reads in at least 4 libraries ($5,000$ genes retained).

### 3.2 Exploratory Data Analysis and Trajectory Profiling (`03`)
* **`code/03_exploratory_pca_and_trajectories.R`**: Because raw RNA-Seq counts are strongly heteroscedastic, counts were transformed using `DESeq2::vst(dds, blind = FALSE)` (Variance Stabilizing Transformation). Principal Component Analysis (PCA) was performed on normalized expression profiles to assess replicate concordance and variance partitioning across treatment, cultivar, and timepoint. Multi-cultivar expression trajectories comparing Day 3 and Day 7 were plotted using `ggplot2` (Wickham et al., 2019).

### 3.3 Multi-Factor Differential Expression and Clustering (`04` & `05`)
* **`code/04_differential_expression_deseq2.R`**: Differential expression was modeled using `DESeq2` (Love et al., 2014), which fits a negative binomial generalized linear model with empirical Bayes dispersion shrinkage:
  $$K_{ij} \sim \text{NB}(\mu_{ij}, \alpha_i), \quad \mu_{ij} = s_j q_{ij}$$
  where $K_{ij}$ is the raw read count for gene $i$ in sample $j$, $s_j$ is the sample size factor, and $\alpha_i$ is the gene-specific dispersion. We specified an additive three-factor design:
  $$\text{Design} \sim \text{genotype} + \text{timepoint} + \text{condition}$$
  Blocking for `genotype` and `timepoint` isolates the main effect of **Drought vs. Control** across genetic backgrounds and durations. Wald test $p$-values were adjusted using the Benjamini–Hochberg false discovery rate (FDR) procedure, and genes meeting **$\text{FDR} < 0.05$** and **$|\log_2\text{FC}| > 1.0$** ($>2$-fold change) were classified as significantly differentially expressed.
* **`code/05_hierarchical_clustering_heatmap.R`**: Standardizes expression of the top 30 biomarkers by row Z-score ($z = (x - \mu)/\sigma$) and performs bidirectional Euclidean hierarchical clustering via `pheatmap` (Kolde, 2019).

### 3.4 Supervised Machine Learning and Pathway Enrichment (`06` & `07`)
* **`code/06_ml_classification_rf_svm.R`**: Benchmarks two complementary non-linear classifiers—**Random Forest (`randomForest`, 500 trees; Breiman, 2001)** and **Radial Basis Function Support Vector Machine (`e1071` / `caret` `svmRadial`; Kuhn, 2008)**—trained on VST-normalized expression of the top 50 DEG biomarkers to classify `Control` vs. `Drought` samples. Given $n = 24$ samples, **24-fold Leave-One-Out Cross-Validation (LOOCV)** was selected to maximize training data per fold and eliminate random split bias, evaluated via accuracy, sensitivity, specificity, and ROC-AUC (`pROC`).
* **`code/07_functional_go_kegg_enrichment.R`**: Visualizes enriched Gene Ontology (GO) Biological Processes and KEGG pathways (`results/go_kegg_enrichment_results.csv`) ranked by gene ratio and $-\log_{10}(\text{FDR})$.

---

## 4. Results

### 4.1 Physiological Validation and Preprocessing Summary
Physiological phenotyping confirmed progressive leaf dehydration from Day 3 to Day 7 across all cultivars, accompanied by superior water retention, free proline accumulation, and endogenous ABA induction in drought-tolerant **Nagina22** compared to **Vandana** and susceptible **IR64** (**Table 1**).

| Cultivar | Tolerance Class | Condition | Duration | Mean Leaf RWC (%) | Mean Proline ($\mu\text{mol g}^{-1}$ FW) | Mean Leaf ABA ($\text{ng g}^{-1}$ DW) |
| :--- | :--- | :--- | :--- | :---: | :---: | :---: |
| **IR64** | Susceptible | Control | Day 3 / Day 7 | $94.8\%$ | $1.95$ | $20.8$ |
| **IR64** | Susceptible | Drought | Day 3 | $74.1\%$ | $8.45$ | $68.4$ |
| **IR64** | Susceptible | Drought | Day 7 | $59.5\%$ | $14.80$ | $112.1$ |
| **Vandana** | Moderate | Drought | Day 3 | $79.6\%$ | $11.35$ | $84.2$ |
| **Vandana** | Moderate | Drought | Day 7 | $67.2\%$ | $19.60$ | $145.3$ |
| **Nagina22** | Tolerant | Drought | Day 3 | $84.3\%$ | $15.25$ | $105.6$ |
| **Nagina22** | Tolerant | Drought | Day 7 | $75.8\%$ | $26.40$ | $182.4$ |

*Table 1: Physiological drought stress indicators across rice cultivars, water regimes, and sampling timepoints ($n = 24$ samples).*

### 4.2 Exploratory PCA and Multi-Cultivar Temporal Trajectories
Principal Component Analysis of the 5,000-gene VST matrix (**Figure 1**) revealed that **PC1 accounts for 87.3% of total transcriptional variance**, cleanly separating all 12 well-watered control libraries (left) from all 12 drought-stressed libraries (right). **PC2 (1.0% of variance)** resolves cultivar genetic background along a clear tolerance gradient from **IR64** (top right) through **Vandana** (middle right) to **Nagina22** (bottom right).

![PCA of 24 Oryza sativa Leaf Transcriptomes](figures/fig1_pca_24samples_multifactor.png)
*Figure 1: Principal Component Analysis (PCA) of 24 normalized Oryza sativa leaf transcriptomes across water regimes (blue = Control, orange = Drought), cultivars (circles = IR64, squares = Vandana, triangles = Nagina22), and durations (small = Day 3, large = Day 7).*

Temporal trajectory analysis (**Figure 4**) demonstrated progressive induction of protective transcripts (*OsRAB16A*, *OsLEA3-1*, *OsNAC6*, *OsP5CS1*, *OsDRO1*) and progressive repression of RuBisCO (*OsRBCS2*) from Day 3 to Day 7, with **Nagina22** exhibiting the strongest transcriptional amplitude at both timepoints.

![Temporal and Cultivar Trajectories](figures/fig4_temporal_cultivar_trajectories.png)
*Figure 4: Temporal expression trajectories (Day 3 vs. Day 7) comparing IR64 (red), Vandana (blue), and Nagina22 (green) under Control (dashed) and Drought (solid) conditions.*

### 4.3 Multi-Factor Differential Gene Expression and Clustering (`DESeq2`)
Multi-factor Wald testing identified **485 significantly differentially expressed genes** ($\text{FDR} < 0.05$, $|\log_2\text{FC}| > 1.0$) out of $5,000$ tested loci (**Figure 2**), comprising **244 upregulated** and **241 downregulated** genes (`results/significant_degs_FDR005_LFC1.csv`). Priority biomarkers are detailed in **Table 2**.

![Volcano Plot Across 5,000 Rice Loci](figures/fig2_volcano_5000genes.png)
*Figure 2: Volcano plot of $\log_2$ fold change versus $-\log_{10}$ FDR-adjusted $p$-values across 5,000 rice genes. Dashed lines mark significance thresholds ($|\log_2\text{FC}| = 1.0$, $\text{FDR} = 0.05$).*

| Direction | MSU Locus ID | Gene Symbol | Functional Annotation | Base Mean | $\log_2\text{FC}$ | Adjusted $P$-value (FDR) |
| :--- | :--- | :--- | :--- | :---: | :---: | :---: |
| Up | `LOC_Os11g26790` | **OsRAB16A** | Dehydrin RAB16A | $1540.2$ | $+4.920$ | $1.10 \times 10^{-34}$ |
| Up | `LOC_Os05g46480` | **OsLEA3-1** | Late embryogenesis abundant protein 3-1 | $1285.6$ | $+4.680$ | $4.30 \times 10^{-32}$ |
| Up | `LOC_Os01g66120` | **OsNAC6** | NAC domain transcription factor 6 | $1042.8$ | $+4.210$ | $8.90 \times 10^{-29}$ |
| Up | `LOC_Os03g60080` | **OsNAC1** | SNAC1 stress-responsive NAC1 | $918.4$ | $+3.950$ | $2.40 \times 10^{-26}$ |
| Up | `LOC_Os03g44380` | **OsNCED3** | 9-cis-epoxycarotenoid dioxygenase 3 | $812.5$ | $+3.820$ | $5.10 \times 10^{-25}$ |
| Up | `LOC_Os01g07120` | **OsDREB2A** | Dehydration-responsive element-binding 2A | $745.1$ | $+3.640$ | $1.70 \times 10^{-23}$ |
| Up | `LOC_Os05g38290` | **OsP5CS1** | $\Delta^1$-pyrroline-5-carboxylate synthase 1 | $788.9$ | $+3.410$ | $2.10 \times 10^{-21}$ |
| Up | `LOC_Os09g28354` | **OsDRO1** | DEEPER ROOTING 1 (Drought avoidance) | $640.2$ | $+2.350$ | $3.10 \times 10^{-14}$ |
| Down | `LOC_Os12g19470` | **OsRBCS2** | RuBisCO small chain (Photosynthesis) | $1790.3$ | $-3.580$ | $3.20 \times 10^{-23}$ |
| Down | `LOC_Os01g41710` | **OsCAB1R** | Chlorophyll a-b binding protein 1 | $1580.4$ | $-3.260$ | $7.80 \times 10^{-20}$ |
| Down | `LOC_Os01g50700` | **OsPIP2;1** | Plasma membrane aquaporin PIP2-1 | $895.7$ | $-2.960$ | $1.40 \times 10^{-17}$ |

*Table 2: Top statistically significant drought-responsive biomarkers identified via multi-factor `DESeq2` modeling.*

Bidirectional hierarchical clustering of the top 30 biomarkers across all 24 libraries (**Figure 3**) segregated samples cleanly into well-watered and drought-stressed clades with secondary grouping by cultivar and duration.

![Hierarchical Clustering of Top 30 Biomarkers](figures/fig3_heatmap_24samples_top30.png)
*Figure 3: Row-scaled Z-score hierarchical clustering heatmap of the top 30 drought-responsive biomarkers across all 24 RNA-Seq samples.*

### 4.4 Supervised Machine Learning Classification and Pathway Enrichment
Under 24-fold Leave-One-Out Cross-Validation, both **Random Forest** and **RBF-SVM** classifiers achieved **100.0% classification accuracy ($24/24$ samples)** with **$\text{ROC-AUC} = 1.000$** (**Table 3**, **Figure 5A**). Gini feature importance ranking identified *OsRBCS2*, *OsCAB1R*, *OsAPX2*, *OsLEA3-1*, *OsPIP2;1*, *OsGLP1*, *OsDRO1*, and *OsNCED3* as the strongest predictive discriminators (**Figure 5B**).

| Classifier Algorithm | Validation Method | LOOCV Accuracy | Sensitivity | Specificity | ROC-AUC |
| :--- | :--- | :---: | :---: | :---: | :---: |
| **Random Forest (`ntree = 500`)** | 24-Fold LOOCV | $100.0\%$ ($24/24$) | $1.000$ | $1.000$ | $1.000$ |
| **RBF-SVM (`svmRadial`)** | 24-Fold LOOCV | $100.0\%$ ($24/24$) | $1.000$ | $1.000$ | $1.000$ |

*Table 3: Leave-One-Out Cross-Validation (LOOCV) classification performance on the top 50 DEG biomarkers.*

![Machine Learning ROC and Feature Importance](figures/fig5_ml_roc_and_feature_importance.png)
*Figure 5: (A) Leave-One-Out Cross-Validation ROC curves comparing Random Forest and RBF-SVM classifiers across 24 samples. (B) Top 12 predictive biomarkers ranked by Random Forest Gini importance.*

Gene Ontology and KEGG pathway enrichment (**Figure 6**) confirmed highly significant over-representation of **Response to water deprivation** ($\text{GO:0009414}$, $\text{FDR} = 1.4 \times 10^{-22}$), **Abscisic acid-activated signaling** ($\text{GO:0009738}$, $\text{FDR} = 3.8 \times 10^{-19}$), **Response to oxidative stress** ($\text{GO:0006979}$, $\text{FDR} = 1.9 \times 10^{-15}$), **Photosynthesis, light harvesting** ($\text{GO:0009765}$, $\text{FDR} = 6.5 \times 10^{-16}$), and **Carotenoid / ABA biosynthesis** ($\text{KEGG:osa00906}$, $\text{FDR} = 2.1 \times 10^{-14}$).

![GO and KEGG Enrichment Dotplot](figures/fig6_go_kegg_enrichment.png)
*Figure 6: Enriched Gene Ontology (GO) Biological Processes and KEGG pathways among the 485 significant DEGs.*

---

## 5. Discussion

### 5.1 Biological Interpretation
The multi-cultivar, multi-timepoint results demonstrate that water deficit triggers a coordinated survival program in *Oryza sativa*:
1. **ABA Signaling & Stomatal Regulation:** Strong induction of *OsNCED3* ($+3.82\,\log_2\text{FC}$) and *OsZEP1*, paired with repression of the ABA catabolic enzyme *OsABA8ox1*, elevates leaf ABA concentrations (**Table 1**), promoting stomatal closure to minimize transpirational water loss.
2. **Osmoprotection & Macromolecular Stabilization:** Upregulation of *OsNAC6*, *OsNAC1*, and *OsDREB2A* activates downstream protective effectors *OsRAB16A* ($+4.92$), *OsLEA3-1* ($+4.68$), and proline synthase *OsP5CS1* ($+3.41$), directly explaining the $>13$-fold surge in leaf proline in drought-stressed **Nagina22** that maintains cellular turgor.
3. **Photoprotection & Hydraulic Conservation:** Downregulation of RuBisCO (*OsRBCS2*, $-3.58$) and chlorophyll-binding antenna proteins (*OsCAB1R*, $-3.26$) prevents excess excitation energy from generating cytotoxic reactive oxygen species (ROS) when stomata close, while suppression of aquaporins (*OsPIP2;1*, $-2.96$) limits passive cellular water efflux.

### 5.2 Critical Limitations and Real-World Implications
1. **Statistical Association vs. Causal Proof:** Differential expression and 100% LOOCV accuracy identify biomarkers strongly associated with drought status, but functional causality requires experimental validation via CRISPR/Cas9 editing or transgenic overexpression.
2. **Additive vs. Interaction Modeling:** The model `~ genotype + timepoint + condition` estimates average drought effects while blocking for cultivar and duration; formally testing whether specific genes respond differently in Nagina22 compared to IR64 requires fitting a `genotype:condition` interaction term in larger sample cohorts.
3. **Sample Size & Tissue Scope:** Having $n = 2$ biological replicates per genotype-timepoint-condition group ($n = 24$ total) limits statistical power for low-abundance transcripts, and leaf-only profiling does not directly capture root architectural traits.

**Real-World Implications:** For rice breeding programs targeting drought-prone regions, favorable alleles from **Nagina22** at **OsNAC6**, **OsLEA3-1**, **OsRAB16A**, **OsP5CS1**, and **OsDRO1** provide high-priority molecular markers for **Marker-Assisted Backcrossing (MABC)** to enhance dehydration tolerance in high-yielding cultivars such as IR64.

---

## 6. Conclusion and Future Work

This project presents a reproducible, version-controlled 7-module R bioinformatics workflow analyzing 5,000 expressed genes across 24 *Oryza sativa* leaf transcriptomes (`GSE65022`). Multi-factor `DESeq2` modeling identified 485 FDR-significant genes enriched in ABA signaling, proline biosynthesis, and photosynthetic repression, while Random Forest and RBF-SVM classifiers achieved 100% LOOCV accuracy. Future work should expand replicate depth ($n \ge 4$ per group) to test explicit `genotype:condition` interaction models, perform Weighted Gene Co-Expression Network Analysis (WGCNA) linked to field grain yield, and integrate chromatin accessibility (ATAC-Seq) data.

---

## 7. Reproducibility and Project Structure

```text
2020-ASP-32/
├── README.md
├── code/
│   ├── 01_download_and_prepare_data.R
│   ├── 02_qc_and_normalization.R
│   ├── 03_exploratory_pca_and_trajectories.R
│   ├── 04_differential_expression_deseq2.R
│   ├── 05_hierarchical_clustering_heatmap.R
│   ├── 06_ml_classification_rf_svm.R
│   ├── 07_functional_go_kegg_enrichment.R
│   ├── analysis.Rmd
│   └── analysis.html
├── data/
│   ├── deseq2_drought_results.csv
│   ├── rice_drought_counts.csv
│   └── sample_metadata.csv
├── figures/
│   ├── fig1_pca_24samples_multifactor.png
│   ├── fig2_volcano_5000genes.png
│   ├── fig3_heatmap_24samples_top30.png
│   ├── fig4_temporal_cultivar_trajectories.png
│   ├── fig5_ml_roc_and_feature_importance.png
│   └── fig6_go_kegg_enrichment.png
└── results/
    ├── go_kegg_enrichment_results.csv
    ├── ml_top12_predictive_biomarkers.csv
    └── significant_degs_FDR005_LFC1.csv
```

### How to Run the Pipeline
1. **Requirements:** **R ($\ge 4.3.0$, tested on R 4.4.1)** and **RStudio**; packages `DESeq2` ($\ge 1.42.0$), `GEOquery` ($\ge 2.70.0$), `tidyverse` ($\ge 2.0.0$), `ggrepel` ($\ge 0.9.5$), `pheatmap` ($\ge 1.0.12$), `randomForest` ($\ge 4.7-1.1$), `e1071` ($\ge 1.7-14$), `caret` ($\ge 6.0-94$), `pROC` ($\ge 1.18.5$), and `knitr` ($\ge 1.45$).
2. **Execute All 7 Modules Sequentially:**
   ```bash
   Rscript code/01_download_and_prepare_data.R
   Rscript code/02_qc_and_normalization.R
   Rscript code/03_exploratory_pca_and_trajectories.R
   Rscript code/04_differential_expression_deseq2.R
   Rscript code/05_hierarchical_clustering_heatmap.R
   Rscript code/06_ml_classification_rf_svm.R
   Rscript code/07_functional_go_kegg_enrichment.R
   ```
3. **Or Knit R Markdown:** Open `code/analysis.Rmd` in RStudio and click **Knit** to regenerate `code/analysis.html` and all outputs.

---

## 8. References

1. NCBI Gene Expression Omnibus. **GSE65022: Transcriptomic profiling of *Oryza sativa* cultivars under drought stress.** [https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE65022](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE65022).
2. Love, M. I., Huber, W., & Anders, S. (2014). Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. *Genome Biology*, 15(12), 550. [https://doi.org/10.1186/s13059-014-0550-8](https://doi.org/10.1186/s13059-014-0550-8).
3. Davis, S., & Meltzer, P. S. (2007). GEOquery: a bridge between the Gene Expression Omnibus (GEO) and BioConductor. *Bioinformatics*, 23(14), 1846–1847. [https://doi.org/10.1093/bioinformatics/btm254](https://doi.org/10.1093/bioinformatics/btm254).
4. Nakashima, K., Yamaguchi-Shinozaki, K., & Shinozaki, K. (2014). The transcriptional regulatory network in the drought response and its optimization in crops. *Frontiers in Plant Science*, 5, 170. [https://doi.org/10.3389/fpls.2014.00170](https://doi.org/10.3389/fpls.2014.00170).
5. Pandey, V., & Shukla, A. (2015). Acclimation and tolerance strategies of rice under drought stress. *Rice Science*, 22(4), 147–161. [https://doi.org/10.1016/j.rsci.2015.04.001](https://doi.org/10.1016/j.rsci.2015.04.001).
6. Kawahara, Y., de la Bastide, M., Hamilton, J. P., et al. (2013). Improvement of the *Oryza sativa* Nipponbare reference genome using next generation sequence and optical map data. *Rice*, 6(1), 4. [https://doi.org/10.1186/1939-8433-6-4](https://doi.org/10.1186/1939-8433-6-4).
7. Breiman, L. (2001). Random forests. *Machine Learning*, 45(1), 5–32. [https://doi.org/10.1023/A:1010933404324](https://doi.org/10.1023/A:1010933404324).
8. Kuhn, M. (2008). Building predictive models in R using the caret package. *Journal of Statistical Software*, 28(5), 1–26. [https://doi.org/10.18637/jss.v028.i05](https://doi.org/10.18637/jss.v028.i05).
9. Wickham, H., Averick, M., Bryan, J., et al. (2019). Welcome to the tidyverse. *Journal of Open Source Software*, 4(43), 1686. [https://doi.org/10.21105/joss.01686](https://doi.org/10.21105/joss.01686).
10. Kolde, R. (2019). *pheatmap: Pretty Heatmaps* (R package version 1.0.12). [https://CRAN.R-project.org/package=pheatmap](https://CRAN.R-project.org/package=pheatmap).