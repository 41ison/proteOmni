# proteOmni

**proteOmni** is a comprehensive, unified Shiny-based dashboard for visual quality control (QC), diagnostics, and differential abundance analysis of proteomics results from multiple search engines and acquisition strategies. It centralizes eight specialized modules covering DDA, DIA, and *de novo* sequencing workflows into a single interactive application.

<p align="center">
<img src="https://github.com/41ison/proteOmni/blob/main/images/proteomni_interface.png" width="700">
</p>

---

## Table of Contents

1. [Benefits](#benefits-of-using-proteomni)
2. [System Requirements](#system-requirements)
3. [Installation & Launch](#installation--launch)
4. [Project Structure](#project-structure)
5. [Modules](#modules-included)
   - [PSManalyst](#1-psmanalyst--fragpipe--dda)
   - [QC4DIANN](#2-qc4diann--dia-nn--dia)
   - [PwrQuant](#3-pwrquant--limma--stats)
   - [Casanovo](#4-casanovo--de-novo)
   - [InstaNovo](#5-instanovo--de-novo)
   - [EncyclopeDIA](#6-encyclopedia--encyclopedia--dia)
   - [Sage](#7-sage--sage--ddadia)
   - [MaxQuant](#8-maxquant--maxquant--dda)
6. [Input File Reference](#input-file-reference)
7. [PwrQuant — Analytical Pipeline](#pwrquant--analytical-pipeline)
8. [Session Logging](#session-logging)
9. [Troubleshooting](#troubleshooting)
10. [Citation](#citation)

---

## Benefits of using proteOmni

- **Unified Interface** — One central hub for evaluating proteomics results across Data-Dependent Acquisition (DDA), Data-Independent Acquisition (DIA), and *de novo* sequencing — from FragPipe, DIA-NN, MaxQuant, Sage, EncyclopeDIA, Casanovo, and InstaNovo.
- **Deep QC Insights** — Generate detailed metrics including protease fingerprints, sequence logos, mass accuracy (ppm), retention time prediction errors, charge state and peptide length distributions, missed cleavages, GRAVY index, and isoelectric point (pI) profiles.
- **Interactive Visualization** — Explore 3D QuantUMS score distributions, interactive PCA plots, sample correlation matrices, cosine/Euclidean/Jaccard similarity heatmaps, and annotated MS/MS fragmentation spectra directly in the browser.
- **Peptide-to-Protein Mapping** — Map identified peptides onto user-provided FASTA sequences with a colour-coded protein sequence viewer.
- **Differential Abundance Analysis** — Full limma-based workflow with normalization, batch correction (ComBat/SVA), flexible missing value imputation (KNN, MinProb, or missForest), MA plots, volcano plots, and statistical power mapping.
- **Functional Enrichment** — Over-representation analysis (ORA) with gene-set annotation using clusterProfiler and enrichplot.
- **Publication-Ready Output** — Download filtered result matrices and universally formatted plots (PNG/ZIP) ready for reporting and publication.
- **Session Logging** — Automatic session and console logging with a one-click log download for full reproducibility.

<p align="center">
<img src="https://github.com/41ison/proteOmni/blob/main/images/video.gif" width="700">
</p>

---

## System Requirements

| Requirement | Minimum | Recommended |
|---|---|---|
| **R** | 4.2.0 | 4.5.0+ |
| **OS** | Windows 10, macOS 12, Ubuntu 20.04 | Latest stable release |
| **RAM** | 8 GB | 16 GB+ (large DIA datasets) |
| **Browser** | Any modern browser | Chrome / Firefox |

### R package dependencies

proteOmni auto-installs all required packages on first launch. The key dependencies are:

**CRAN packages:** `shiny`, `shinydashboard`, `shinyjs`, `fresh`, `tidyverse`, `tidytext`, `janitor`, `ggpointdensity`, `ggtext`, `ggrepel`, `ggseqlogo`, `lsa`, `vegan`, `plotly`, `viridis`, `ggfortify`, `seqinr`, `zip`, `DT`, `colourpicker`, `R6`, `gridExtra`, `scales`, `lavaan`, `naniar`, `patchwork`, `pwr`, `missForest`, `data.table`, `GGally`, `arrow`

**Bioconductor packages:** `limma`, `Biostrings`, `sva`, `impute`, `ComplexHeatmap`, `clusterProfiler`, `GO.db`, `enrichplot`

**GitHub packages:** `diann` ([vdemichev/diann-rpackage](https://github.com/vdemichev/diann-rpackage))

---

## Installation & Launch

### Option 1 — Launcher scripts (recommended)

**macOS / Linux**

```bash
# Make executable (first time only)
chmod +x proteOmni_MacOS.command

# Double-click the file, or run from terminal:
./proteOmni_MacOS.command
```

**Windows**

Double-click `proteOmni_Windows.bat`. R must be on your system `PATH` (see [Troubleshooting](#troubleshooting) if you get a `Rscript not recognized` error).

### Option 2 — From R / RStudio / Positron

```r
setwd("/path/to/proteOmni/app")
shiny::runApp("proteOmni.R", launch.browser = TRUE)
```

### Option 3 — From the terminal

```bash
cd /path/to/proteOmni/app
Rscript -e "shiny::runApp('proteOmni.R', launch.browser = TRUE)"
```

> **First launch note:** On the first run, `proteOmni.R` will automatically check for and install any missing packages. This may take several minutes. Subsequent launches are fast.

---

## Project Structure

```
proteOmni/
├── proteOmni.R                  # Main app entry point (UI + server + package bootstrap)
├── proteOmni_MacOS.command      # One-click launcher for macOS / Linux
├── proteOmni_Windows.bat        # One-click launcher for Windows
├── modules/
│   ├── mod_PSManalyst.r         # PSManalyst module (FragPipe / DDA)
│   ├── mod_QC4DIANN.r           # QC4DIANN module (DIA-NN / DIA)
│   ├── mod_PwrQuant.r           # PwrQuant module (limma / differential abundance)
│   ├── dash_deNovo.r            # Casanovo de novo module
│   ├── mod_InstaNovo.r          # InstaNovo de novo module
│   ├── mod_EncyclopeDIA.r       # EncyclopeDIA module
│   ├── mod_Sage.r               # Sage module
│   ├── mod_MaxQuantMSMS.r       # MaxQuant module
│   ├── utils_fasta.r            # Shared FASTA parsing utilities
│   └── mod_TEMPLATE.r           # Template for adding new modules
├── www/
│   └── favicon.svg              # App favicon
├── images/                      # Screenshots and GIFs for README
└── README.md
```

Each module follows the standard Shiny module pattern with three exported functions:
- `<Module>_sidebar_ui(id)` — sidebar controls
- `<Module>_body_ui(id)` — main panel tabs and plots
- `<Module>_server(id, ...)` — reactive server logic

---

## Modules Included

### 1. PSManalyst — *FragPipe / DDA*

Visual QC for FragPipe DDA results. Requires `psm.tsv` and `combined_protein.tsv`; optionally accepts a FASTA file for peptide-to-protein mapping.

**Tabs:**

| Tab | Content |
|---|---|
| **PSM Viewer** | Protease fingerprint heatmap, N/C-termini sequence logos, m/z over retention time, mass error distributions (ppm and Da), charge state and peptide length distributions, amino acid frequencies, missed cleavage analysis, and pairwise sample scatter plots |
| **MS/MS Spectrum Viewer** | Annotated b/y fragment ion spectra for any selected PSM with colour-coded ion series |
| **Protein Viewer** | Peptide sequence coverage mapped onto FASTA sequences with a colour-coded viewer; sample-to-sample comparison via cosine and Jaccard similarity matrices |

---

### 2. QC4DIANN — *DIA-NN / DIA*

Diagnostics for DIA-NN `.parquet` report files. Optionally integrates a FASTA file for sequence-level coverage.

**Tabs:**

| Tab | Content |
|---|---|
| **QC Filters & Distributions** | XIC reconstruction quality, ion density in m/z–RT space, RT prediction error, charge state and peptide length distributions, missed cleavages, FASTA sequence coverage |
| **Interactive Viewer** | Sample correlation heatmap, cosine/Euclidean/Jaccard similarity matrices, 3D QuantUMS score distribution (interactive Plotly), PCA plot, and full pairwise sample correlation matrix |

---

### 3. PwrQuant — *limma / stats*

End-to-end differential abundance and statistical power analysis pipeline. Accepts any protein abundance matrix in `.tsv` or `.csv` format (proteins × samples). See [PwrQuant — Analytical Pipeline](#pwrquant--analytical-pipeline) for full details.

**Tabs:**

| Tab | Content |
|---|---|
| **Metadata Mapping** | Editable table for assigning samples to conditions and batches |
| **Sparsity** | Missing-value heatmap (`naniar::vis_miss`) |
| **Pre-processing QC** | CV distributions per condition, mean–variance relationship (loess trend), raw and normalized abundance boxplots |
| **Differential Abundance** | MA/Bland-Altman plots, volcano plots, top-20 DAP bar mirror chart, raw p-value histograms per contrast |
| **Correlation** | Inter-contrast logFC scatter with Spearman ρ and concordant/inverse/mismatch classification |
| **Power Statistics** | Sigma vs. |logFC| reliability map — proteins above the minimum detectable fold-change at 80% power are flagged as reliable |
| **Enrichment** | ORA dotplots generated with `clusterProfiler::enrichGO`, supporting 13 organism databases |
| **UpSet Plots** | Visualize intersections of proteins across multiple sample groups |

---

### 4. Casanovo — *de novo*

Visualiser for [Casanovo](https://github.com/Noble-Lab/casanovo) *de novo* sequencing output. Loads all `.mztab` files from a user-specified directory.

**Features:** score and per-amino-acid score filtering, peptide length and score distributions, N/C-termini sequence logos, amino acid frequency heatmap.

---

### 5. InstaNovo — *de novo*

Visualiser for [InstaNovo](https://github.com/instadeepai/InstaNovo) *de novo* sequencing results. Accepts a `.csv` results file; optionally integrates a FASTA file.

**Tabs:**

| Tab | Content |
|---|---|
| **Overview** | Score distribution, peptide length, charge state, mass error in ppm, PSM retention vs. score threshold curve |
| **Peptide Analysis** | Median score by peptide length, ppm error vs. score, GRAVY hydrophobicity index, pI distribution, N/C-termini sequence logos |

---

### 6. EncyclopeDIA — *EncyclopeDIA / DIA*

Aggregates and explores EncyclopeDIA DIA results. Reads all `.txt` result files from a user-specified directory.

**Tabs:**

| Tab | Content |
|---|---|
| **Overview** | Protein and peptide identifications per file, score distribution, posterior error probability (PEP), q-value, peptide yield vs. FDR curve |
| **Peptide Properties** | Charge state distribution, post-translational modifications, peptide length, GRAVY index, pI distribution, amino acid frequencies |

---

### 7. Sage — *Sage / DDA/DIA*

QC dashboard for [Sage](https://github.com/lazear/sage) search engine results. Accepts `results.sage.tsv` or `.parquet` format.

**Tabs:**

| Tab | Content |
|---|---|
| **Overview** | PSM counts, unique proteins and peptides per file, LDA discriminant score distribution |
| **Peptide Properties** | Charge state, length density, missed cleavages, GRAVY hydrophobicity, pI distribution |
| **Mass Errors** | RT vs. mass error scatter, fragment error in Da and ppm, RT vs. precursor error, precursor mass error density |
| **Scoring & Validation** | Peptide and protein q-value distributions, peptide yield vs. FDR curve |

---

### 8. MaxQuant — *MaxQuant / DDA*

QC module for MaxQuant results. Requires `msms.txt` and `evidence.txt` from a MaxQuant output directory.

**Features:** Annotated MS/MS fragmentation spectrum viewer (b/y ions colour-coded by series) for any peptide in `msms.txt`; evidence-level QC metrics from `evidence.txt` including mass error distributions, charge states, PTM profiles, missed cleavages, and more.

---

## Input File Reference

| Module | Required files | Optional |
|---|---|---|
| **PSManalyst** | `psm.tsv`, `combined_protein.tsv` | FASTA file |
| **QC4DIANN** | `report.parquet` (DIA-NN output) | FASTA file |
| **PwrQuant** | Protein abundance matrix (`.tsv` / `.csv`, proteins × samples) | — |
| **Casanovo** | Directory path containing `.mztab` files | — |
| **InstaNovo** | InstaNovo results `.csv` file | FASTA file |
| **EncyclopeDIA** | Directory path containing EncyclopeDIA `.txt` result files | FASTA file |
| **Sage** | `results.sage.tsv` or `results.sage.parquet` | FASTA file |
| **MaxQuant** | `msms.txt`, `evidence.txt` | — |

### PwrQuant abundance matrix format

The matrix must have proteins as rows and samples as columns. The first column is used as the row identifier (protein/gene names). Values should be **raw intensities** or **log2-transformed intensities** (proteOmni applies `log2(x + 1)` internally if values appear untransformed). Example:

```
ProteinID    Sample_A1    Sample_A2    Sample_B1    Sample_B2
TP53         1.2e7        1.4e7        2.1e7        2.3e7
EGFR         NA           8.3e6        9.1e6        NA
...
```

Accepted delimiters: tab (`.tsv`, `.txt`) or comma (`.csv`). Duplicate row IDs are resolved automatically with `make.unique()`.

---

## PwrQuant — Analytical Pipeline

The PwrQuant module implements a complete quantitative proteomics workflow in the following sequential steps:

### Step 1 — Metadata assignment
Assign each sample column a **Condition** and **Batch** label using the editable table in the *Metadata Mapping* tab. Condition labels are used to build the design matrix; Batch labels are used for ComBat correction.

### Step 2 — Missing value filter
Proteins with fewer than `min_valid_pct`% of valid (non-missing) values per group are excluded before imputation and modelling. Set to 0 to disable filtering.

### Step 3 — Imputation (robust mode only)
Imputation is performed when **Limma Regression Method** is set to `robust`. Three strategies are available (selectable in the sidebar):

| Method | Speed | Missing-data model | Best for |
|---|---|---|---|
| **KNN** (default) | ⚡ ~0.2 s | MAR / MCAR — borrows information from k=5 nearest proteins | General use; moderate missingness |
| **MinProb** | ⚡⚡ ~0.03 s | MNAR — Gaussian draw at the detection limit (mean − 1.8 SD per column) | MNAR-dominated datasets, large matrices |
| **missForest** | 🐢 ~40 min | MAR / MCAR — random-forest multiple imputation | Maximum accuracy; small matrices are preferred |

> **Why missForest is slow:** after transposition, missForest receives a P-sample × N-protein matrix and builds one random forest per protein to predict missing values. For a typical proteomics dataset with 3,000 proteins across 3 groups, this means ~9,000 forests — reducing `ntree` or `maxiter` does not help because the bottleneck is the number of trees, not their depth. KNN is **~10,000× faster** and MinProb is **~73,000× faster** for equivalent datasets in our tests.

When `ls` (ordinary least squares) regression is selected, imputation is skipped entirely and the raw log2 matrix is passed directly to limma.

### Step 4 — Batch correction
If more than one unique batch label is present, `sva::ComBat` is applied to the imputed matrix using empirical Bayes priors.

### Step 5 — Normalization
Three between-array normalization methods are available via `limma::normalizeBetweenArrays`:

| Method | When to use |
|---|---|
| `cyclicloess` (default) | General purpose; robust to composition effects |
| `quantile` | When identical distributions across samples is a valid assumption |
| `scale` | Per-sample mean/variance scaling |

### Step 6 — Linear modelling and eBayes
A `~ 0 + condition` design matrix is built and contrasts are constructed from user-specified pairs (e.g. `Treatment-Control`). `limma::lmFit` is called with the selected regression method (`ls` or `robust`), followed by `limma::contrasts.fit` and `limma::eBayes`.

The **eBayes trend** parameter (sidebar toggle) controls whether the prior variance is allowed to vary with the mean log-intensity. Use the *Mean–Variance* plot in the Pre-processing QC tab to guide this choice:
- **Flat trend line** → `trend = FALSE`
- **Positively sloped trend** → `trend = TRUE`

### Step 7 — Power analysis
`pwr::pwr.t.test` is used with the median replicate count per group, α = 0.05, and power = 0.80 to compute the minimum detectable effect size. Proteins with |logFC| ≥ this threshold are flagged as **reliable**, and the power map is visualised in the *Power Statistics* tab.

### Step 8 — Results and significance calling
A protein is called **significant** if:
- `adj.P.Val ≤ 0.05`,
- `Is_reliable == TRUE` (above the power-adjusted fold-change threshold), and
- it is not imputation-driven (i.e. not completely absent in one of the contrast groups).

Proteins that are fully missing in one condition are flagged as `imputation_driven` and classified as *Not significant*, avoiding false positives driven by structural zeros.

### Step 9 — Functional enrichment (ORA)
`clusterProfiler::enrichGO` is run on the significant reliable proteins from each contrast × direction combination. Supported organism databases include Human, Mouse, Rat, Zebrafish, *C. elegans*, Drosophila, Yeast, Arabidopsis, Bovine, Canine, Pig, Chicken, and Macaque.

Protein identifiers are auto-parsed from three common formats:
1. UniProt FASTA header with `GN=SYMBOL` tag (e.g. `sp|P04637|P53_HUMAN GN=TP53`)
2. Pipe-separated UniProt ID with gene prefix (e.g. `sp|P04637|TP53_HUMAN`)
3. Plain gene symbols (e.g. `TP53`)

---

## Session Logging

Every proteOmni session automatically writes a log file (`Session_Info_log.txt`) in the working directory containing:

- Session timestamp and R version
- Platform information
- Installed package versions for all key dependencies
- All R `message()` and `warning()` calls with timestamps

The log file can be downloaded at any time using the **Download Log History** button at the top of the sidebar. This is useful for reproducing analyses and reporting issues.

---

## Troubleshooting

<details>

<summary><b>Rscript is not recognized (Windows)</b></summary>

When trying to execute the application using the `proteOmni_Windows.bat` file for the first time on Windows, you might encounter the following error:

> `'rscript' is not recognized as an internal or external command, operable program or batch file.`

This happens because Windows doesn't know where the R executable (`Rscript.exe`) is located. To fix this, you must add the R `bin` folder to your system's Environment Variables path.

**How to add R to your PATH:**

1. Open the Windows **Start Menu**, search for **"Environment Variables"**, and click on **"Edit the system environment variables"**.

   <p align="center">
   <img src="https://github.com/41ison/proteOmni/blob/main/images/paste-7.png" width="300">
   </p>

2. In the System Properties window, click the **"Environment Variables..."** button near the bottom.

   <p align="center">
   <img src="https://github.com/41ison/proteOmni/blob/main/images/paste-8.png" width="300">
   </p>

3. In the new window, find the **"Path"** variable under the *System variables* list (or *User variables*), select it, and click **"Edit..."**.

   <p align="center">
   <img src="https://github.com/41ison/proteOmni/blob/main/images/paste-9.png" width="300">
   </p>

4. Click **"New"** and paste the folder path to your R `bin` directory. This path usually looks like: `C:\Program Files\R\R-4.x.x\bin` *(Replace `4.x.x` with your specific R version)*.

   <p align="center">
   <img src="https://github.com/41ison/proteOmni/blob/main/images/paste-10.png" width="300">
   <img src="https://github.com/41ison/proteOmni/blob/main/images/paste-12.png" width="300">
   </p>

5. Click **"OK"** on all windows to save the changes.

6. Open a new Command Prompt (or just double-click the `.bat` file again) to run proteOmni successfully.

</details>

<details>

<summary><b>Error in loadNamespace (Windows / macOS)</b></summary>

When launching proteOmni via the `.bat` or `.command` file, you may see an error like:

> `Error in loadNamespace(j <- i[[1L]], c(lib.loc, .libPaths()), versionCheck = vI[[j]]) :`
> `namespace 'promises' 1.3.3 is already loaded, but >= 1.5.0 is required`

This means one or more R packages in your library are outdated. To fix this, update the affected package(s) from within R or RStudio:

```r
install.packages("promises")
```

Replace `"promises"` with whatever package name appears in your error. After installation, re-launch proteOmni. You can check the currently installed version with:

```r
packageVersion("promises")
```

</details>

<details>

<summary><b>App takes a long time to start (first launch)</b></summary>

On the first launch, `proteOmni.R` checks for and installs all missing dependencies (CRAN + Bioconductor + GitHub). This is expected and may take 5–15 minutes depending on your internet speed and machine. Subsequent launches are fast.

To pre-install all dependencies manually before the first launch, run from R:

```r
source("/path/to/proteOmni/app/proteOmni.R")
```

</details>

<details>

<summary><b>PwrQuant imputation is very slow</b></summary>

If you selected **missForest** as the imputation method, imputation can take 30–40 minutes for a typical 3,000-protein dataset with 3 conditions. This is expected — missForest builds one random forest per protein per group.

**Solution:** switch to **KNN** (default, ~10,000× faster) or **MinProb** (~73,000× faster) in the *Imputation Method* dropdown, which appears in the sidebar when `robust` regression is selected. For most proteomics datasets with MNAR-type missingness, **MinProb** is the recommended choice.

</details>

<details>

<summary><b>ORA returns no enriched terms</b></summary>

This can happen for three reasons:

1. **Protein identifiers are not gene symbols.** ORA requires gene symbols (e.g. `TP53`, `EGFR`). proteOmni attempts to parse them from UniProt FASTA headers (`GN=` tag) and pipe-separated IDs automatically, but if your matrix uses accession numbers only, no mapping will be found.

2. **Too few significant proteins.** If there are fewer than ~5 significant reliable proteins in a contrast, the hypergeometric test is underpowered. Try relaxing the minimum valid value threshold or checking whether the experiment is adequately powered.

3. **Wrong organism database.** Ensure the selected database in the sidebar matches your sample organism.

</details>

<details>

<summary><b>ComBat error / batch correction skipped</b></summary>

ComBat requires at least 2 distinct batch labels. If all samples are assigned to the same batch (default `"1"`), batch correction is automatically skipped with a notification. If ComBat crashes (e.g. due to rank-deficient design or excessive NAs), proteOmni falls back silently to the uncorrected matrix and logs the error — check the session log for details.

</details>

---

## Citation

If you use proteOmni in your research, please cite the following:

- Chaves AFA. PSManalyst: A Dashboard for Visual Quality Control of FragPipe Results. *J Proteome Res.* 2025;24(9):4344-4346. doi: [10.1021/acs.jproteome.5c00557](https://doi.org/10.1021/acs.jproteome.5c00557).
- Moschem JDC, de Barros BCSC, Serrano SMT, Chaves AFA. Decoding the Impact of Isolation Window Selection and QuantUMS Filtering in DIA-NN for DIA Quantification of Peptides and Proteins. *J Proteome Res.* 2025;24(8):3860-3873. doi: [10.1021/acs.jproteome.5c00009](https://doi.org/10.1021/acs.jproteome.5c00009).

---

<p align="center">
Made by the proteOmni team &nbsp;·&nbsp; <a href="https://github.com/41ison/proteOmni">github.com/41ison/proteOmni</a>
</p>

<p align="center">
Alison FA Chaves | Pedro G Castro | André A Tchakerian
</p>