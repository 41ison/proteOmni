[![DOI](https://zenodo.org/badge/1205443879.svg)](https://doi.org/10.5281/zenodo.21195113)
![GitHub Release](https://img.shields.io/github/v/release/41ison/proteOmni)
![GitHub License](https://img.shields.io/github/license/41ison/proteOmni)

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
- **Differential Abundance Analysis** — Full limma-based workflow with normalization, batch correction (ComBat/SVA), flexible missing value imputation (KNN, MinProb, BPCA or missForest), MA plots, volcano plots, and a simulation-based prospective power analysis.
- **Prospective Power / Sensitivity Analysis** — Minimum detectable difference (MDD) estimated by simulating fresh datasets from the fitted empirical Bayes priors, spiking in known fold changes, and re-running the real design through `lmFit` → `eBayes` → BH — optionally including intensity-dependent missingness and imputation. A replicate sweep answers "how many replicates would I need?".
- **Functional Enrichment** — GO over-representation analysis (ORA) with `clusterProfiler::enrichGO`, run per contrast and split by regulation direction, using Bioconductor OrgDb annotation packages for 20 supported organisms.
- **Protein-Protein Interaction Networks** — STRING interaction networks for the significant proteins of any contrast via the `STRINGdb` package, with nodes halo-coloured by up/down regulation and confidence-score edge filtering.
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
| **OS** | Windows 10, macOS, Linux | Latest stable release |
| **RAM** | 8 GB | 16 GB+ (large DIA datasets) |
| **Browser** | Any modern browser | Chrome / Firefox |

> [!IMPORTANT]
> If you are using Windows, you will need to install [Rtools](https://cran.r-project.org/bin/windows/Rtools/). See troubleshooting for adding [Rscript] to your Environment Variables.

### R package dependencies

proteOmni auto-installs all required packages on first launch. The key dependencies are:

**CRAN packages:** `shiny`, `shinydashboard`, `shinyjs`, `fresh`, `devtools`, `tidyverse`, `tidytext`, `janitor`, `ggpointdensity`, `ggtext`, `ggrepel`, `ggseqlogo`, `ggsci`, `lsa`, `vegan`, `plotly`, `viridis`, `RColorBrewer`, `ggfortify`, `seqinr`, `zip`, `DT`, `colourpicker`, `R6`, `gridExtra`, `scales`, `lavaan`, `naniar`, `patchwork`, `missForest`, `data.table`, `GGally`, `arrow`, `httr`, `jsonlite`, `BiocManager`

**Bioconductor packages:** `limma`, `Biostrings`, `sva`, `impute`, `pcaMethods`, `ComplexHeatmap`, `clusterProfiler`, `GO.db`, `enrichplot`, `AnnotationDbi`, `STRINGdb`

> The power analysis is implemented directly against limma's empirical Bayes objects and the non-central *t* distribution in base R, so the `pwr` package is **not** a dependency anymore. See [Step 7](#step-7--power-analysis--minimum-detectable-difference-mdd) for why a textbook power calculation does not apply here.

**OrgDb annotation packages (auto-installed on demand):** the organism selected in the PwrQuant *GO Enrichment* dropdown determines which Bioconductor `org.*` package is required (e.g. `org.Hs.eg.db` for human, `org.Mm.eg.db` for mouse). proteOmni installs the matching package via `BiocManager` the first time that organism is used. 20 organisms are supported — see [Step 9](#step-9--functional-enrichment-ora-with-clusterprofiler).

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
source("run.R")
```

### Option 3 — From the terminal

```bash
cd /path/to/proteOmni/app
Rscript run.R
```

> **First launch note:** On the first run, `run.R` sources `bootstrap.R`, which checks for and installs any missing packages (including Shiny itself) before the app starts. This may take several minutes. Subsequent launches are fast.
>
> Always start proteOmni through `run.R`. Calling `shiny::runApp("proteOmni.R")` directly fails on a fresh R installation, because R has to load the `shiny` package before it ever reads the app file. If you only want to install dependencies without launching the app, run `Rscript -e 'source("bootstrap.R")'`.

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

**Observation:** From version 2.6.0 DIA-NN have a nice `Analyze` module to help you process your data.

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
| **Power Statistics** | Prospective MDD power curve from simulated data, an MDD-vs-replicates design sweep, and a conditional per-protein sensitivity map (diagnostic only) |
| **Enrichment** | GO over-representation analysis with `clusterProfiler::enrichGO`, run per contrast and split by up/down direction; results shown as a dotplot (top terms per contrast) and a Manhattan plot (all enriched terms), across 20 supported OrgDb organisms |
| **Interaction Network** | STRING protein-protein interaction network (`STRINGdb`) for a selected contrast, with nodes halo-coloured by regulation and edges filtered by a confidence score; unmapped proteins are listed for transparency |
| **UpSet Plots** | Visualize intersections of proteins across multiple sample groups; you can download the table |

---

### 4. Casanovo — *de novo*

**Observation:** Starting from July 2026, Casanovo have a dedicated GUI built in Java language. Please, use the [CasanovoGUI app](https://github.com/Noble-Lab/CasanovoGUI).

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

### 7. Sage — *Sage / DDA*

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
| **PSManalyst** | `path/to/search_results`, `combined_protein.tsv` | FASTA file |
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
| **KNN** (default) | ~0.2 s | MAR / MCAR — borrows information from k=5 nearest proteins | General use; moderate missingness |
| **MinProb** | ~0.03 s | MNAR — Gaussian draw at the detection limit (mean − 1.8 SD per column) | MNAR-dominated datasets, large matrices |
| **missForest** | ~40 min | MAR / MCAR — random-forest multiple imputation | Maximum accuracy; small matrices are preferred |
| **bPCA** | ~2 min | MAR / MCAR - BPCA borrows signal across all samples and conditions | Flexible statistical framework |

> **Why missForest is slow:** after transposition, missForest receives a P-sample × N-protein matrix and builds one random forest per protein to predict missing values. For a typical proteomics dataset with 3,000 proteins across 3 groups, this means ~9,000 forests — reducing `ntree` or `maxiter` does not help because the bottleneck is the number of trees, not their depth. KNN is **~10,000× faster** and MinProb is **~73,000× faster** for equivalent datasets in our tests.

When `ls` (ordinary least squares) regression is selected, imputation is skipped entirely and the raw log2 matrix is passed directly to limma.

### Step 4 — Batch correction
If more than one unique batch label is present, `sva::ComBat` is applied to the imputed matrix using empirical Bayes priors.

### Step 5 — Normalization
Three between-array normalization methods are available via `limma::normalizeBetweenArrays`:

| Method | When to use |
|---|---|
| `none` | When the data behaves very well |
| `cyclicloess` (default) | General purpose; robust to composition effects |
| `quantile` | When identical distributions across samples is a valid assumption |
| `scale` | Per-sample mean/variance scaling |

### Step 6 — Linear modelling and eBayes
A `~ 0 + condition` design matrix is built and contrasts are constructed from user-specified pairs (e.g. `Treatment-Control`). `limma::lmFit` is called with the selected regression method (`ls` or `robust`), followed by `limma::contrasts.fit` and `limma::eBayes`.

The **eBayes trend** parameter (sidebar toggle) controls whether the prior variance is allowed to vary with the mean log-intensity. Use the *Mean–Variance* plot in the Pre-processing QC tab to guide this choice:
- **Flat trend line** → `trend = FALSE`
- **Positively sloped trend** → `trend = TRUE`

### Step 7 — Power analysis / minimum detectable difference (MDD)

The *Power Statistics* tab reports **two distinct quantities**. They answer different questions and should not be confused.

#### 7a — Prospective MDD by simulation (the actual power statement)

A textbook two-sample calculation (`pwr.t.test`, `power.t.test`) at α = 0.05 does not describe this pipeline, for three independent reasons:

1. **The alpha is not the operating threshold.** Significance is called on BH-adjusted p-values across thousands of proteins. At the BH boundary the largest rejected raw p-value is `FDR × R / m` — with 5,000 proteins and 1,000 rejections that is 0.01, and far smaller with few rejections. Using 0.05 is anti-conservative.
2. **The degrees of freedom are wrong.** A two-sample t-test has `2n − 2` df; the moderated statistic actually used has `2n − 2 + d₀`. Pairing a shrunken variance with unshrunken df is internally inconsistent, and conservative.
3. **It is circular.** Comparing each protein's MDD against that same protein's observed fold change is the observed-power fallacy (see 7b).

Errors (1) and (2) push in opposite directions, so a naive number can look plausible while being unjustifiable — the magnitude is not the real problem, the interpretation is.

Instead, proteOmni simulates the pipeline end-to-end:

1. Per-protein variances are drawn from the empirical Bayes prior estimated by `eBayes()` — $$\sigma^2_g \sim d_0 s_0^2 / \chi^2_{d_0}$$ (with `trend = TRUE`, `s2.prior` stays a per-protein function of abundance, preserving the mean–variance trend).
2. Fresh abundance matrices are built at the protein mean intensities, and a known two-sided fold change is spiked into a random subset of proteins (proportion π₁).
3. Optionally, values are knocked out using an **intensity-dependent (MNAR) dropout model** — a binomial GLM of missing-value counts on mean log2 abundance fitted to *your* matrix — and the same imputation used in Step 3 is applied.
4. The **real** design and contrast matrices are pushed through `lmFit` → `contrasts.fit` → `eBayes` → BH, and power is measured as the fraction of spiked proteins recovered at the same FDR cutoff used on the real data.

The **MDD** is the true log2 fold change at which the power curve crosses the target power, obtained by monotone interpolation. The realised false discovery proportion is also tracked as a sanity check on the BH machinery.

Sidebar controls (under *Prospective MDD (simulation)*; runs on demand via **Run MDD simulation**):

| Control | Default | Meaning |
|---|---|---|
| **Assumed proportion differential (π₁)** | 0.10 | Fraction of proteins truly changing — affects the BH threshold and therefore power. Use your observed hit rate as a guide; report the value, and check sensitivity at 0.01 and 0.30 if it sits near a decision boundary |
| **Target power** | 0.80 | Convention. Raise to 0.90 for confirmatory work; the MDD rises accordingly |
| **Max log₂FC tested** | 2 | Upper end of the grid. Raise if the summary reports the target power was never reached |
| **log₂FC grid step** | 0.2 | Curve resolution. The crossing region is steep, so 0.1 gives a materially tighter estimate; 0.2 is adequate for a first look |
| **Simulations per grid point** | 5 | Monte Carlo replicates: 5 for exploration, ~20 for a number you intend to publish. Cost is linear |
| **Random seed** | 4817 | Fixed for reproducibility. Change it to confirm a conclusion is not seed-specific |
| **Include missingness + imputation** | on | Leave on for realism. Turn off to isolate how much sensitivity the missing-data handling costs |

> **Double-counting caveat on the missingness arm.** The variance prior was itself estimated from already-imputed data, so the gap between the missingness and complete-data arms measures the *marginal* cost of a second imputation round — treat it as a lower bound on the true cost.

**Outputs.** A power curve per contrast with a shaded 95% Monte Carlo band (a wide band means more simulations are needed) and a dashed vertical line at the MDD, plus a summary table:

| Column | Meaning |
|---|---|
| `MDD (log2FC)` | Fold change, in log2 units, detectable at the target power |
| `MDD (fold change)` | The same number on a linear scale — the one to quote in a manuscript |
| `Mean realised FDP` | Ground-truth false discovery proportion. Because the simulation knows the truth, it can verify BH is delivering the requested FDR; a value meaningfully above your cutoff is evidence that a pipeline step (usually imputation) is breaking FDR control |
| `Note` | Flags a target power never reached, or already reached at the smallest fold change tested |

**Replicate sweep** (**Run replicate sweep**; comma-separated grid in the sidebar, minimum 2, default `3,4,5,6,8,10,12,15` — include your current *n* so the plot has a reference point) rebuilds a balanced design at each replicate count and reruns the simulation with the variance prior held fixed, answering the design question directly. Expect roughly $$1/\sqrt{n}$$ behaviour, degrading at the low end where the df penalty bites, with returns typically flattening past about eight replicates. Because only the crossing point is needed, the sweep brackets and bisects rather than scanning the grid.

Two limits: the prior is held fixed across the sweep (defensible — it is a property of the assay, not the sample size — but it was estimated at your current residual df, so the small-*n* end is optimistic about how well eBayes would shrink), and when the requested replicate count equals the min-valid floor, no dropout can be simulated, so that row shows 0% missingness by construction.

> [!TIP]
> A defensible way to report the result: *"At 80% power and a 5% BH-FDR, and assuming 10% of proteins are differentially abundant, this design resolves a 1.4-fold change (0.49 log2). The estimate comes from parametric simulation using the empirical Bayes variance prior fitted to these data, including the observed intensity-dependent missingness."*

#### 7b — Conditional per-protein sensitivity (diagnostic only)

Shown in the *Conditional sensitivity* box and exported as `Conditional_MDD_Log2FC` in the downloaded results. It fixes errors (1) and (2) above but not (3), which is why it is labelled a diagnostic.

For each protein, the log2 fold change that would have been needed to reach the target power **given that protein's posterior variance** is computed from the non-central *t* distribution at:

- the moderated degrees of freedom, `df.residual + df.prior`, and
- the per-test alpha implied by the BH criterion, `FDR × R / m` (falling back to the Bonferroni-like floor `FDR / m` when there are no rejections).

Unbalanced group sizes and non-pairwise contrasts are handled correctly because `stdev.unscaled` already encodes the design.

> **This quantity is not used to call significance.** Filtering on an observed-variance MDD is the observed-power fallacy, and is algebraically vacuous: since `MDD_g = d · s_g` and `t_g = logFC_g / (s_g · √(2/n))`, the test `|logFC| ≥ MDD` reduces to `|t| ≥ d·√(n/2)` — a fixed cutoff on the *t* statistic, with the variance cancelling entirely. It adds no power information; it is just a second, hidden significance threshold. It is shown as a diagnostic map of posterior SD vs. |logFC| only.

### Step 8 — Results and significance calling
A protein is called **significant** if:
- `adj.P.Val ≤ 0.05`,
- log2FC threshold user-defined,
- Safeguard for not imputation-driven (i.e. not completely absent in one of the contrast groups).

Proteins that are fully missing in one condition are flagged as `imputation_driven` and classified as *Not significant*, avoiding false positives driven by structural zeros. You still can extract valuable information from the proteins detected only in one condition by using the UpSet View in PwrQuant.

### Step 9 — Functional enrichment (ORA with clusterProfiler)
GO over-representation analysis is performed with `clusterProfiler::enrichGO()`. For every selected contrast, the significant proteins are split by regulation direction (**Increased** / **Decreased**) and each gene set is tested separately against the shared background (*universe*) of all quantified proteins in the matrix. Enriched terms across contrasts are stacked into a single result and rendered as:
- a **dotplot** — top terms per contrast (count on the x-axis, size = gene count, colour = adjusted p-value), and
- a **Manhattan plot** — every enriched term per contrast, height = −log10 adjusted p-value.

Annotations come from a Bioconductor OrgDb (`org.*`) package rather than a live web lookup, so no per-run network fetch is required once the package is installed. The **Organism (OrgDb)** dropdown selects the package and its NCBI taxon ID; the matching package is auto-installed via `BiocManager` on first use. 20 organisms are supported, including:

| Organism | OrgDb package | Taxon |
|---|---|---|
| *Homo sapiens* (human) | `org.Hs.eg.db` | 9606 |
| *Mus musculus* (mouse) | `org.Mm.eg.db` | 10090 |
| *Rattus norvegicus* (rat) | `org.Rn.eg.db` | 10116 |
| *Drosophila melanogaster* | `org.Dm.eg.db` | 7227 |
| *Danio rerio* (zebrafish) | `org.Dr.eg.db` | 7955 |
| *Saccharomyces cerevisiae* | `org.Sc.sgd.db` | 559292 |
| *Arabidopsis thaliana* | `org.At.tair.db` | 3702 |
| *Escherichia coli* K-12 | `org.EcK12.eg.db` | 83333 |

*(showing 8 of 20; the remaining organisms include cow, pig, dog, chicken, rhesus, chimpanzee, frog, mosquito, worm, and others.)*

Sidebar controls let you choose the contrasts to test, the **GO Category** (ALL / BP / CC / MF), the **GO term FDR cutoff**, and how many terms are shown per contrast. Protein identifiers are auto-detected as either UniProt accessions or gene symbols (setting the `enrichGO` keyType to `UNIPROT` or `SYMBOL`); a coverage notification reports how many significant IDs are annotated in the chosen OrgDb.

Protein identifiers are auto-parsed from three common formats:
1. UniProt accession, plain or pipe-separated (e.g. `P04637` or `sp|P04637|TP53_HUMAN`)
2. UniProt FASTA header with `GN=SYMBOL` tag (e.g. `sp|P04637|P53_HUMAN GN=TP53`)
3. Plain gene symbols (e.g. `TP53`)

### Step 10 — Protein-protein interaction network (STRING)
The *Interaction Network* tab maps the significant proteins of a chosen contrast onto the [STRING](https://string-db.org/) database (v12.0) with the `STRINGdb` package and renders the network via `plot_network()`. The organism reuses the taxon selected for GO enrichment in Step 9.

Controls:
- **Protein Set** — build the network from all significant proteins (*Combined*), or restrict to *Increased only* or *Decreased only*.
- **Min. combined score (0–1000)** — filters STRING edges by confidence (STRING's default is 400).

Nodes are halo-coloured **red** for increased and **blue** for decreased proteins, mirroring the regulation flagging used elsewhere in the module. Before building, proteOmni checks the selected taxon against STRING's official species list; if the organism is not covered it reports closely related covered organisms as suggestions instead of failing silently. Proteins that cannot be mapped to a STRING identifier are listed separately for transparency and excluded from the plotted network.

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

<summary><b>MDD simulation reports "target power not reached"</b></summary>

The power curve never crossed the target within the tested fold-change range. Options, in rough order of preference:

1. **Raise "Max log₂FC tested".** The default upper bound of 2 (a 4-fold change) may simply be below this design's sensitivity.
2. **Check π₁.** A very small assumed proportion of differential proteins makes the BH threshold stringent and suppresses power. If π₁ is near a decision boundary, check sensitivity at 0.01 and 0.30.
3. **Lower the target power**, or accept the honest conclusion that the design is underpowered for effects of the size you care about — the replicate sweep will tell you how many replicates would be needed.

The mirror-image note, *"target power already reached at the smallest LFC tested"*, means the MDD is below the grid floor; lower the **log₂FC grid step** to resolve it.

</details>

<details>

<summary><b>MDD simulation is slow</b></summary>

Cost scales as (grid points) × (simulations per grid point) × (contrasts). Enabling **Include missingness + imputation** makes each evaluation roughly ten times more expensive, because `lmFit` falls back to a per-row fit as soon as NAs are present, and the chosen imputer runs on every simulated matrix.

For a quick first pass, reduce **Simulations per grid point** to 3, coarsen the grid step to 0.5, and turn missingness off; then rerun at full settings for the number you intend to report. The replicate sweep already uses bisection instead of a full grid, so it is cheaper per design point than the main curve.

</details>

<details>

<summary><b>ORA returns no enriched terms</b></summary>

This can happen for a few reasons:

1. **Protein identifiers aren't annotated in the OrgDb.** `enrichGO` needs UniProt accessions or gene symbols (e.g. `TP53`, `EGFR`) that exist in the selected `org.*` package. proteOmni parses them from UniProt FASTA headers (`GN=` tag), pipe-separated IDs, or accessions automatically, but non-standard identifiers won't map. Check the annotation-coverage notification shown when ORA runs — a low percentage means most IDs weren't recognised.

2. **Wrong organism selected.** Pick the organism that matches your data in the **Organism (OrgDb)** dropdown. A mismatched OrgDb will annotate few or no proteins.

3. **Too few significant proteins.** If a contrast × direction gene set has only a handful of proteins, the over-representation test is underpowered and may return nothing. Try relaxing the minimum valid-value threshold, loosening the GO term FDR cutoff, or checking whether the experiment is adequately powered.

4. **OrgDb package couldn't be installed.** The first use of a new organism triggers a `BiocManager::install()` of its `org.*` package — check your internet connection and the session log if enrichment fails to start.

</details>

<details>

<summary><b>STRING network is empty or the organism is unsupported</b></summary>

The *Interaction Network* tab reuses the organism selected for GO enrichment. If STRING (v12.0) does not cover that taxon, proteOmni reports it explicitly and suggests closely related organisms that *are* covered — this is a coverage limitation of STRING, not an app error. Try a more commonly studied organism, or use a species-level taxon ID instead of a strain-specific one.

An empty network with a covered organism usually means too few proteins mapped to STRING identifiers (see the "could not be mapped" list under the plot) or that the **Min. combined score** threshold is too high — lower it toward STRING's default of 400 to retain more edges.

</details>

<details>

<summary><b>ComBat error / batch correction skipped</b></summary>

ComBat requires at least 2 distinct batch labels. If all samples are assigned to the same batch (default `"1"`), batch correction is automatically skipped with a notification. If ComBat crashes (e.g. due to rank-deficient design or excessive NAs), proteOmni falls back silently to the uncorrected matrix and logs the error — check the session log for details.

</details>

<details>

<summary><b>Linux: Arrow package not compiling</b></summary>

With you have problems to compile the arrow-apache in proteOmni inicialization, install through your package manager, for example:

```bash
# Arch
sudo pacman -S arrow

# Debian/Ubuntu
sudo apt install -y -V ca-certificates lsb-release wget
wget https://packages.apache.org/artifactory/arrow/$(lsb_release --id --short | tr 'A-Z' 'a-z')/apache-arrow-apt-source-latest-$(lsb_release --codename --short).deb

sudo apt install -y -V ./apache-arrow-apt-source-latest-$(lsb_release --codename --short).deb

```

 System dependecies for arrow compiling and work is: `gcc`, `curl`, `openssl`.
</details>

<details>

<summary><b>Linux: Error in utils::browseURL(appUrl)</b></summary>

The complete error is:
```bash
Error in utils::browseURL(appUrl) : 
  'browser' should not be a empty string
Calls: <Anonymous> -> .setupShinyApp -> <Anonymous>
```

R in terminal don't have a default browser for open, you have to configure manually

```bash
# In: /home/yourusername
nano .Rprofile # <- Can use vim or emacs too

# Insert this line
options(browser = "<your default browser>")
```

This line tells R which browser to open, without it, it doesn't know where to open proteOmni

</details>

---

## Citation

If you use proteOmni in your research, consider to cite the following:
- **PSManalyst:**
> Chaves AFA. PSManalyst: A Dashboard for Visual Quality Control of FragPipe Results. *J Proteome Res.* 2025;24(9):4344-4346. doi: [10.1021/acs.jproteome.5c00557](https://doi.org/10.1021/acs.jproteome.5c00557).
- **QC4DIANN:**
> Moschem JDC, de Barros BCSC, Serrano SMT, Chaves AFA. Decoding the Impact of Isolation Window Selection and QuantUMS Filtering in DIA-NN for DIA Quantification of Peptides and Proteins. *J Proteome Res.* 2025;24(8):3860-3873. doi: [10.1021/acs.jproteome.5c00009](https://doi.org/10.1021/acs.jproteome.5c00009).