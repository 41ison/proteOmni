# proteOmni tests

This folder holds three `testthat` suites, for the PwrQuant module
(`app/modules/mod_PwrQuant.r`), the FASTA / iBAQ helpers
(`app/modules/utils_fasta.r`) and the Sage module's iBAQ calculation
(`app/modules/mod_Sage.r`), and three standalone scripts that are run by
hand.

| Path | What it is | Runtime |
|------|------------|---------|
| `PwrQuant/` | `testthat` unit tests for `mod_PwrQuant.r` | ~25 s |
| `utils_fasta/` | `testthat` unit tests for `utils_fasta.r` | ~1 s |
| `Sage/` | `testthat` unit tests for `compute_ibaq_sage()` in `mod_Sage.r` | ~1 s |
| `imputation_strategies_sim.R` | Simulation behind the FDR / power table in the module's *Imputation* documentation | 15-25 min |
| `test_mdd.R` | Smoke test of the minimum-detectable-difference machinery on the training matrix | ~20 s |
| `pwr_diag.R` | Temporary diagnostic of the Shiny reactive graph; not a test | ~10 s |

## Unit tests (`PwrQuant/`)

### Requirements

- R with the app's packages installed (run `source("bootstrap.R")` from `app/`
  once if you have not already).
- `testthat` (3rd edition) and `withr`.
- The training matrix `app/training/mtx_maxLFQ.tsv`. It is read in place and
  never modified.

The tests source the module to *define* its functions. No Shiny app is started
and no network access is needed.

### Running the suite

From a terminal, in the folder that contains `app/`:

```sh
Rscript app/tests/PwrQuant/run_tests.R
```

or from inside `app/`:

```sh
Rscript tests/PwrQuant/run_tests.R
```

From an R session whose working directory is either of those folders:

```r
source("app/tests/PwrQuant/run_tests.R")   # or "tests/PwrQuant/run_tests.R"
```

To run a single file, pass a `filter` matching the part of the file name after
`test-`:

```r
testthat::test_dir("app/tests/PwrQuant", filter = "imputation")
```

`helper-pwrquant.R` is sourced automatically before every test file. It loads
the module, reads the training matrix (3 conditions `m0`, `m1`, `m2` x 3
replicates) and builds the shared fixtures (`pq_log2`, `pq_log2_67`,
`pq_fit`, ...), so test files should not be `source()`d on their own.

### Expected result

A passing run ends like this:

```
✔ | F W  S  OK | Context
✔ |        149 | contrasts
✔ |          8 | cv
✔ |         12 | id-helpers
✔ |         59 | imputation
✔ |         21 | ingestion
✔ |         19 | limma-pipeline
✔ |         16 | missingness
✔ |        376 | plot-helpers
✔ |         15 | power-helpers
✔ |         24 | server [22.3s]
✔ |         22 | simulation
✔ |         18 | tmm

══ Results ══════════════════════════════════════════════════════════════
Duration: 24.6 s

[ FAIL 0 | WARN 0 | SKIP 0 | PASS 739 ]
```

How to read it:

- **`FAIL 0` is the requirement.** `run_tests.R` stops with an error (non-zero
  exit status under `Rscript`) if any expectation fails.
- **`WARN 0`** is expected. A warning usually means an imputation fallback was
  triggered or a package changed behaviour; investigate it.
- **`SKIP`** is above zero only when an optional package is missing. The
  skipped tests are listed with the package name, and the `PASS` count drops
  accordingly:

  | Missing package | Tests skipped |
  |-----------------|---------------|
  | `impute` | KNN imputation and k selection |
  | `missForest` | missForest imputation |
  | `pcaMethods` | Bayesian PCA imputation |
  | `edgeR` | all of `test-tmm.R` and the TMM pipeline test |
  | `shiny`, `shinyjs`, `shinydashboard`, `tidyverse`, `DT`, `data.table` | `test-server.R` |

- **`PASS`** counts expectations, not tests, so it changes whenever tests are
  added. 739 is the count at the time of writing, with every optional package
  installed.
- The duration depends on the machine. `server` is the slowest context because
  it drives the real upload reactives through `shiny::testServer()`.

### What each file covers

| File | Covers |
|------|--------|
| `test-ingestion.R` | Shape, IDs and missingness profile of the training matrix; `log2(x + 1)`; the groupwise valid-value filter |
| `test-contrasts.R` | `build_contrast_matrix()`, `contrast_pair_labels()` and `contrast_group_table()` with free-form condition names (regex metacharacters, operators, hyphens); rejection of unresolvable or unsafe contrasts without evaluating them |
| `test-cv.R` | `compute_cv_mtx()` and its compiled twin |
| `test-id-helpers.R` | UniProt accession parsing, GO key-type detection, organism tables |
| `test-imputation.R` | `groupwise_imputation()` for MinProb, KNN, missForest and BPCA; `select_knn_k()` |
| `test-missingness.R` | The intensity-dependent dropout model and simulated missingness |
| `test-tmm.R` | TMM normalization, its fallbacks and input validation |
| `test-limma-pipeline.R` | The limma fit on the training matrix and `extract_power_stats()` |
| `test-power-helpers.R` | `bh_effective_alpha()` and `solve_ncp()` |
| `test-simulation.R` | `simulate_pipeline_mdd()` (grid and bisection) and the replicate sweep |
| `test-plot-helpers.R` | Heatmap clustering, palettes and STRING colour helpers |
| `test-server.R` | Upload, log2 toggle, sample exclusion, the end-to-end limma run with free-form condition and display names, and the per-run parameter log in `PwrQuant_server()` |

The imputation tests check that each method fills every gap it should, leaves
observed values unchanged, and leaves proteins that are entirely missing in a
group as `NA`. For KNN they check that a protein observed only once in a group
is filled from neighbours of similar abundance rather than with the sample
mean. They also check that MinProb draws come from each sample's own
left tail and do not move when the protein's own values are shifted; that
property holds for MinProb only, since the other three methods use the
protein's observed values by design.

### If a test fails

1. Read the failure block: it names the file, the test and the expectation.
2. Re-run just that file with `filter =` (see above).
3. Tests that draw random numbers set their own seed, so a failure is
   reproducible. If it is not, suspect a package update; compare your package
   versions with `app/proteOmni.lock`.
4. `Could not locate modules/mod_PwrQuant.r` or `Cannot find the PwrQuant test
   folder` means the working directory is wrong; use one of the two folders
   named under *Running the suite*.

## Unit tests (`utils_fasta/`)

These cover the FASTA and iBAQ helpers in `app/modules/utils_fasta.r`, which
are shared by the QC4DIANN, Sage and PSManalyst modules. They need only
`testthat`, `withr`, `dplyr`, `tidyr` and `tibble`: every fixture is built in
`helper-fasta.R` (four small proteins whose tryptic digests can be counted by
hand, written to a temporary FASTA when a test needs a file), so no data file,
Shiny session or network access is involved.

Run them the same way as the PwrQuant suite:

```sh
Rscript app/tests/utils_fasta/run_tests.R
```

```r
source("app/tests/utils_fasta/run_tests.R")   # or "tests/utils_fasta/run_tests.R"
testthat::test_dir("app/tests/utils_fasta", filter = "ibaq")
```

A passing run ends like this:

```
✔ | F W  S  OK | Context
✔ |         17 | coverage
✔ |         44 | digest
✔ |         27 | fasta-parsing
✔ |         46 | ibaq
✔ |        130 | seqprofile

══ Results ══════════════════════════════════════════════════════════════
[ FAIL 0 | WARN 0 | SKIP 0 | PASS 264 ]
```

`FAIL 0`, `WARN 0` and `SKIP 0` are all expected. The `seqprofile` tests need
the `Peptides` package (pinned in `proteOmni.lock`), plus `shiny`, `plotly`,
`DT` and `ggplot2` for the module-server and plot tests; they are skipped when
a package is missing.

| File | Covers |
|------|--------|
| `test-fasta-parsing.R` | `read_fasta_custom()`, including the regression test for non-UniProt headers whose first word contains "sp" or "tr"; header helpers (`fasta_accessions()`, `fasta_header_ids()`, `fasta_header_genes()`, `fasta_header_tokens()`); whole-token lookups (`fasta_index_exact()`, `protein_ids_contain()`); `fasta_seq_normalize()`, `fasta_sequences_by_header_id()`, `match_fasta_sequences()` and `ibaq_fasta_sequences()` |
| `test-digest.R` | `count_theoretical_peptides()` and its agreement with `in_silico_digest()`; `build_theoretical_peptide_counts()`; `protease_regex_from_rule()` for C- and N-terminal cleavage; `classify_peptides()`; `compute_protein_abundance()` (sum and top3) |
| `test-coverage.R` | `peptide_starts()` and `protein_coverage()` |
| `test-ibaq.R` | `filter_diann_report()` (each filter on its own, inclusive cutoffs, absent columns) and `compute_ibaq()` (wide and long output, log2, unmatched protein groups, digest settings, error messages, and the regression test for `id_col = "Genes"` / `"Protein.Names"`) |
| `test-seqprofile.R` | Sequence profiling: `clean_protein_sequence()` (raw vs FASTA, lower case / digits / spaces, several entries, non-canonical and illegal characters, empty input), `seqprofile_windows()`, `extract_profile_regions()`, `run_seq_profile()` for every registered algorithm (data contract, agreement with direct `Peptides` calls and with the in-repo `GRAVY()`, masking of windows with non-canonical residues, classed errors for short sequences and bad parameters), `seqprofile_global_descriptors()`, `plot_seq_profile()`, and `SeqProfile_server()` in `mod_SeqProfile.r` (computes only on **Predict**) |

Not covered: the peptide property and plotting helpers in the same file
(`GRAVY()`, `calculate_pI()`, `calculate_MW()`, `plot_2d_gel()`),
`extract_diann_matrix()` and the `Fasta_server()` Shiny module.

## Unit tests (`Sage/`)

These cover `compute_ibaq_sage()` in `app/modules/mod_Sage.r`. The helper
sources `utils_fasta.r` and then `mod_Sage.r`, as the app does, so every
package the Sage module attaches at the top of the file must be installed
(`shiny`, `shinydashboard`, `shinyjs`, `arrow`, `ggpointdensity`, `ggtext`,
`DT`, ...); a missing one stops the suite with an error rather than a skip.
Fixtures are built in `helper-sage.R`: a five-protein FASTA keyed by header ID,
as Sage reports proteins, and a small filtered LFQ table. No Sage output
folder, Shiny session or network access is involved.

```sh
Rscript app/tests/Sage/run_tests.R
```

```r
source("app/tests/Sage/run_tests.R")   # or "tests/Sage/run_tests.R"
```

A passing run ends like this:

```
✔ | F W  S  OK | Context
✔ |         50 | ibaq-sage

══ Results ══════════════════════════════════════════════════════════════
[ FAIL 0 | WARN 0 | SKIP 0 | PASS 50 ]
```

| File | Covers |
|------|--------|
| `test-ibaq-sage.R` | iBAQ values per protein and file; accession, entry name and gene annotation; row order; proteins absent from the FASTA or without an observable peptide; the accession fallback when headers differ from the search database; contaminant exclusion and its prefixes; digest settings; the empty-input error; a run on the output of `filter_sage_lfq()` |

Not covered: the rest of `mod_Sage.r` (file discovery, the protein roll-up,
`read_sage_digest_params()`, the plots and the Shiny server).

## Standalone scripts

These are not picked up by `testthat`. Run them from the folder that contains
`app/`, because they use paths relative to it. All three read the training
matrix `app/training/mtx_maxLFQ.tsv` (3 conditions `m0`, `m1`, `m2` x 3
replicates) in place and source `app/modules/mod_PwrQuant.r`; no other data
file is needed.

### `imputation_strategies_sim.R`

Reproduces the realised FDR / power table shown in the module documentation
and in `app/README.md`. The full run takes 15-25 min (about 9 min without
missForest).

```sh
Rscript app/tests/imputation_strategies_sim.R
# quick look: fewer simulations, no missForest
PQ_SIM_NSIM=3 PQ_SIM_SKIP_MISSFOREST=1 Rscript app/tests/imputation_strategies_sim.R
```

Expect a progress log followed by a table with one row per dropout level and
true log2FC, each cell being `realised FDR / power (%)`. With the full
settings the numbers should agree with the documented table to within
simulation noise; with `PQ_SIM_NSIM=3` they are noticeably noisier.

It needs `limma`, `dplyr`, `tidyr` and `data.table`, plus the package behind
each imputation arm (`impute` for KNN, `pcaMethods` for BPCA, `missForest`).

### `test_mdd.R`

```sh
Rscript app/tests/test_mdd.R
```

Loads only the statistical helpers of the module, not the Shiny server. It
keeps proteins with at least 2 of 3 valid values in every condition (4055 of
6746, 1.6% missing), imputes with MinProb, applies cyclic loess normalization
and fits the contrasts `m1 − m0` and `m2 − m0`. It then prints, in order:

1. the conditional MDD per contrast;
2. the prospective MDD from a grid search on complete data;
3. the same from bisection, with the largest difference from the grid;
4. a replicate sweep (3 to 15 per group) for `m1 − m0`, on complete data and
   with simulated missingness followed by MinProb imputation, the penalty
   from missingness, and a comparison with 1/sqrt(n) scaling.

The script stops with an error if a helper is missing from the module, if the
missingness model cannot be fitted, if the MDD increases with the number of
replicates, or if a seeded rerun of the sweep gives different values. A clean
run ends with `Monotonicity in n: OK` and `Seeded reproducibility: OK`.
The agreement between grid and bisection is printed but not asserted.

It needs `limma`, `dplyr`, `purrr`, `tibble`, `ggplot2`, `stringr`, `tidyr`
and `data.table`.

### `pwr_diag.R`

```sh
Rscript app/tests/pwr_diag.R
```

Sources the whole module, so it needs the packages the module attaches
(`shiny`, `shinydashboard`, `shinyjs`, `tidyverse`, `DT`, `ComplexHeatmap`,
...; the full list is at the top of the script). The `reactlog` package is not
needed: the log is read with `shiny::reactlog()`.

It uploads the training matrix into `PwrQuant_server()` through
`shiny::testServer()`, assigns the conditions, applies an empty sample
exclusion, then excludes two `m2` replicates so that `m2` becomes a
single-replicate condition. It prints the upload time, the matrix dimensions,
the replicate counts, the thin and emptied conditions, and the 25 reactives
that ran most often. Each of those is shown with its label and `reactId`;
unnamed reactives are labelled with the start of their code, and Shiny's
internal output observers (`observe({ prom <- NULL ...`) cannot be traced back
to an output from the log.

There are no assertions, so there is no pass or fail; the `geom_line()` and
`Too few points to calculate an ellipse` messages are expected with three
replicates per condition.
