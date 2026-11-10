# proteOmni tests

This folder holds the unit tests for the PwrQuant module
(`app/modules/mod_PwrQuant.r`) and three standalone scripts that are run by
hand.

| Path | What it is | Runtime |
|------|------------|---------|
| `PwrQuant/` | `testthat` unit tests (the suite described below) | ~10 s |
| `imputation_strategies_sim.R` | Simulation behind the FDR / power table in the module's *Imputation* documentation | 15-25 min |
| `test_mdd.R` | Smoke test of the minimum-detectable-difference machinery on a real matrix | minutes |
| `pwr_diag.R` | Temporary diagnostic of the Shiny reactive graph; not a test | seconds |

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
✔ |          8 | cv
✔ |         12 | id-helpers
✔ |         59 | imputation
✔ |         21 | ingestion
✔ |         19 | limma-pipeline
✔ |         16 | missingness
✔ |        376 | plot-helpers
✔ |         15 | power-helpers
✔ |          6 | server [6.3s]
✔ |         22 | simulation
✔ |         18 | tmm

══ Results ══════════════════════════════════════════════════════════════
Duration: 8.6 s

[ FAIL 0 | WARN 0 | SKIP 0 | PASS 572 ]
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
  added. 572 is the count at the time of writing, with every optional package
  installed.
- The duration depends on the machine. `server` is the slowest context because
  it drives the real upload reactives through `shiny::testServer()`.

### What each file covers

| File | Covers |
|------|--------|
| `test-ingestion.R` | Shape, IDs and missingness profile of the training matrix; `log2(x + 1)`; the groupwise valid-value filter |
| `test-cv.R` | `compute_cv_mtx()` and its compiled twin |
| `test-id-helpers.R` | UniProt accession parsing, GO key-type detection, organism tables |
| `test-imputation.R` | `groupwise_imputation()` for MinProb, KNN, missForest and BPCA; `select_knn_k()` |
| `test-missingness.R` | The intensity-dependent dropout model and simulated missingness |
| `test-tmm.R` | TMM normalization, its fallbacks and input validation |
| `test-limma-pipeline.R` | The limma fit on the training matrix and `extract_power_stats()` |
| `test-power-helpers.R` | `bh_effective_alpha()` and `solve_ncp()` |
| `test-simulation.R` | `simulate_pipeline_mdd()` (grid and bisection) and the replicate sweep |
| `test-plot-helpers.R` | Heatmap clustering, palettes and STRING colour helpers |
| `test-server.R` | Upload, log2 toggle and sample exclusion in `PwrQuant_server()` |

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

## Standalone scripts

These are not picked up by `testthat` and take much longer. Run them from the
folder that contains `app/`, because they use paths relative to it.

**`imputation_strategies_sim.R`** reproduces the realised FDR / power table
shown in the module documentation and in `app/README.md`.

```sh
Rscript app/tests/imputation_strategies_sim.R
# quick look: fewer simulations, no missForest
PQ_SIM_NSIM=3 PQ_SIM_SKIP_MISSFOREST=1 Rscript app/tests/imputation_strategies_sim.R
```

Expect a progress log followed by a table with one row per dropout level and
true log2FC, each cell being `realised FDR / power (%)`. With the full
settings the numbers should agree with the documented table to within
simulation noise; with `PQ_SIM_NSIM=3` they are noticeably noisier.

**`test_mdd.R`** needs `matrix_human.tsv` in the working directory. It prints
its intermediate results and stops with an error if a check fails; a clean run
ends with `Monotonicity in n: OK` and `Seeded reproducibility: OK`.

**`pwr_diag.R`** needs `mtx_prot_abundance.tsv` in the working directory and
the `reactlog` package. It prints timings and reactive execution counts for an
upload and makes no assertions, so there is no pass or fail.
