# deZombotti_sleep_estimate.R

A standardized analysis pipeline for evaluating the performance of sleep-tracking devices and algorithms. This script applies the de Zambotti analytical framework (Menghini et al., 2021) to epoch-by-epoch actigraphy data from both GGIR-based and machine learning (ML)-based sleep detection methods.

---

## Reference

Menghini, L., Cellini, N., Goldstone, A., Baker, F. C., & de Zambotti, M. (2021). A standardized framework for testing the performance of sleep-tracking technology: step-by-step guidelines and open-source code. *Sleep, 44*(2), zsaa170.
[https://doi.org/10.1093/sleep/zsaa170](https://doi.org/10.1093/sleep/zsaa170)

The framework functions (in the `Functions/` folder) are from the original repository:
[https://github.com/SRI-human-sleep/sleep-trackers-performance](https://github.com/SRI-human-sleep/sleep-trackers-performance)

---

## What the Script Does

The script takes epoch-by-epoch sleep/wake classifications from a device (or model) and a reference standard, then computes a comprehensive set of performance metrics across subjects.

It runs two parallel analysis pipelines:

**GGIR Algorithm Analysis** evaluates three GGIR-based sleep detection algorithms — Sadeh, van Hees, and Cole-Kripke — comparing their output against a machine learning reference. For each algorithm, the script matches subjects' GGIR output files with their ML reference prediction files, merges the epoch-level labels, and runs the full de Zambotti pipeline.

**Machine Learning Model Analysis** evaluates seven ML sleep classification models — LSTM, Random Forest, Neural Network (MLP), k-Nearest Neighbor, Extreme Gradient Boosting, Logistic Regression, and Decision Tree — using their cross-validated epoch-by-epoch predictions. The same analytical pipeline is applied to each model.

For each algorithm or model, the pipeline computes:

- **Sleep metrics per subject** (TIB, TST, SE, SOL, WASO) for both the device and the reference
- **Individual discrepancies** — the difference between device and reference for each metric per subject
- **Group-level Bland-Altman statistics** — bias, standard deviation, and limits of agreement (LOA) with confidence intervals, with automatic testing for proportional bias and heteroscedasticity
- **Bland-Altman plots** — one SVG figure per sleep metric (TST, SE, SOL, WASO)
- **Epoch-by-epoch error matrix** — summed confusion matrix across all subjects
- **Individual EBE metrics** — accuracy, sensitivity, specificity, PPV, NPV, F1 score, and balanced accuracy per subject
- **Group-level EBE summary** — averaged EBE metrics with confidence intervals, plus advanced metrics (Cohen's kappa, PABAK, bias index, prevalence index)

Results for all algorithms or models are aggregated and saved as combined CSV files. A publication-ready version of the group discrepancy table (with full algorithm names and merged LOA/CI columns) is also exported.

---

## Stage Coding Convention

This script uses the de Zambotti two-state (sleep/wake) convention:

| State | Code |
|-------|------|
| Wake  | 10   |
| Sleep | 5    |

GGIR text labels (`"Sleep"`, `"Wake"`) and ML numeric labels (`1` = Sleep, `0` = Wake) are automatically recoded to this convention during processing.

---

## Input Files

All input files must be in **CSV format** (`.csv`).

### GGIR Epoch-by-Epoch Files

Located in: `<GGIR_PATH>/<algorithm_name>/Epoch_b_Epoch/`

Naming convention: `<prefix>_label_<subjectID>.csv`

Required columns:

| Column       | Type   | Description                                      |
|--------------|--------|--------------------------------------------------|
| `GGIR_label` | string | GGIR sleep/wake classification (`Sleep` or `Wake`) |
| `label`      | string | Original reference label (`Sleep` or `Wake`)     |
| `Epoch`      | string | Epoch timestamp (`YYYY-MM-DD HH:MM:SS`)          |

Example file: [`examples/example_ggir_label_S001.csv`](examples/example_ggir_label_S001.csv)

### ML Reference Prediction Files (for GGIR section)

Located in: `<ML_REF_PATH>/` (e.g., `ML_outputs_rf/subject_fold/`)

Naming convention: `pred_<subjectID>_base<suffix>.csv`

Required columns:

| Column    | Type    | Description                          |
|-----------|---------|--------------------------------------|
| `label_n` | integer | True reference label (1=Sleep, 0=Wake) |

Example file: [`examples/example_pred_S001_base_fold1.csv`](examples/example_pred_S001_base_fold1.csv)

### ML Model Prediction Files (for ML section)

Located in: `<DATA_PATH>/ML_outputs_<model_name>/subject_fold/`

Naming convention: `pred_<subjectID>_<suffix>.csv`

Required columns:

| Column        | Type    | Description                             |
|---------------|---------|-----------------------------------------|
| `.pred_class` | integer | Model-predicted class (1=Sleep, 0=Wake) |
| `label_n`     | integer | True reference label (1=Sleep, 0=Wake)  |

Example file: [`examples/example_pred_S001_rf.csv`](examples/example_pred_S001_rf.csv)

---

## Output Files

All output CSV files (`.csv`) and SVG plots (`.svg`) are saved to `<OUTPUT_PATH>/GGIR/` and `<OUTPUT_PATH>/ML/` respectively.

### Per-Algorithm Outputs (saved under `results/<algorithm_name>/`)

| File                         | Description                                                   |
|------------------------------|---------------------------------------------------------------|
| `sleep_metrics.csv`          | Per-subject TIB, TST, SE, SOL, WASO for device and reference |
| `ind_discr.csv`              | Per-subject device minus reference discrepancies              |
| `group_discr.csv`            | Group Bland-Altman statistics (bias, SD, LOAs, CIs)           |
| `error_metrics.csv`          | Epoch-by-epoch confusion matrix (summed)                      |
| `group_ebe_sleep.csv`        | Group-level averaged EBE metrics with advanced statistics     |
| `ind_ebe_sleep.csv`          | Per-subject EBE metrics for sleep and wake combined           |
| `ebe_stat_summary.csv`       | Summary of EBE metrics (mean ± SD across subjects)            |
| `ebe_stat.csv`               | Per-subject EBE metrics (sleep stage)                         |
| `tst_ba_plot.svg`            | Bland-Altman plot: Total Sleep Time                           |
| `se_ba_plot.svg`             | Bland-Altman plot: Sleep Efficiency                           |
| `sol_ba_plot.svg`            | Bland-Altman plot: Sleep Onset Latency                        |
| `waso_ba_plot.svg`           | Bland-Altman plot: Wake After Sleep Onset                     |
| `logtf_boot_tst_ba_plot.svg` | Bland-Altman plot: TST (bootstrap CI variant)                 |

### Aggregated Outputs (saved directly under `GGIR/` or `ML/`)

| File                                   | Description                                              |
|----------------------------------------|----------------------------------------------------------|
| `final_slp_estimate_all.csv`           | Sleep metrics combined across all algorithms/models      |
| `final_ind_disc_dt_all.csv`            | Individual discrepancies combined across all algorithms  |
| `final_group_disc_all.csv`             | Group discrepancy statistics combined                    |
| `final_group_disc_all_publication.csv` | Publication-ready version with full algorithm names      |
| `final_ref_time_all.csv`              | Reference sleep onset timestamps *(GGIR only)*           |
| `final_ggir_time_all.csv`             | GGIR-detected sleep onset timestamps *(GGIR only)*       |
| `ggir_summary_sleeponset.csv`         | Sleep onset hour summary by algorithm *(GGIR only)*      |

---

## Setup and Usage

### Prerequisites

R (≥ 4.0) must be installed. Required packages:

```r
install.packages(c(
  "data.table", "tidyverse", "caret", "zeallot", "moments",
  "reshape2", "ggplot2", "BlandAltmanLeh", "ggExtra", "epiR", "DescTools"
))
```

### Step 1: Clone or download this repository

```bash
git clone https://github.com/SRI-human-sleep/sleep-trackers-performance.git
```

Place `deZombotti_sleep_estimate.R` in the repository root alongside the `Functions/` folder.

### Step 2: Prepare your data

Organize your data directory as follows:

```
<DATA_PATH>/
├── GGIR_result/
│   ├── GGIR_334_Sadeh/
│   │   └── Epoch_b_Epoch/
│   │       ├── prefix_label_S001.csv
│   │       └── ...
│   ├── GGIR_334_vH/
│   └── GGIR_334_CK/
├── ML_outputs_rf/
│   └── subject_fold/
│       ├── pred_S001_base_fold1.csv
│       └── ...
├── ML_outputs_lstm/
├── ML_outputs_nn_mlp/
├── ML_outputs_knn/
├── ML_outputs_xgb/
├── ML_outputs_log_reg/
└── ML_outputs_tree/
```

### Step 3: Set your paths

Open `deZombotti_sleep_estimate.R` and update **Section 2** with your local paths:

```r
DATA_PATH      <- "/path/to/your/data"
GGIR_PATH      <- file.path(DATA_PATH, "GGIR_result")
ML_REF_PATH    <- file.path(DATA_PATH, "ML_outputs_rf", "subject_fold")
OUTPUT_PATH    <- file.path(DATA_PATH, "deZombotti_sleep_estimate")
FUNCTIONS_PATH <- file.path(dirname(rstudioapi::getActiveDocumentContext()$path), "Functions")
```

### Step 4: Run the script

Open `deZombotti_sleep_estimate.R` in RStudio and click **Source**, or run from the R console:

```r
source("deZombotti_sleep_estimate.R")
```

---

## Example Files

The `examples/` folder contains three small synthetic CSV files (200 epochs each, subject ID `S001`) that illustrate the expected input format:

| File                               | Represents                             |
|------------------------------------|----------------------------------------|
| `example_ggir_label_S001.csv`      | GGIR epoch-by-epoch output             |
| `example_pred_S001_base_fold1.csv` | ML reference prediction (GGIR section) |
| `example_pred_S001_rf.csv`         | ML model prediction (ML section)       |

These files are **not real participant data** — all values are randomly generated for illustration purposes only.

---

## Running Unit Tests

Unit tests are provided in `test_deZombotti_pipeline.R` and use the `testthat` framework. No real data is needed — all tests use synthetic data generated internally.

```r
install.packages("testthat")
source("test_deZombotti_pipeline.R")
```

The tests cover: `ebe2sleep`, `indDiscr`, `errorMatrix`, `indEBE`, and the label recoding logic. Each test checks a specific property (row counts, column names, value bounds, arithmetic correctness, and edge cases).

Expected output:
```
[ FAIL 0 | WARN 0 | SKIP 0 | PASS 20 ]
```

---

## Security Notes

This script has been reviewed for public sharing. The following types of information have been removed:

- Institutional server and network share paths
- Personal usernames embedded in file paths
- Internal study identifiers and version-dated folder names
- Commented-out paths referencing unpublished manuscript titles

If you adapt this script for your own study, ensure that no personal or institutional identifiers remain in your copy before pushing to a public repository.

---

## License

This repository follows the license of the upstream [sleep-trackers-performance](https://github.com/SRI-human-sleep/sleep-trackers-performance) repository. Please refer to the `LICENSE` file for details.
