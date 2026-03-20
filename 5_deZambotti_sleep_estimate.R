# ============================================================
# 5_deZombotti_sleep_estimate.R
#
# Title:       Sleep Tracking Performance Analysis Pipeline
# Description: Reorganizes epoch-by-epoch actigraphy data to apply
#              the de Zambotti analytical framework for evaluating
#              sleep tracking performance. Computes sleep metrics,
#              discrepancy analyses, Bland-Altman statistics, and
#              epoch-by-epoch (EBE) metrics for both GGIR-based and
#              machine learning (ML)-based sleep detection methods.
#              The script is designed to compare both GGIR and ML.
#              The script serve as an example how to use the sleep
#              tracker performance functions. 
#
# Reference:
#   Menghini, L., Cellini, N., Goldstone, A., Baker, F. C., &
#   de Zambotti, M. (2021). A standardized framework for testing the
#   performance of sleep-tracking technology: step-by-step guidelines
#   and open-source code. Sleep, 44(2), zsaa170.
#   https://doi.org/10.1093/sleep/zsaa170
#
# Input files (CSV format):
#   - GGIR epoch-by-epoch files:
#       Named as: <prefix>_label_<subjectID>.csv
#       Required columns: GGIR_label (Sleep/Wake), label (Sleep/Wake),
#                         Epoch (timestamp)
#   - ML prediction files:
#       Named as: pred_<subjectID>_<suffix>.csv
#       Required columns: .pred_class (0=Wake, 1=Sleep),
#                         label_n (0=Wake, 1=Sleep)
#   - ML reference files (for GGIR section):
#       Named as: pred_<subjectID>_<suffix>.csv
#       Required columns: label_n (0=Wake, 1=Sleep)
#
# Output files (CSV and SVG format):
#   - sleep_metrics.csv         : Per-subject sleep metrics (TIB, TST, SE, SOL, WASO)
#   - ind_discr.csv             : Individual-level device-reference discrepancies
#   - group_discr.csv           : Group-level Bland-Altman discrepancy statistics
#   - error_metrics.csv         : Epoch-by-epoch confusion matrix
#   - group_ebe_sleep.csv       : Group-level EBE performance metrics
#   - ind_ebe_sleep.csv         : Individual EBE metrics (sleep + wake combined)
#   - ebe_stat_summary.csv      : Summary of EBE metrics across subjects
#   - ebe_stat.csv              : Per-subject EBE metrics
#   - *_ba_plot.svg             : Bland-Altman plots per sleep metric
#   - final_*.csv               : Aggregated results across all algorithms/models
# ============================================================


# ============================================================
# SECTION 1: REQUIRED PACKAGES
# ============================================================
# Install any missing packages before loading
required_packages <- c(
  "data.table", "tidyverse", "caret", "zeallot", "moments",
  "reshape2", "ggplot2", "BlandAltmanLeh", "ggExtra", "epiR", "DescTools"
)
new_packages <- required_packages[!(required_packages %in% installed.packages()[, "Package"])]
if (length(new_packages) > 0) install.packages(new_packages)

library(data.table)
library(tidyverse)
library(caret)
library(zeallot)
library(moments)
library(reshape2)
library(ggplot2)
library(BlandAltmanLeh)
library(ggExtra)
library(epiR)
library(DescTools)

# Clear the workspace before running
rm(list = ls())


# ============================================================
# SECTION 2: USER CONFIGURATION
# ============================================================

# Root path to your data directory (parent folder containing all input data)
# Example: DATA_PATH <- "C:/Users/YourName/Documents/SleepStudyData"
DATA_PATH <- "/path/to/your/data"

# Path to the GGIR results folder (contains algorithm subfolders with epoch-by-epoch CSVs)
# Expected structure: <GGIR_PATH>/<algorithm_name>/Epoch_b_Epoch/<subjectID_files>.csv
GGIR_PATH <- file.path(DATA_PATH, "GGIR_result")

# Path to the ML reference prediction folder (subject-fold cross-validated predictions)
# The files this path pointed to serve as the subject numbers that you want to include
# The goal is to make sure you only include subjects you want to analyze
# ml_ids will be extracted from these files.
# Expected file names: pred_<subjectID>_<suffix>.csv
ML_REF_PATH <- file.path(DATA_PATH, "ML_outputs_rf", "subject_fold")

# Root output directory where all de Zambotti analysis results will be saved
# Subdirectories (GGIR/ and ML/) will be created automatically
OUTPUT_PATH <- file.path(DATA_PATH, "deZambotti_sleep_estimate")

# Path to the Functions folder containing the de Zambotti framework .R files
# If running from the cloned repository root, this relative path should work:
FUNCTIONS_PATH <- file.path("PATH/to/Functions/folder", "Functions")
# For example:
# FUNCTIONS_PATH <- "/path/to/sleep-trackers-performance/Functions"


# ============================================================
# SECTION 3: CREATE OUTPUT DIRECTORIES
# ============================================================
# Create output subdirectories for GGIR and ML results.
# showWarnings = FALSE suppresses messages if directories already exist.
GGIR_OUTPUT_PATH <- file.path(OUTPUT_PATH, "GGIR")
ML_OUTPUT_PATH   <- file.path(OUTPUT_PATH, "ML")

dir.create(GGIR_OUTPUT_PATH, recursive = TRUE, showWarnings = FALSE)
dir.create(ML_OUTPUT_PATH,   recursive = TRUE, showWarnings = FALSE)


# ============================================================
# SECTION 4: LOAD DE ZAMBOTTI FRAMEWORK FUNCTIONS
# ============================================================
# Source all .R function files from the Functions folder.
# These implement the standardized analytical pipeline from
# Menghini et al. (2021) 
r_files <- list.files(FUNCTIONS_PATH, pattern = "\\.R$", full.names = TRUE)
if (length(r_files) == 0) {
  stop("No .R function files found at FUNCTIONS_PATH: ", FUNCTIONS_PATH,
       "\nPlease set FUNCTIONS_PATH to the 'Functions' folder of the sleep-trackers-performance repo.")
}
invisible(lapply(r_files, source))
message("Loaded ", length(r_files), " function(s) from: ", FUNCTIONS_PATH)


# ============================================================
# SECTION 5: HELPER FUNCTION — RUN DE ZAMBOTTI PIPELINE
# ============================================================
# run_dezombotti_pipeline() encapsulates the complete analysis workflow
# for a single algorithm or model. It is called for each GGIR algorithm
# and each ML model in the loops below.
#
# Arguments:
#   all_data   : data.table with columns: subject, device, reference
#                device and reference must be numeric: wake = 10, sleep = 5
#   output_dir : character, path to save CSVs and SVG plots for this run
#   label      : character, human-readable label for progress messages
#
# Returns:
#   A named list: sleep_metrics, ind_discr, group_discr
#   (used for cross-algorithm aggregation in Sections 6 and 7)

run_dezombotti_pipeline <- function(all_data, output_dir, label) {

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  # ----------------------------------------------------------
  # Step 5a: Convert epoch-by-epoch data to per-subject sleep metrics
  # ----------------------------------------------------------
  # ebe2sleep() computes:
  #   TIB  = Time in Bed (total recording duration, minutes)
  #   TST  = Total Sleep Time (minutes)
  #   SE   = Sleep Efficiency (% of TIB scored as sleep)
  #   SOL  = Sleep Onset Latency (minutes of wake before first sleep epoch)
  #   WASO = Wake After Sleep Onset (minutes of wake after first sleep epoch)
  # Input CSV format: columns subject, device, reference (numeric: 5/10)
  sleep_metrics <- ebe2sleep(
    data        = all_data,
    idCol       = "subject",
    RefCol      = "reference",
    deviceCol   = "device",
    epochLenght = 30,            # epoch length in seconds
    staging     = FALSE,         # two-state classification (sleep/wake only)
    stages      = c(wake = 10, sleep = 5),
    digits      = 2
  )

  # ----------------------------------------------------------
  # Step 5b: Compute individual-level discrepancies (device minus reference)
  # ----------------------------------------------------------
  # indDiscr() subtracts reference values from device values for each
  # sleep metric (TST_diff, SE_diff, SOL_diff, WASO_diff).
  # Positive values indicate overestimation by the device.
  # Output CSV: ind_discr.csv
  ind_discr <- indDiscr(
    data    = sleep_metrics,
    staging = FALSE,
    digits  = 2,
    doPlot  = FALSE
  )

  # ----------------------------------------------------------
  # Step 5c: Compute group-level discrepancy statistics (Bland-Altman)
  # ----------------------------------------------------------
  # groupDiscr() computes bias, SD, and limits of agreement (LOAs) for
  # each sleep metric across all subjects, following the Bland-Altman (1999)
  # framework. Tests for proportional bias and heteroscedasticity.
  # Output CSV: group_discr.csv
  group_discr <- rbind(
    groupDiscr(data = sleep_metrics, measures = c("TST_device",  "TST_ref"),
               CI.type = "classic", CI.level = .95, digits = 2),
    groupDiscr(data = sleep_metrics, measures = c("SE_device",   "SE_ref"),
               CI.type = "classic", boot.type = "basic", CI.level = .95, digits = 2),
    groupDiscr(data = sleep_metrics, measures = c("SOL_device",  "SOL_ref"),
               CI.type = "classic", CI.level = .95, digits = 2),
    groupDiscr(data = sleep_metrics, measures = c("WASO_device", "WASO_ref"),
               CI.type = "classic", boot.type = "basic", CI.level = .95, digits = 2)
  )

  # ----------------------------------------------------------
  # Step 5d: Generate Bland-Altman plots for each sleep metric
  # ----------------------------------------------------------
  # BAplot() creates a scatter plot of (device - reference) differences
  # against reference values, with bias and LOA lines overlaid.
  # Output SVG files (21 x 10 cm): *_ba_plot.svg
  tst_ba_plot       <- BAplot(data = sleep_metrics, measures = c("TST_device",  "TST_ref"),
                               xaxis = "reference", CI.type = "classic",
                               CI.level = .95, ylim = c(-300, 200))
  se_ba_plot        <- BAplot(data = sleep_metrics, measures = c("SE_device",   "SE_ref"),
                               xaxis = "reference", CI.type = "classic",
                               CI.level = .95, ylim = c(-40, 40))
  sol_ba_plot       <- BAplot(data = sleep_metrics, measures = c("SOL_device",  "SOL_ref"),
                               xaxis = "reference", CI.type = "classic",
                               CI.level = .95, ylim = c(-200, 200))
  waso_ba_plot      <- BAplot(data = sleep_metrics, measures = c("WASO_device", "WASO_ref"),
                               xaxis = "reference", CI.type = "classic",
                               CI.level = .95, ylim = c(-200, 350))
  logtf_tst_ba_plot <- BAplot(data = sleep_metrics, measures = c("TST_device",  "TST_ref"),
                               xaxis = "reference", CI.type = "classic",
                               boot.type = "basic", CI.level = .95, ylim = c(-200, 350))

  # ----------------------------------------------------------
  # Step 5e: Compute epoch-by-epoch error matrix
  # ----------------------------------------------------------
  # errorMatrix() tallies how many epochs the device scored as wake/sleep
  # against how many the reference scored as wake/sleep, producing a
  # confusion matrix summed across all subjects.
  # matrixType = "sum" returns absolute epoch counts.
  # Output CSV: error_metrics.csv
  err_matrix <- errorMatrix(
    data       = all_data,
    idCol      = "subject",
    RefCol     = "reference",
    deviceCol  = "device",
    staging    = FALSE,
    stages     = c(wake = 10, sleep = 5),
    matrixType = "sum",
    CI.type    = "classic",
    boot.type  = "basic",
    CI.level   = .95,
    digits     = 2
  )

  # ----------------------------------------------------------
  # Step 5f: Compute per-subject epoch-by-epoch (EBE) performance metrics
  # ----------------------------------------------------------
  # indEBE() computes confusion-matrix-based metrics per subject:
  #   accuracy, sensitivity, specificity, PPV, NPV, F1, balanced_acc
  # stage = 5 targets the "sleep" class (positive class).
  # Output CSV: ebe_stat.csv
  ebe_stat <- indEBE(
    data       = all_data,
    stage      = 5,
    stageLabel = "sleep",
    doPlot     = FALSE,
    digits     = 2
  )
  setDT(ebe_stat)

  # ----------------------------------------------------------
  # Step 5g: Summarize EBE metrics across subjects (mean +/- SD)
  # ----------------------------------------------------------
  # Collapses per-subject EBE metrics into group-level summary statistics.
  # Output CSV: ebe_stat_summary.csv
  ebe_stat_summary <- ebe_stat[, .(
    acc_mean          = mean(accuracy),         acc_sd          = sd(accuracy),
    balanced_acc_mean = mean(balanced_acc),     balanced_acc_sd = sd(balanced_acc),
    sens_mean         = mean(sensitivity),      sens_sd         = sd(sensitivity),
    spec_mean         = mean(specificity),      spec_sd         = sd(specificity),
    prevalence_mean   = mean(prevalence),       prevalence_sd   = sd(prevalence),
    ppv_mean          = mean(ppv),              ppv_sd          = sd(ppv),
    npv_mean          = mean(npv),              npv_sd          = sd(npv),
    f1_mean           = mean(f1),               f1_sd           = sd(f1)
  )]

  # ----------------------------------------------------------
  # Step 5h: Compute combined sleep/wake individual EBE table
  # ----------------------------------------------------------
  # Runs indEBE() for both sleep (stage=5) and wake (stage=10) and
  # concatenates the results side by side.
  # Output CSV: ind_ebe_sleep.csv
  ebe_combined <- cbind(
    indEBE(data = all_data, stage = 5,  stageLabel = "sleep", doPlot = FALSE, digits = 2),
    indEBE(data = all_data, stage = 10, stageLabel = "wake",  doPlot = FALSE, digits = 2)
  )

  # ----------------------------------------------------------
  # Step 5i: Compute group-level averaged EBE metrics
  # ----------------------------------------------------------
  # groupEBE() averages individual confusion matrices across subjects
  # (metricsType = "avg") and computes advanced metrics including
  # Cohen's kappa, PABAK, PPV, and NPV.
  # Output CSV: group_ebe_sleep.csv
  grp_ebe_sleep <- groupEBE(
    data            = all_data,
    idCol           = "subject",
    stage           = 5,
    stageLabel      = "sleep",
    metricsType     = "avg",
    doROC           = FALSE,
    CI.type         = "classic",
    advancedMetrics = TRUE
  )

  # ---- Save all CSV outputs ----
  # Output format: CSV (.csv)
  fwrite(sleep_metrics,    file.path(output_dir, "sleep_metrics.csv"))
  fwrite(ind_discr,        file.path(output_dir, "ind_discr.csv"))
  fwrite(group_discr,      file.path(output_dir, "group_discr.csv"))
  fwrite(err_matrix,       file.path(output_dir, "error_metrics.csv"))
  fwrite(grp_ebe_sleep,    file.path(output_dir, "group_ebe_sleep.csv"))
  fwrite(ebe_combined,     file.path(output_dir, "ind_ebe_sleep.csv"))
  fwrite(ebe_stat_summary, file.path(output_dir, "ebe_stat_summary.csv"))
  fwrite(ebe_stat,         file.path(output_dir, "ebe_stat.csv"))

  # ---- Save all plot outputs ----
  # Output format: SVG (.svg), 21 x 10 cm
  ggsave(file.path(output_dir, "tst_ba_plot.svg"),
         plot = tst_ba_plot,       width = 21, height = 10, units = "cm")
  ggsave(file.path(output_dir, "se_ba_plot.svg"),
         plot = se_ba_plot,        width = 21, height = 10, units = "cm")
  ggsave(file.path(output_dir, "sol_ba_plot.svg"),
         plot = sol_ba_plot,       width = 21, height = 10, units = "cm")
  ggsave(file.path(output_dir, "logtf_boot_tst_ba_plot.svg"),
         plot = logtf_tst_ba_plot, width = 21, height = 10, units = "cm")
  ggsave(file.path(output_dir, "waso_ba_plot.svg"),
         plot = waso_ba_plot,      width = 21, height = 10, units = "cm")

  message("  [Done] ", label, " -> results saved to: ", output_dir)

  # Return key results for cross-algorithm aggregation in Sections 6-7
  list(
    sleep_metrics = sleep_metrics,
    ind_discr     = ind_discr,
    group_discr   = group_discr
  )
}


# ============================================================
# SECTION 6: GGIR ALGORITHM ANALYSIS
# ============================================================
# Evaluate three GGIR-based sleep/wake detection algorithms:
#   - Sadeh     (GGIR 3.3.4 with Sadeh scoring)
#   - van Hees  (GGIR built-in van Hees algorithm)
#   - CK        (Cole-Kripke algorithm)
#
# For each algorithm:
#   1. Match GGIR epoch-by-epoch output files with ML reference predictions
#   2. Merge device and reference labels into a single data.table
#   3. Recode labels from text (Sleep/Wake) to numeric (5/10)
#   4. Run the full de Zambotti pipeline via run_dezombotti_pipeline()
#   5. Aggregate results across algorithms for combined export

ggir_algorithms <- c("GGIR_334_Sadeh", "GGIR_334_vH", "GGIR_334_CK")

# Initialise aggregation lists
ls_group_disc_ggir   <- list()
ls_slp_estimate_ggir <- list()
ls_ind_disc_dt_ggir  <- list()
ls_ref_time_all      <- list()
ls_ggir_time_all     <- list()

# match_ids is populated in the first GGIR iteration and reused across
# all subsequent iterations (and by the ML section below), because all
# algorithms share the same set of subjects.
match_ids <- NULL

for (algo_idx in seq_along(ggir_algorithms)) {
  algo_name    <- ggir_algorithms[algo_idx]
  # This is where you save the ebe files from GGIR extraction.
  ggir_label_p <- file.path(GGIR_PATH, algo_name, "Epoch_b_Epoch")

  message("Processing GGIR algorithm: ", algo_name)

  # ----------------------------------------------------------
  # Step 6a: Identify subjects present in both GGIR and ML reference outputs
  # ----------------------------------------------------------
  # GGIR files:   CSV, named as <prefix>_label_<subjectID>.csv
  # ML ref files: CSV, named as pred_<subjectID>_base<suffix>.csv
  ls_ggir_files   <- list.files(ggir_label_p, pattern = "\\.csv$")
  ls_ml_ref_files <- list.files(ML_REF_PATH)

  ggir_ids  <- word(word(ls_ggir_files, 2, sep = "_"), 1, sep = ".csv")
  ml_ids    <- unique(word(ls_ml_ref_files, 3, sep = "_"))
  match_ids <- intersect(ggir_ids, ml_ids)

  message("  Matched subjects: ", length(match_ids))

  # ----------------------------------------------------------
  # Step 6b: Load and merge epoch-by-epoch data for each subject
  # ----------------------------------------------------------
  ls_all_data  <- list()
  ls_ref_time  <- list()
  ls_ggir_time <- list()

  for (subj_idx in seq_along(match_ids)) {
    ID <- match_ids[subj_idx]

    # Select matching files for this subject
    ggir_file <- ls_ggir_files[grep(paste0("label_", ID, "\\.csv"), ls_ggir_files)]
    ml_file   <- ls_ml_ref_files[grep(paste0("pred_", ID, "_base"), ls_ml_ref_files)]

    # Load GGIR epoch-by-epoch labels (CSV)
    # Required columns: GGIR_label (Sleep/Wake), label (Sleep/Wake), Epoch (timestamp)
    dt <- fread(file.path(ggir_label_p, ggir_file))

    # Load ML reference predictions (CSV) and convert label to Sleep/Wake string
    # Required columns: label_n (1 = Sleep, 0 = Wake)
    dt_ref <- fread(file.path(ML_REF_PATH, ml_file))
    dt_ref[, ml_ref := ifelse(label_n == 1, "Sleep", "Wake")]

    # Record the epoch timestamp of the first sleep epoch for sleep-onset analysis
    first_ggir_sleep <- which(dt$GGIR_label %in% "Sleep")[1]
    first_ref_sleep  <- which(dt$label %in% "Sleep")[1]

    dt_ggir_time <- data.table(
      time_ggir = dt[first_ggir_sleep, Epoch],
      subject   = word(ls_ggir_files[subj_idx], 2, sep = "_")
    )
    dt_ref_time <- data.table(
      time_ref = dt[first_ref_sleep, Epoch],
      subject  = word(ls_ggir_files[subj_idx], 2, sep = "_")
    )

    # Merge GGIR output and ML reference by epoch index
    dt[, index := seq_len(.N)]
    dt_ref[, index := seq_len(.N)]
    dt[, ml_ref := dt_ref$ml_ref]

    # Select the columns needed for the pipeline:
    # device       = GGIR algorithm label (Sleep/Wake)
    # reference_or = original PSG/reference label from GGIR file
    # reference    = ML-derived reference label (used as the analysis reference)
    dt_new <- dt[, .(
      subject      = ID,
      device       = GGIR_label,
      reference_or = label,
      reference    = ml_ref
    )]

    ls_all_data[[subj_idx]]  <- dt_new
    ls_ref_time[[subj_idx]]  <- dt_ref_time
    ls_ggir_time[[subj_idx]] <- dt_ggir_time
    rm(dt, dt_new)
  }

  # ----------------------------------------------------------
  # Step 6c: Combine all subjects into a single data.table
  # ----------------------------------------------------------
  all_data      <- setDT(rbindlist(ls_all_data))
  all_ref_time  <- setDT(rbindlist(ls_ref_time))
  all_ggir_time <- setDT(rbindlist(ls_ggir_time))

  # ----------------------------------------------------------
  # Step 6d: Recode text labels to numeric stage codes
  # ----------------------------------------------------------
  # Per the de Zambotti framework: wake = 10, sleep = 5.
  # Empty reference strings (NA epochs) are treated as Wake.
  all_data[reference %in% "", reference := "Wake"]
  all_data[device    == "Wake",  device    := 10]
  all_data[device    == "Sleep", device    := 5]
  all_data[reference == "Wake",  reference := 10]
  all_data[reference == "Sleep", reference := 5]

  # ----------------------------------------------------------
  # Step 6e: Run the de Zambotti pipeline for this GGIR algorithm
  # ----------------------------------------------------------
  results_dir  <- file.path(GGIR_OUTPUT_PATH, "results", algo_name)
  pipeline_out <- run_dezombotti_pipeline(all_data, results_dir, algo_name)

  # ----------------------------------------------------------
  # Step 6f: Tag results with algorithm name and store for aggregation
  # ----------------------------------------------------------
  setDT(pipeline_out$group_discr)[,   analysis_type := algo_name]
  setDT(pipeline_out$sleep_metrics)[, analysis_type := algo_name]
  setDT(pipeline_out$ind_discr)[,     analysis_type := algo_name]
  all_ref_time[,  analysis_type := algo_name]
  all_ggir_time[, analysis_type := algo_name]

  ls_group_disc_ggir[[algo_idx]]   <- pipeline_out$group_discr
  ls_slp_estimate_ggir[[algo_idx]] <- pipeline_out$sleep_metrics
  ls_ind_disc_dt_ggir[[algo_idx]]  <- pipeline_out$ind_discr
  ls_ref_time_all[[algo_idx]]      <- all_ref_time
  ls_ggir_time_all[[algo_idx]]     <- all_ggir_time
}

# ----------------------------------------------------------
# Step 6g: Aggregate GGIR results across all algorithms
# ----------------------------------------------------------
# Combine per-algorithm results into single data.tables
# and save as aggregated CSV files (CSV format).
group_disc_ggir_all   <- rbindlist(ls_group_disc_ggir)
slp_estimate_ggir_all <- rbindlist(ls_slp_estimate_ggir)
ind_disc_ggir_all     <- rbindlist(ls_ind_disc_dt_ggir)
ref_time_all          <- rbindlist(ls_ref_time_all)
ggir_time_all         <- rbindlist(ls_ggir_time_all)

fwrite(slp_estimate_ggir_all, file.path(GGIR_OUTPUT_PATH, "final_slp_estimate_all.csv"))
fwrite(ind_disc_ggir_all,     file.path(GGIR_OUTPUT_PATH, "final_ind_disc_dt_all.csv"))
fwrite(group_disc_ggir_all,   file.path(GGIR_OUTPUT_PATH, "final_group_disc_all.csv"))
fwrite(ref_time_all,          file.path(GGIR_OUTPUT_PATH, "final_ref_time_all.csv"))
fwrite(ggir_time_all,         file.path(GGIR_OUTPUT_PATH, "final_ggir_time_all.csv"))

# ----------------------------------------------------------
# Step 6h: Compute sleep onset timing summary per algorithm
# ----------------------------------------------------------
# Convert Epoch timestamps to POSIXct (Eastern Time) and extract the
# clock hour of the first sleep epoch. Summarize mean, SD, min, max
# across subjects per algorithm, then combine with the reference summary.
ref_time_all[,  time_ref  := as.POSIXct(time_ref,  format = "%Y-%m-%d %H:%M:%OS", tz = "America/New_York")]
ggir_time_all[, time_ggir := as.POSIXct(time_ggir, format = "%Y-%m-%d %H:%M:%OS", tz = "America/New_York")]

ref_time_all[,  hour_ref  := hour(time_ref)]
ggir_time_all[, hour_ggir := hour(time_ggir)]

ggir_onset_summary <- ggir_time_all[, .(
  avg_hour = mean(hour_ggir), sd_hour = sd(hour_ggir),
  max_hour = max(hour_ggir),  min_hour = min(hour_ggir)
), by = "analysis_type"]
ggir_onset_summary[, analysis_type := word(analysis_type, 4, sep = "_")]

ref_onset_summary <- ref_time_all[, .(
  avg_hour = mean(hour_ref),  sd_hour = sd(hour_ref),
  max_hour = max(hour_ref),   min_hour = min(hour_ref)
), by = "analysis_type"][1]
ref_onset_summary[, analysis_type := "Reference"]

onset_summary_all <- rbind(ggir_onset_summary, ref_onset_summary)
fwrite(onset_summary_all, file.path(GGIR_OUTPUT_PATH, "ggir_summary_sleeponset.csv"))

# ----------------------------------------------------------
# Step 6i: Format GGIR group discrepancy table for publication
# ----------------------------------------------------------
# Merges LOA values with their confidence intervals into a single column,
# removes internal metadata columns, and maps algorithm codes to full names.
group_disc_ggir_pub <- copy(group_disc_ggir_all)
group_disc_ggir_pub[, LOA.lower := paste0(LOA.lower, " ", LOA.lower_CI)]
group_disc_ggir_pub[, LOA.upper := paste0(LOA.upper, " ", LOA.upper_CI)]
group_disc_ggir_pub[, `:=`(
  LOA.upper_CI = NULL, LOA.lower_CI = NULL,
  CI.type      = NULL, CI.level     = NULL,
  logTransf    = NULL
)]
group_disc_ggir_pub[, Algorithm := word(analysis_type, 4, sep = "_")]
group_disc_ggir_pub[Algorithm == "CK",      Algorithm := "Cole-Kripke"]
group_disc_ggir_pub[Algorithm == "vanHees", Algorithm := "van Hees"]

fwrite(group_disc_ggir_pub, file.path(GGIR_OUTPUT_PATH, "final_group_disc_all_publication.csv"))
message("GGIR analysis complete. Results saved to: ", GGIR_OUTPUT_PATH)


# ============================================================
# SECTION 7: MACHINE LEARNING MODEL ANALYSIS
# ============================================================
# Evaluate seven ML sleep detection models:
#   lstm    = Long Short-Term Memory neural network
#   rf      = Random Forest
#   nn_mlp  = Neural Network (Multi-Layer Perceptron)
#   knn     = k-Nearest Neighbor
#   xgb     = Extreme Gradient Boosting
#   log_reg = Logistic Regression
#   tree    = Decision Tree
#
# For each model:
#   1. Load epoch-by-epoch prediction files (predicted class + true label)
#   2. Recode labels: 0 -> Wake (10), 1 -> Sleep (5)
#   3. Run the full de Zambotti pipeline via run_dezombotti_pipeline()
#   4. Aggregate results across models for combined export
#
# NOTE: match_ids computed in Section 6 is reused here. Both GGIR and
# ML analyses operate on the same set of matched subjects.

ml_algorithms <- c("lstm", "rf", "nn_mlp", "knn", "xgb", "log_reg", "tree")

# Initialise aggregation lists
ls_group_disc_ml   <- list()
ls_slp_estimate_ml <- list()
ls_ind_disc_dt_ml  <- list()

for (algo_idx in seq_along(ml_algorithms)) {
  algo_name <- ml_algorithms[algo_idx]
  ml_path   <- file.path(DATA_PATH, paste0("ML_outputs_", algo_name), "subject_fold")
  ml_files  <- list.files(ml_path, pattern = "_pred_")

  message("Processing ML model: ", algo_name)

  # ----------------------------------------------------------
  # Step 7a: Load epoch-by-epoch ML predictions for each subject
  # ----------------------------------------------------------
  # ML prediction files: CSV format
  # Required columns:
  #   .pred_class : integer, model-predicted class (0 = Wake, 1 = Sleep)
  #   label_n     : integer, true reference label  (0 = Wake, 1 = Sleep)
  ls_all_data <- list()

  for (file_idx in seq_along(match_ids)) {
    ID      <- match_ids[file_idx]
    ml_file <- ml_files[grep(paste0("pred_", ID, "_"), ml_files)]

    # Load ML prediction CSV
    dt <- fread(file.path(ml_path, ml_file))

    # Select subject, device prediction, and reference label
    dt_new <- dt[, .(
      subject   = ID,
      device    = .pred_class,
      reference = label_n
    )]

    ls_all_data[[file_idx]] <- dt_new
    rm(dt, dt_new)
  }

  # ----------------------------------------------------------
  # Step 7b: Combine all subjects into a single data.table
  # ----------------------------------------------------------
  all_data <- setDT(rbindlist(ls_all_data))

  # ----------------------------------------------------------
  # Step 7c: Recode numeric ML labels to de Zambotti stage codes
  # ----------------------------------------------------------
  # ML outputs: 0 = Wake, 1 = Sleep
  # de Zambotti codes: wake = 10, sleep = 5
  all_data[reference == "",  reference := 0]   # treat empty as Wake
  all_data[device    == 0L,  device    := 10]
  all_data[device    == 1L,  device    := 5]
  all_data[reference == 0L,  reference := 10]
  all_data[reference == 1L,  reference := 5]

  # ----------------------------------------------------------
  # Step 7d: Run the de Zambotti pipeline for this ML model
  # ----------------------------------------------------------
  results_dir  <- file.path(ML_OUTPUT_PATH, "results", algo_name)
  pipeline_out <- run_dezombotti_pipeline(all_data, results_dir, algo_name)

  # ----------------------------------------------------------
  # Step 7e: Tag results with model name and store for aggregation
  # ----------------------------------------------------------
  setDT(pipeline_out$group_discr)[,   analysis_type := algo_name]
  setDT(pipeline_out$sleep_metrics)[, analysis_type := algo_name]
  setDT(pipeline_out$ind_discr)[,     analysis_type := algo_name]

  ls_group_disc_ml[[algo_idx]]   <- pipeline_out$group_discr
  ls_slp_estimate_ml[[algo_idx]] <- pipeline_out$sleep_metrics
  ls_ind_disc_dt_ml[[algo_idx]]  <- pipeline_out$ind_discr
}

# ----------------------------------------------------------
# Step 7f: Aggregate ML results across all models
# ----------------------------------------------------------
# Combine per-model results and save as aggregated CSV files.
group_disc_ml_all   <- rbindlist(ls_group_disc_ml)
slp_estimate_ml_all <- rbindlist(ls_slp_estimate_ml)
ind_disc_ml_all     <- rbindlist(ls_ind_disc_dt_ml)

fwrite(slp_estimate_ml_all, file.path(ML_OUTPUT_PATH, "final_slp_estimate_all.csv"))
fwrite(ind_disc_ml_all,     file.path(ML_OUTPUT_PATH, "final_ind_disc_dt_all.csv"))
fwrite(group_disc_ml_all,   file.path(ML_OUTPUT_PATH, "final_group_disc_all.csv"))

# ----------------------------------------------------------
# Step 7g: Format ML group discrepancy table for publication
# ----------------------------------------------------------
# Merges LOA values with CIs, removes metadata columns,
# and maps short algorithm codes to full descriptive names.
group_disc_ml_pub <- copy(group_disc_ml_all)
group_disc_ml_pub[, LOA.lower := paste0(LOA.lower, " ", LOA.lower_CI)]
group_disc_ml_pub[, LOA.upper := paste0(LOA.upper, " ", LOA.upper_CI)]
group_disc_ml_pub[, `:=`(
  LOA.upper_CI = NULL, LOA.lower_CI = NULL,
  CI.type      = NULL, CI.level     = NULL,
  logTransf    = NULL
)]

# Map short model names to publication-ready full names
algo_name_map <- c(
  rf      = "Random Forest",
  nn_mlp  = "Neural Network (MLP)",
  knn     = "k-Nearest Neighbor",
  xgb     = "Extreme Gradient Boosting",
  log_reg = "Logistic Regression",
  tree    = "Decision Tree",
  lstm    = "Long Short-Term Memory"
)
group_disc_ml_pub[, Algorithm := algo_name_map[analysis_type]]

fwrite(group_disc_ml_pub, file.path(ML_OUTPUT_PATH, "final_group_disc_all_publication.csv"))

message("ML analysis complete. Results saved to: ", ML_OUTPUT_PATH)
message("All analyses finished. All outputs in: ", OUTPUT_PATH)
