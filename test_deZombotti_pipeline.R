# ============================================================
# test_deZombotti_pipeline.R
#
# Unit Tests for deZombotti_sleep_estimate.R
#
# Description:
#   Tests the core de Zambotti framework functions using
#   synthetic epoch-by-epoch data. No real participant data
#   is required to run these tests.
#
# Framework: testthat (https://testthat.r-lib.org/)
#
# Usage:
#   1. Install testthat if not already installed:
#        install.packages("testthat")
#   2. Set FUNCTIONS_PATH below to your local Functions/ folder.
#   3. Run the entire script, or run individual test blocks.
#      In RStudio, you can also use: testthat::test_file("test_deZombotti_pipeline.R")
#
# Input:  None (synthetic data is generated internally)
# Output: Pass/fail messages printed to the R console
# ============================================================

library(testthat)
library(data.table)

# ---- Load Framework Functions ----
# Update FUNCTIONS_PATH to point to the Functions/ folder of the
# sleep-trackers-performance repository.
FUNCTIONS_PATH <- file.path(dirname(rstudioapi::getActiveDocumentContext()$path), "Functions")
# Or set manually:
# FUNCTIONS_PATH <- "/path/to/sleep-trackers-performance/Functions"

r_files <- list.files(FUNCTIONS_PATH, pattern = "\\.R$", full.names = TRUE)
if (length(r_files) == 0) stop("Functions not found at: ", FUNCTIONS_PATH)
invisible(lapply(r_files, source))


# ============================================================
# HELPER: make_test_data()
# ============================================================
# Generates a minimal synthetic epoch-by-epoch dataset.
#
# Arguments:
#   n_subjects : number of simulated subjects
#   n_epochs   : number of 30-second epochs per subject
#   seed       : random seed for reproducibility
#   sleep_prop : proportion of epochs classified as Sleep (stage = 5)
#
# Returns:
#   data.table with columns: subject (int), device (int), reference (int)
#   Stage codes: 5 = Sleep, 10 = Wake (de Zambotti convention)
make_test_data <- function(n_subjects = 2,
                           n_epochs   = 120,
                           seed       = 42,
                           sleep_prop = 0.70) {
  set.seed(seed)
  data.table(
    subject   = rep(seq_len(n_subjects), each = n_epochs),
    device    = sample(c(5L, 10L), n_subjects * n_epochs,
                       replace = TRUE, prob = c(sleep_prop, 1 - sleep_prop)),
    reference = sample(c(5L, 10L), n_subjects * n_epochs,
                       replace = TRUE, prob = c(sleep_prop, 1 - sleep_prop))
  )
}

# Convenience: a dataset where every epoch is Sleep (edge case)
make_all_sleep_data <- function(n_subjects = 2, n_epochs = 60) {
  data.table(
    subject   = rep(seq_len(n_subjects), each = n_epochs),
    device    = rep(5L, n_subjects * n_epochs),
    reference = rep(5L, n_subjects * n_epochs)
  )
}

# Convenience: a dataset where every epoch is Wake (edge case)
make_all_wake_data <- function(n_subjects = 2, n_epochs = 60) {
  data.table(
    subject   = rep(seq_len(n_subjects), each = n_epochs),
    device    = rep(10L, n_subjects * n_epochs),
    reference = rep(10L, n_subjects * n_epochs)
  )
}


# ============================================================
# TEST SUITE 1: ebe2sleep()
# ============================================================
# Tests that ebe2sleep correctly converts epoch-by-epoch data
# into per-subject sleep metric summaries.

test_that("ebe2sleep: returns one row per subject", {
  td <- make_test_data(n_subjects = 3)
  result <- ebe2sleep(data = td, idCol = "subject", RefCol = "reference",
                      deviceCol = "device", epochLenght = 30,
                      staging = FALSE, stages = c(wake = 10, sleep = 5), digits = 2)
  expect_equal(nrow(result), 3)
})

test_that("ebe2sleep: output contains all required columns (non-staging mode)", {
  td <- make_test_data()
  result <- ebe2sleep(data = td, idCol = "subject", RefCol = "reference",
                      deviceCol = "device", epochLenght = 30,
                      staging = FALSE, stages = c(wake = 10, sleep = 5), digits = 2)
  expected_cols <- c("subject", "TIB", "TST_ref", "TST_device",
                     "SE_ref", "SE_device", "SOL_ref", "SOL_device",
                     "WASO_ref", "WASO_device")
  expect_true(all(expected_cols %in% colnames(result)))
})

test_that("ebe2sleep: TIB equals n_epochs * epoch_length / 60 (minutes)", {
  # 120 epochs * 30 s / 60 = 60 minutes
  td <- make_test_data(n_subjects = 1, n_epochs = 120)
  result <- ebe2sleep(data = td, idCol = "subject", RefCol = "reference",
                      deviceCol = "device", epochLenght = 30,
                      staging = FALSE, stages = c(wake = 10, sleep = 5), digits = 2)
  expect_equal(result$TIB, 60)
})

test_that("ebe2sleep: SE is between 0 and 100 for both device and reference", {
  td <- make_test_data(n_subjects = 5)
  result <- ebe2sleep(data = td, idCol = "subject", RefCol = "reference",
                      deviceCol = "device", epochLenght = 30,
                      staging = FALSE, stages = c(wake = 10, sleep = 5), digits = 2)
  expect_true(all(result$SE_ref    >= 0 & result$SE_ref    <= 100))
  expect_true(all(result$SE_device >= 0 & result$SE_device <= 100))
})

test_that("ebe2sleep: TST equals TIB when all epochs are Sleep", {
  td <- make_all_sleep_data(n_subjects = 1, n_epochs = 60)
  result <- ebe2sleep(data = td, idCol = "subject", RefCol = "reference",
                      deviceCol = "device", epochLenght = 30,
                      staging = FALSE, stages = c(wake = 10, sleep = 5), digits = 2)
  expect_equal(result$TST_ref,    result$TIB)
  expect_equal(result$TST_device, result$TIB)
})

test_that("ebe2sleep: TST equals 0 when all epochs are Wake", {
  td <- make_all_wake_data(n_subjects = 1, n_epochs = 60)
  result <- ebe2sleep(data = td, idCol = "subject", RefCol = "reference",
                      deviceCol = "device", epochLenght = 30,
                      staging = FALSE, stages = c(wake = 10, sleep = 5), digits = 2)
  expect_equal(result$TST_ref,    0)
  expect_equal(result$TST_device, 0)
})

test_that("ebe2sleep: SOL equals 0 when recording starts with Sleep", {
  td <- data.table(
    subject   = rep(1L, 60),
    device    = rep(5L, 60),
    reference = rep(5L, 60)
  )
  result <- ebe2sleep(data = td, idCol = "subject", RefCol = "reference",
                      deviceCol = "device", epochLenght = 30,
                      staging = FALSE, stages = c(wake = 10, sleep = 5), digits = 2)
  expect_equal(result$SOL_ref,    0)
  expect_equal(result$SOL_device, 0)
})


# ============================================================
# TEST SUITE 2: indDiscr()
# ============================================================
# Tests that individual discrepancy values are computed correctly.

test_that("indDiscr: returns one row per subject", {
  td     <- make_test_data(n_subjects = 4)
  sm     <- ebe2sleep(data = td, idCol = "subject", RefCol = "reference",
                      deviceCol = "device", epochLenght = 30,
                      staging = FALSE, stages = c(wake = 10, sleep = 5), digits = 2)
  result <- indDiscr(data = sm, staging = FALSE, digits = 2, doPlot = FALSE)
  expect_equal(nrow(result), 4)
})

test_that("indDiscr: TST_diff equals TST_device minus TST_ref", {
  td  <- make_test_data(n_subjects = 3)
  sm  <- ebe2sleep(data = td, idCol = "subject", RefCol = "reference",
                   deviceCol = "device", epochLenght = 30,
                   staging = FALSE, stages = c(wake = 10, sleep = 5), digits = 2)
  res <- indDiscr(data = sm, staging = FALSE, digits = 2, doPlot = FALSE)
  expected <- round(sm$TST_device - sm$TST_ref, 2)
  expect_equal(res$TST_diff, expected)
})

test_that("indDiscr: SE_diff equals SE_device minus SE_ref", {
  td  <- make_test_data(n_subjects = 3)
  sm  <- ebe2sleep(data = td, idCol = "subject", RefCol = "reference",
                   deviceCol = "device", epochLenght = 30,
                   staging = FALSE, stages = c(wake = 10, sleep = 5), digits = 2)
  res <- indDiscr(data = sm, staging = FALSE, digits = 2, doPlot = FALSE)
  expected <- round(sm$SE_device - sm$SE_ref, 2)
  expect_equal(res$SE_diff, expected)
})

test_that("indDiscr: discrepancies are zero when device equals reference", {
  n <- 100
  td <- data.table(
    subject   = rep(1:2, each = n),
    device    = rep(c(rep(5L, 70), rep(10L, 30)), 2),
    reference = rep(c(rep(5L, 70), rep(10L, 30)), 2)
  )
  sm  <- ebe2sleep(data = td, idCol = "subject", RefCol = "reference",
                   deviceCol = "device", epochLenght = 30,
                   staging = FALSE, stages = c(wake = 10, sleep = 5), digits = 2)
  res <- indDiscr(data = sm, staging = FALSE, digits = 2, doPlot = FALSE)
  expect_true(all(res$TST_diff == 0))
  expect_true(all(res$SE_diff  == 0))
})


# ============================================================
# TEST SUITE 3: errorMatrix()
# ============================================================
# Tests the confusion matrix computation.

test_that("errorMatrix: device_tot row sums to total epoch count", {
  n_epochs <- 100
  td <- make_test_data(n_subjects = 1, n_epochs = n_epochs)
  result <- errorMatrix(
    data = td, idCol = "subject", RefCol = "reference", deviceCol = "device",
    staging = FALSE, stages = c(wake = 10, sleep = 5), matrixType = "sum",
    CI.type = "classic", boot.type = "basic", CI.level = .95, digits = 2
  )
  device_row   <- result[result$reference == "device_tot", ]
  total_epochs <- as.numeric(device_row$device_wake) + as.numeric(device_row$device_sleep)
  expect_equal(total_epochs, n_epochs)
})

test_that("errorMatrix: reference_tot column sums correctly", {
  n_epochs <- 80
  td <- make_test_data(n_subjects = 1, n_epochs = n_epochs)
  result <- errorMatrix(
    data = td, idCol = "subject", RefCol = "reference", deviceCol = "device",
    staging = FALSE, stages = c(wake = 10, sleep = 5), matrixType = "sum",
    CI.type = "classic", boot.type = "basic", CI.level = .95, digits = 2
  )
  ref_rows  <- result[result$reference != "device_tot", ]
  total_ref <- sum(as.numeric(ref_rows$reference_tot))
  expect_equal(total_ref, n_epochs)
})


# ============================================================
# TEST SUITE 4: indEBE()
# ============================================================
# Tests epoch-by-epoch metrics per subject.

test_that("indEBE: accuracy is between 0 and 100 for all subjects", {
  td     <- make_test_data(n_subjects = 4, n_epochs = 60)
  result <- indEBE(data = td, stage = 5, stageLabel = "sleep",
                   doPlot = FALSE, digits = 2)
  expect_true(all(result$accuracy >= 0 & result$accuracy <= 100))
})

test_that("indEBE: sensitivity and specificity are between 0 and 100", {
  td     <- make_test_data(n_subjects = 3, n_epochs = 80)
  result <- indEBE(data = td, stage = 5, stageLabel = "sleep",
                   doPlot = FALSE, digits = 2)
  expect_true(all(result$sensitivity >= 0 & result$sensitivity <= 100))
  expect_true(all(result$specificity  >= 0 & result$specificity  <= 100))
})

test_that("indEBE: returns one row per subject", {
  td     <- make_test_data(n_subjects = 5, n_epochs = 60)
  result <- indEBE(data = td, stage = 5, stageLabel = "sleep",
                   doPlot = FALSE, digits = 2)
  expect_equal(nrow(result), 5)
})

test_that("indEBE: perfect accuracy when device equals reference", {
  n <- 100
  td <- data.table(
    subject   = rep(1:3, each = n),
    device    = rep(c(rep(5L, 70), rep(10L, 30)), 3),
    reference = rep(c(rep(5L, 70), rep(10L, 30)), 3)
  )
  result <- indEBE(data = td, stage = 5, stageLabel = "sleep",
                   doPlot = FALSE, digits = 2)
  expect_true(all(result$accuracy    == 100))
  expect_true(all(result$sensitivity == 100))
  expect_true(all(result$specificity == 100))
})

test_that("indEBE: balanced_acc equals (sensitivity + specificity) / 2", {
  td     <- make_test_data(n_subjects = 3, n_epochs = 80)
  result <- indEBE(data = td, stage = 5, stageLabel = "sleep",
                   doPlot = FALSE, digits = 2)
  expected_balanced <- round((result$sensitivity + result$specificity) / 2, 2)
  expect_equal(result$balanced_acc, expected_balanced)
})


# ============================================================
# TEST SUITE 5: Label recoding validation
# ============================================================
# Tests the Sleep/Wake -> 5/10 recoding logic used in the pipeline.

test_that("Label recoding: Sleep/Wake strings convert to 5/10 correctly", {
  td <- data.table(
    subject   = rep(1:2, each = 10),
    device    = rep(c("Sleep", "Wake"), 10),
    reference = rep(c("Wake", "Sleep"), 10)
  )
  td[device    == "Wake",  device    := 10]
  td[device    == "Sleep", device    := 5]
  td[reference == "Wake",  reference := 10]
  td[reference == "Sleep", reference := 5]
  expect_true(all(td$device    %in% c("5", "10")))
  expect_true(all(td$reference %in% c("5", "10")))
})

test_that("Label recoding: numeric 0/1 converts to 10/5 correctly", {
  td <- data.table(
    subject   = rep(1:2, each = 10),
    device    = rep(c(0L, 1L), 10),
    reference = rep(c(1L, 0L), 10)
  )
  td[device    == 0L, device    := 10L]
  td[device    == 1L, device    := 5L]
  td[reference == 0L, reference := 10L]
  td[reference == 1L, reference := 5L]
  expect_true(all(td$device    %in% c(5L, 10L)))
  expect_true(all(td$reference %in% c(5L, 10L)))
})

test_that("Label recoding: empty string reference defaults to Wake (10)", {
  td <- data.table(
    subject   = 1:4,
    device    = c(5L, 10L, 5L, 10L),
    reference = c("Sleep", "", "Wake", "")
  )
  td[reference %in% "", reference := "Wake"]
  td[reference == "Wake",  reference := 10]
  td[reference == "Sleep", reference := 5]
  expect_true(all(td$reference %in% c("5", "10")))
  expect_equal(td$reference[c(2, 4)], c("10", "10"))
})


# ============================================================
# Run all tests and report
# ============================================================
message("\n============================================")
message("Running all unit tests for deZombotti pipeline...")
message("============================================\n")

test_results <- testthat::test_file(
  path     = rstudioapi::getActiveDocumentContext()$path,
  reporter = "summary"
)

message("\nAll tests complete.")
