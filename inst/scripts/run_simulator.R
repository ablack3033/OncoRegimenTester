#!/usr/bin/env Rscript
# Command-line entrypoint for OncoRegimenTester
#
# Usage:
#   Rscript inst/scripts/run_simulator.R [options]
#
# Options:
#   --param_dir DIR     Directory with CSV parameter pack (default: inst/example_params/mCRC)
#   --output_dir DIR    Output directory (default: output)
#   --n_patients N      Number of patients (default: 100)
#   --seed N            Random seed (default: 42)
#   --no_omop           Skip OMOP table rendering
#   --noise_prob SD     Probability jitter SD (default: 0)
#   --noise_dur DAYS    Duration jitter days (default: 0)
#   --noise_cycle DAYS  Cycle jitter days (default: 0)
#   --noise_miss RATE   Missingness rate (default: 0)
#   --cohort_id ID      Filter to specific cohort

args <- commandArgs(trailingOnly = TRUE)

# Simple argument parser
get_arg <- function(flag, default = NULL) {
  idx <- which(args == flag)
  if (length(idx) == 0) return(default)
  if (idx >= length(args)) return(default)
  args[idx + 1]
}

has_flag <- function(flag) {
  flag %in% args
}

library(OncoRegimenTester)

noise <- noise_config(
  probability_jitter_sd = as.numeric(get_arg("--noise_prob", "0")),
  duration_jitter_days = as.integer(get_arg("--noise_dur", "0")),
  cycle_jitter_days = as.integer(get_arg("--noise_cycle", "0")),
  missingness_rate = as.numeric(get_arg("--noise_miss", "0"))
)

config <- simulator_config(
  param_dir = get_arg("--param_dir", "inst/example_params/mCRC"),
  output_dir = get_arg("--output_dir", "output"),
  seed = as.integer(get_arg("--seed", "42")),
  n_patients = as.integer(get_arg("--n_patients", "100")),
  noise = noise,
  omop_output = !has_flag("--no_omop"),
  cohort_id = get_arg("--cohort_id")
)

result <- run_simulation(config)

message("\nDone. Output in: ", config$output_dir)
