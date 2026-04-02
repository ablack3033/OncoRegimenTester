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

getArg <- function(flag, default = NULL) {
  idx <- which(args == flag)
  if (length(idx) == 0) return(default)
  if (idx >= length(args)) return(default)
  args[idx + 1]
}

hasFlag <- function(flag) {
  flag %in% args
}

library(OncoRegimenTester)

noise <- noiseConfig(
  probabilityJitterSd = as.numeric(getArg("--noise_prob", "0")),
  durationJitterDays = as.integer(getArg("--noise_dur", "0")),
  cycleJitterDays = as.integer(getArg("--noise_cycle", "0")),
  missingnessRate = as.numeric(getArg("--noise_miss", "0"))
)

config <- simulatorConfig(
  paramDir = getArg("--param_dir", "inst/example_params/mCRC"),
  outputDir = getArg("--output_dir", "output"),
  seed = as.integer(getArg("--seed", "42")),
  nPatients = as.integer(getArg("--n_patients", "100")),
  noise = noise,
  omopOutput = !hasFlag("--no_omop"),
  cohortId = getArg("--cohort_id")
)

result <- runSimulation(config)

message("\nDone. Output in: ", config$outputDir)
