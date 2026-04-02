#' Configuration for the OncoRegimenTester simulator
#'
#' @description
#' Creates configuration objects that control simulation behavior.
#' All settings are explicit and interpretable.

#' Create a noise configuration
#'
#' Controls generation-time stochastic jitter for realism.
#' This is NOT differential privacy -- it is simulation variability
#' so that outputs are not perfectly deterministic copies of the parameter tables.
#'
#' @param probability_jitter_sd SD on logit scale applied to sampling probabilities
#' @param duration_jitter_days Max +/- days added to sampled regimen durations
#' @param cycle_jitter_days Max +/- days shift on admin days within a cycle
#' @param transition_jitter_sd SD on logit scale for transition probabilities
#' @param missingness_rate Probability of dropping an individual drug exposure
#' @return A list with class "noise_config"
#' @export
noise_config <- function(probability_jitter_sd = 0,
                         duration_jitter_days = 0L,
                         cycle_jitter_days = 0L,
                         transition_jitter_sd = 0,
                         missingness_rate = 0) {
  cfg <- list(
    probability_jitter_sd = probability_jitter_sd,
    duration_jitter_days = as.integer(duration_jitter_days),
    cycle_jitter_days = as.integer(cycle_jitter_days),
    transition_jitter_sd = transition_jitter_sd,
    missingness_rate = missingness_rate
  )
  class(cfg) <- "noise_config"
  cfg
}


#' Create a simulator configuration
#'
#' @param param_dir Directory containing the CSV parameter pack
#' @param output_dir Directory for generated outputs
#' @param seed Master random seed for reproducibility
#' @param n_patients Number of synthetic patients to generate
#' @param calendar_start Earliest possible index date (ISO format)
#' @param calendar_end Latest possible index date (ISO format)
#' @param noise A noise_config object
#' @param generate_conditions Whether to generate condition records
#' @param generate_procedures Whether to generate procedure records
#' @param generate_measurements Whether to generate measurement records
#' @param omop_output Whether to render OMOP-like output tables
#' @param cohort_id If set, restrict to this cohort only
#' @return A list with class "simulator_config"
#' @export
simulator_config <- function(param_dir = "inst/example_params/mCRC",
                             output_dir = "output",
                             seed = 42L,
                             n_patients = 100L,
                             calendar_start = "2018-01-01",
                             calendar_end = "2023-12-31",
                             noise = noise_config(),
                             generate_conditions = TRUE,
                             generate_procedures = TRUE,
                             generate_measurements = TRUE,
                             omop_output = TRUE,
                             cohort_id = NULL) {
  cfg <- list(
    param_dir = param_dir,
    output_dir = output_dir,
    seed = as.integer(seed),
    n_patients = as.integer(n_patients),
    calendar_start = as.Date(calendar_start),
    calendar_end = as.Date(calendar_end),
    noise = noise,
    generate_conditions = generate_conditions,
    generate_procedures = generate_procedures,
    generate_measurements = generate_measurements,
    omop_output = omop_output,
    cohort_id = cohort_id
  )
  class(cfg) <- "simulator_config"
  cfg
}
