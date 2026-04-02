#' OMOP-like output renderer
#'
#' @description
#' Renders simulation outputs as OMOP-like CSV tables.
#' Not full OMOP fidelity -- simplified for method testing.

#' @import data.table


#' Render all output tables and write to disk
#'
#' @param patients data.table of synthetic patients
#' @param therapy list from generate_therapy (line_plans, drug_exposures, ground_truth_episodes)
#' @param supporting list from generate_supporting_events (conditions, procedures, measurements)
#' @param detected list from detect_regimens (episodes, lines)
#' @param config simulator_config
#' @return Invisibly returns a list of all output tables
#' @export
render_omop <- function(patients, therapy, supporting, detected, config) {
  dir.create(config$output_dir, recursive = TRUE, showWarnings = FALSE)

  # -- person.csv --
  person <- data.table(
    person_id = patients$patient_id,
    gender_concept_id = ifelse(patients$sex == "M", 8507L, 8532L),
    year_of_birth = as.integer(format(config$calendar_start, "%Y")) -
      as.integer(gsub("\\D.*", "", patients$age_group)),
    cohort_id = patients$cohort_id,
    profile_id = patients$profile_id
  )

  # -- observation_period.csv --
  observation_period <- data.table(
    person_id = patients$patient_id,
    observation_period_start_date = config$calendar_start + patients$index_day,
    observation_period_end_date = config$calendar_start + patients$index_day +
      patients$obs_length_days
  )

  # -- drug_exposure.csv --
  de <- therapy$drug_exposures
  drug_exposure <- data.table()
  if (nrow(de) > 0) {
    drug_exposure <- data.table(
      drug_exposure_id = seq_len(nrow(de)),
      person_id = de$patient_id,
      drug_concept_id = de$drug_feature_id,
      drug_exposure_start_date = config$calendar_start +
        patients$index_day[match(de$patient_id, patients$patient_id)] +
        de$drug_exposure_start_day,
      drug_exposure_end_date = config$calendar_start +
        patients$index_day[match(de$patient_id, patients$patient_id)] +
        de$drug_exposure_end_day,
      regimen_id = de$regimen_id,
      line_number = de$line_number,
      is_supportive = de$is_supportive
    )
  }

  # -- condition_occurrence.csv --
  condition_occurrence <- data.table()
  if (nrow(supporting$conditions) > 0) {
    co <- supporting$conditions
    condition_occurrence <- data.table(
      condition_occurrence_id = seq_len(nrow(co)),
      person_id = co$patient_id,
      condition_concept_id = co$condition_feature_id,
      condition_start_date = config$calendar_start +
        patients$index_day[match(co$patient_id, patients$patient_id)] +
        co$condition_start_day,
      condition_end_date = NA
    )
  }

  # -- procedure_occurrence.csv --
  procedure_occurrence <- data.table()
  if (nrow(supporting$procedures) > 0) {
    po <- supporting$procedures
    procedure_occurrence <- data.table(
      procedure_occurrence_id = seq_len(nrow(po)),
      person_id = po$patient_id,
      procedure_concept_id = po$procedure_feature_id,
      procedure_date = config$calendar_start +
        patients$index_day[match(po$patient_id, patients$patient_id)] +
        po$procedure_day
    )
  }

  # -- measurement.csv --
  measurement <- data.table()
  if (nrow(supporting$measurements) > 0) {
    me <- supporting$measurements
    measurement <- data.table(
      measurement_id = seq_len(nrow(me)),
      person_id = me$patient_id,
      measurement_concept_id = me$measurement_feature_id,
      measurement_date = config$calendar_start +
        patients$index_day[match(me$patient_id, patients$patient_id)] +
        me$measurement_day,
      value_bin = me$value_bin
    )
  }

  # Write all tables
  fwrite(person, file.path(config$output_dir, "person.csv"))
  fwrite(observation_period, file.path(config$output_dir, "observation_period.csv"))
  if (nrow(drug_exposure) > 0) {
    fwrite(drug_exposure, file.path(config$output_dir, "drug_exposure.csv"))
  }
  if (nrow(condition_occurrence) > 0) {
    fwrite(condition_occurrence, file.path(config$output_dir, "condition_occurrence.csv"))
  }
  if (nrow(procedure_occurrence) > 0) {
    fwrite(procedure_occurrence, file.path(config$output_dir, "procedure_occurrence.csv"))
  }
  if (nrow(measurement) > 0) {
    fwrite(measurement, file.path(config$output_dir, "measurement.csv"))
  }

  # Ground truth and derived tables (always written)
  fwrite(therapy$line_plans, file.path(config$output_dir, "ground_truth_lines.csv"))
  fwrite(therapy$ground_truth_episodes,
         file.path(config$output_dir, "ground_truth_regimen_episodes.csv"))

  if (nrow(detected$episodes) > 0) {
    fwrite(detected$episodes, file.path(config$output_dir, "derived_regimen_episodes.csv"))
  }
  if (nrow(detected$lines) > 0) {
    fwrite(detected$lines, file.path(config$output_dir, "derived_lines.csv"))
  }

  message(sprintf("Output written to %s/", config$output_dir))

  invisible(list(
    person = person,
    observation_period = observation_period,
    drug_exposure = drug_exposure,
    condition_occurrence = condition_occurrence,
    procedure_occurrence = procedure_occurrence,
    measurement = measurement,
    ground_truth_lines = therapy$line_plans,
    ground_truth_regimen_episodes = therapy$ground_truth_episodes,
    derived_regimen_episodes = detected$episodes,
    derived_lines = detected$lines
  ))
}
