#' OMOP-like output renderer
#'
#' @description
#' Renders simulation outputs as OMOP-like CSV tables.
#' Not full OMOP fidelity -- simplified for method testing.
#'
#' Can write to disk as CSVs or insert into a CDMConnector CDM reference
#' object for direct use with OHDSI tooling.

#' @import data.table


#' Render all output tables and write to disk
#'
#' @param patients data.table of synthetic patients
#' @param therapy list from generateTherapy (linePlans, drugExposures, groundTruthEpisodes)
#' @param supporting list from generateSupportingEvents (conditions, procedures, measurements)
#' @param detected list from detectRegimens (episodes, lines)
#' @param config simulatorConfig
#' @return Invisibly returns a list of all output data.tables
#' @export
renderOmop <- function(patients, therapy, supporting, detected, config) {
  dir.create(config$outputDir, recursive = TRUE, showWarnings = FALSE)

  tables <- buildOmopTables(patients, therapy, supporting, detected, config)

  fwrite(tables$person, file.path(config$outputDir, "person.csv"))
  fwrite(tables$observation_period, file.path(config$outputDir, "observation_period.csv"))
  if (nrow(tables$drug_exposure) > 0) {
    fwrite(tables$drug_exposure, file.path(config$outputDir, "drug_exposure.csv"))
  }
  if (nrow(tables$condition_occurrence) > 0) {
    fwrite(tables$condition_occurrence, file.path(config$outputDir, "condition_occurrence.csv"))
  }
  if (nrow(tables$procedure_occurrence) > 0) {
    fwrite(tables$procedure_occurrence, file.path(config$outputDir, "procedure_occurrence.csv"))
  }
  if (nrow(tables$measurement) > 0) {
    fwrite(tables$measurement, file.path(config$outputDir, "measurement.csv"))
  }

  fwrite(therapy$linePlans, file.path(config$outputDir, "ground_truth_lines.csv"))
  fwrite(therapy$groundTruthEpisodes,
         file.path(config$outputDir, "ground_truth_regimen_episodes.csv"))

  if (nrow(detected$episodes) > 0) {
    fwrite(detected$episodes, file.path(config$outputDir, "derived_regimen_episodes.csv"))
  }
  if (nrow(detected$lines) > 0) {
    fwrite(detected$lines, file.path(config$outputDir, "derived_lines.csv"))
  }

  message(sprintf("Output written to %s/", config$outputDir))

  invisible(tables)
}


#' Build OMOP-like data.tables from simulation output
#'
#' @param patients data.table of synthetic patients
#' @param therapy list from generateTherapy
#' @param supporting list from generateSupportingEvents
#' @param detected list from detectRegimens
#' @param config simulatorConfig
#' @return Named list of OMOP-like data.tables
#' @keywords internal
buildOmopTables <- function(patients, therapy, supporting, detected, config) {
  person <- data.table(
    person_id = patients$patient_id,
    gender_concept_id = ifelse(patients$sex == "M", 8507L, 8532L),
    year_of_birth = as.integer(format(config$calendarStart, "%Y")) -
      as.integer(gsub("\\D.*", "", patients$age_group)),
    cohort_id = patients$cohort_id,
    profile_id = patients$profile_id
  )

  observationPeriod <- data.table(
    person_id = patients$patient_id,
    observation_period_start_date = config$calendarStart + patients$index_day,
    observation_period_end_date = config$calendarStart + patients$index_day +
      patients$obs_length_days
  )

  de <- therapy$drugExposures
  drugExposure <- data.table()
  if (nrow(de) > 0) {
    drugExposure <- data.table(
      drug_exposure_id = seq_len(nrow(de)),
      person_id = de$patient_id,
      drug_concept_id = de$drug_feature_id,
      drug_exposure_start_date = config$calendarStart +
        patients$index_day[match(de$patient_id, patients$patient_id)] +
        de$drug_exposure_start_day,
      drug_exposure_end_date = config$calendarStart +
        patients$index_day[match(de$patient_id, patients$patient_id)] +
        de$drug_exposure_end_day,
      regimen_id = de$regimen_id,
      line_number = de$line_number,
      is_supportive = de$is_supportive
    )
  }

  conditionOccurrence <- data.table()
  if (nrow(supporting$conditions) > 0) {
    co <- supporting$conditions
    conditionOccurrence <- data.table(
      condition_occurrence_id = seq_len(nrow(co)),
      person_id = co$patient_id,
      condition_concept_id = co$condition_feature_id,
      condition_start_date = config$calendarStart +
        patients$index_day[match(co$patient_id, patients$patient_id)] +
        co$condition_start_day,
      condition_end_date = NA
    )
  }

  procedureOccurrence <- data.table()
  if (nrow(supporting$procedures) > 0) {
    po <- supporting$procedures
    procedureOccurrence <- data.table(
      procedure_occurrence_id = seq_len(nrow(po)),
      person_id = po$patient_id,
      procedure_concept_id = po$procedure_feature_id,
      procedure_date = config$calendarStart +
        patients$index_day[match(po$patient_id, patients$patient_id)] +
        po$procedure_day
    )
  }

  meas <- data.table()
  if (nrow(supporting$measurements) > 0) {
    me <- supporting$measurements
    meas <- data.table(
      measurement_id = seq_len(nrow(me)),
      person_id = me$patient_id,
      measurement_concept_id = me$measurement_feature_id,
      measurement_date = config$calendarStart +
        patients$index_day[match(me$patient_id, patients$patient_id)] +
        me$measurement_day,
      value_bin = me$value_bin
    )
  }

  list(
    person = person,
    observation_period = observationPeriod,
    drug_exposure = drugExposure,
    condition_occurrence = conditionOccurrence,
    procedure_occurrence = procedureOccurrence,
    measurement = meas,
    ground_truth_lines = therapy$linePlans,
    ground_truth_regimen_episodes = therapy$groundTruthEpisodes,
    derived_regimen_episodes = detected$episodes,
    derived_lines = detected$lines
  )
}


#' Write simulation output into a CDMConnector CDM reference object
#'
#' Takes a CDMConnector CDM reference and inserts generated OMOP tables
#' into the database connection backing the CDM. This allows the synthetic
#' data to be used directly with OHDSI analysis tools.
#'
#' @param cdm A CDMConnector CDM reference object (from `CDMConnector::cdmFromCon()`)
#' @param patients data.table of synthetic patients
#' @param therapy list from generateTherapy
#' @param supporting list from generateSupportingEvents
#' @param detected list from detectRegimens
#' @param config simulatorConfig
#' @return The CDM reference object with synthetic tables inserted
#' @export
writeToCdm <- function(cdm, patients, therapy, supporting, detected, config) {
  if (!requireNamespace("CDMConnector", quietly = TRUE)) {
    stop("CDMConnector is required for writeToCdm(). ",
         "Install with: install.packages('CDMConnector')", call. = FALSE)
  }
  if (!requireNamespace("DBI", quietly = TRUE)) {
    stop("DBI is required for writeToCdm(). ",
         "Install with: install.packages('DBI')", call. = FALSE)
  }

  tables <- buildOmopTables(patients, therapy, supporting, detected, config)
  con <- attr(cdm, "dbcon")

  # CDMConnector stores the connection as an attribute; fall back to cdm_con
  if (is.null(con) && requireNamespace("CDMConnector", quietly = TRUE)) {
    con <- tryCatch(CDMConnector::cdmCon(cdm), error = function(e) NULL)
  }
  if (is.null(con)) {
    stop("Could not extract database connection from CDM object", call. = FALSE)
  }

  writeSchema <- attr(cdm, "write_schema")
  if (is.null(writeSchema)) {
    writeSchema <- tryCatch(CDMConnector::cdmWriteSchema(cdm), error = function(e) NULL)
  }

  omopTables <- c("person", "observation_period", "drug_exposure",
                   "condition_occurrence", "procedure_occurrence", "measurement")

  for (tblName in omopTables) {
    tbl <- tables[[tblName]]
    if (!is.null(tbl) && nrow(tbl) > 0) {
      targetName <- if (!is.null(writeSchema)) {
        paste0(writeSchema, ".", tblName)
      } else {
        tblName
      }
      DBI::dbWriteTable(con, tblName, as.data.frame(tbl),
                        overwrite = TRUE, temporary = FALSE)
    }
  }

  message("Synthetic OMOP tables written to CDM database connection")
  cdm
}
