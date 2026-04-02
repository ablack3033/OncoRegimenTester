#' Parameter loader: loads and validates CSV parameter files
#'
#' @description
#' Each loader reads one CSV, checks required columns, validates data
#' constraints (e.g. probability normalization), and returns a data.table.
#' The aggregate loader [loadParameterPack()] loads everything at once.

#' @import data.table

# Tolerance for probability sums
.PROB_SUM_TOL <- 0.05


#' Read a CSV parameter file
#' @param path File path
#' @return data.table
#' @keywords internal
readParamCsv <- function(path) {
  if (!file.exists(path)) {
    stop("Parameter file not found: ", path, call. = FALSE)
  }
  dt <- data.table::fread(path, strip.white = TRUE)
  setnames(dt, trimws(names(dt)))
  dt
}


#' Check that required columns exist in a data.table
#' @keywords internal
checkColumns <- function(dt, required, filename) {
  missingCols <- setdiff(required, names(dt))
  if (length(missingCols) > 0) {
    stop(
      filename, " missing required columns: ",
      paste(sort(missingCols), collapse = ", "),
      call. = FALSE
    )
  }
}


#' Check that probabilities sum to ~1 within groups
#' @return Character vector of warnings
#' @keywords internal
checkProbabilityGroups <- function(dt, groupCols, probCol, filename) {
  warningsOut <- character(0)
  groups <- dt[, .(total = sum(get(probCol))), by = groupCols]
  bad <- groups[abs(total - 1.0) > .PROB_SUM_TOL]
  if (nrow(bad) > 0) {
    for (i in seq_len(nrow(bad))) {
      keyStr <- paste(
        paste0(groupCols, "=", as.character(bad[i, ..groupCols])),
        collapse = ", "
      )
      warningsOut <- c(
        warningsOut,
        sprintf("%s: group (%s) has %s sum = %.4f",
                filename, keyStr, probCol, bad$total[i])
      )
    }
  }
  warningsOut
}


#' Parse boolean-like values from CSV
#' @keywords internal
parseBool <- function(x) {
  tolower(trimws(as.character(x))) %in% c("true", "1", "yes", "t")
}


# ---- Individual loaders ----

#' @keywords internal
loadManifest <- function(paramDir) {
  path <- file.path(paramDir, "manifest.csv")
  dt <- readParamCsv(path)
  required <- c("simulator_version", "cohort_id", "cohort_name",
                 "window_version", "min_cell_count", "generation_noise_mode", "notes")
  checkColumns(dt, required, "manifest.csv")
  dt
}


#' @keywords internal
loadCohortPopulation <- function(paramDir) {
  path <- file.path(paramDir, "cohort_population.csv")
  dt <- readParamCsv(path)
  required <- c("cohort_id", "age_group", "sex", "obs_length_bin", "proportion")
  checkColumns(dt, required, "cohort_population.csv")
  dt[, proportion := as.numeric(proportion)]
  warningsOut <- checkProbabilityGroups(dt, "cohort_id", "proportion", "cohort_population.csv")
  list(data = dt, warnings = warningsOut)
}


#' @keywords internal
loadTreatmentProfileMix <- function(paramDir) {
  path <- file.path(paramDir, "treatment_profile_mix.csv")
  dt <- readParamCsv(path)
  required <- c("cohort_id", "age_group", "sex", "profile_id", "profile_label", "proportion")
  checkColumns(dt, required, "treatment_profile_mix.csv")
  dt[, proportion := as.numeric(proportion)]
  warningsOut <- checkProbabilityGroups(
    dt, c("cohort_id", "age_group", "sex"), "proportion", "treatment_profile_mix.csv"
  )
  list(data = dt, warnings = warningsOut)
}


#' @keywords internal
loadDrugFeatureCatalog <- function(paramDir) {
  path <- file.path(paramDir, "drug_feature_catalog.csv")
  dt <- readParamCsv(path)
  required <- c("drug_feature_id", "drug_feature_name", "drug_class",
                 "is_antineoplastic", "is_supportive", "route_group", "concept_set_name")
  checkColumns(dt, required, "drug_feature_catalog.csv")
  dt[, is_antineoplastic := parseBool(is_antineoplastic)]
  dt[, is_supportive := parseBool(is_supportive)]
  dupes <- dt[duplicated(drug_feature_id), drug_feature_id]
  if (length(dupes) > 0) {
    stop("drug_feature_catalog.csv has duplicate drug_feature_id: ",
         paste(dupes, collapse = ", "), call. = FALSE)
  }
  dt
}


#' @keywords internal
loadRegimenCatalog <- function(paramDir) {
  path <- file.path(paramDir, "regimen_catalog.csv")
  dt <- readParamCsv(path)
  required <- c("cohort_id", "regimen_id", "regimen_name", "line_preference",
                 "regimen_type", "is_maintenance", "notes")
  checkColumns(dt, required, "regimen_catalog.csv")
  dt[, is_maintenance := parseBool(is_maintenance)]
  dupes <- dt[duplicated(regimen_id), regimen_id]
  if (length(dupes) > 0) {
    stop("regimen_catalog.csv has duplicate regimen_id: ",
         paste(dupes, collapse = ", "), call. = FALSE)
  }
  dt
}


#' @keywords internal
loadRegimenDrugMap <- function(paramDir) {
  path <- file.path(paramDir, "regimen_drug_map.csv")
  dt <- readParamCsv(path)
  required <- c("regimen_id", "drug_feature_id", "drug_role", "required_flag")
  checkColumns(dt, required, "regimen_drug_map.csv")
  dt[, required_flag := parseBool(required_flag)]
  dt
}


#' @keywords internal
loadProfileRegimenStart <- function(paramDir) {
  path <- file.path(paramDir, "profile_regimen_start.csv")
  dt <- readParamCsv(path)
  required <- c("cohort_id", "profile_id", "line_number", "regimen_id", "probability")
  checkColumns(dt, required, "profile_regimen_start.csv")
  dt[, line_number := as.integer(line_number)]
  dt[, probability := as.numeric(probability)]
  warningsOut <- checkProbabilityGroups(
    dt, c("cohort_id", "profile_id", "line_number"), "probability",
    "profile_regimen_start.csv"
  )
  list(data = dt, warnings = warningsOut)
}


#' @keywords internal
loadRegimenTransition <- function(paramDir) {
  path <- file.path(paramDir, "regimen_transition.csv")
  dt <- readParamCsv(path)
  required <- c("cohort_id", "profile_id", "from_regimen_id",
                 "to_regimen_id", "transition_type", "probability")
  checkColumns(dt, required, "regimen_transition.csv")
  dt[, probability := as.numeric(probability)]
  warningsOut <- checkProbabilityGroups(
    dt, c("cohort_id", "profile_id", "from_regimen_id"), "probability",
    "regimen_transition.csv"
  )
  list(data = dt, warnings = warningsOut)
}


#' @keywords internal
loadRegimenDuration <- function(paramDir) {
  path <- file.path(paramDir, "regimen_duration.csv")
  dt <- readParamCsv(path)
  required <- c("cohort_id", "profile_id", "regimen_id", "duration_bin", "probability")
  checkColumns(dt, required, "regimen_duration.csv")
  dt[, probability := as.numeric(probability)]
  warningsOut <- checkProbabilityGroups(
    dt, c("cohort_id", "profile_id", "regimen_id"), "probability",
    "regimen_duration.csv"
  )
  list(data = dt, warnings = warningsOut)
}


#' @keywords internal
loadDrugSchedule <- function(paramDir) {
  path <- file.path(paramDir, "drug_schedule.csv")
  dt <- readParamCsv(path)
  required <- c("regimen_id", "drug_feature_id", "cycle_length_days",
                 "admin_day_pattern", "exposure_type", "repeat_within_cycle_flag")
  checkColumns(dt, required, "drug_schedule.csv")
  dt[, cycle_length_days := as.integer(cycle_length_days)]
  dt[, repeat_within_cycle_flag := parseBool(repeat_within_cycle_flag)]
  dt
}


#' @keywords internal
loadDrugDropAddRules <- function(paramDir) {
  path <- file.path(paramDir, "drug_drop_add_rules.csv")
  if (!file.exists(path)) return(data.table())
  dt <- readParamCsv(path)
  required <- c("cohort_id", "profile_id", "regimen_id",
                 "modification_type", "drug_feature_id", "probability")
  checkColumns(dt, required, "drug_drop_add_rules.csv")
  dt[, probability := as.numeric(probability)]
  dt
}


#' @keywords internal
loadSupportiveDrugRules <- function(paramDir) {
  path <- file.path(paramDir, "supportive_drug_rules.csv")
  if (!file.exists(path)) return(data.table())
  dt <- readParamCsv(path)
  required <- c("regimen_id", "drug_feature_id", "probability",
                 "ignore_for_regimen_detection_flag")
  checkColumns(dt, required, "supportive_drug_rules.csv")
  dt[, probability := as.numeric(probability)]
  dt[, ignore_for_regimen_detection_flag := parseBool(ignore_for_regimen_detection_flag)]
  dt
}


#' @keywords internal
loadConditionSignal <- function(paramDir) {
  path <- file.path(paramDir, "condition_signal.csv")
  if (!file.exists(path)) return(data.table())
  dt <- readParamCsv(path)
  required <- c("cohort_id", "profile_id", "window_id", "condition_feature_id", "prevalence")
  checkColumns(dt, required, "condition_signal.csv")
  dt[, prevalence := as.numeric(prevalence)]
  dt
}


#' @keywords internal
loadProcedureSignal <- function(paramDir) {
  path <- file.path(paramDir, "procedure_signal.csv")
  if (!file.exists(path)) return(data.table())
  dt <- readParamCsv(path)
  required <- c("cohort_id", "profile_id", "window_id", "procedure_feature_id", "prevalence")
  checkColumns(dt, required, "procedure_signal.csv")
  dt[, prevalence := as.numeric(prevalence)]
  dt
}


#' @keywords internal
loadMeasurementSignal <- function(paramDir) {
  path <- file.path(paramDir, "measurement_signal.csv")
  if (!file.exists(path)) return(data.table())
  dt <- readParamCsv(path)
  required <- c("cohort_id", "profile_id", "window_id", "measurement_feature_id", "prevalence")
  checkColumns(dt, required, "measurement_signal.csv")
  dt[, prevalence := as.numeric(prevalence)]
  dt
}


#' @keywords internal
loadMeasurementValueDist <- function(paramDir) {
  path <- file.path(paramDir, "measurement_value_dist.csv")
  if (!file.exists(path)) return(list(data = data.table(), warnings = character(0)))
  dt <- readParamCsv(path)
  required <- c("cohort_id", "profile_id", "measurement_feature_id",
                 "value_bin", "probability")
  checkColumns(dt, required, "measurement_value_dist.csv")
  dt[, probability := as.numeric(probability)]
  warningsOut <- checkProbabilityGroups(
    dt, c("cohort_id", "profile_id", "measurement_feature_id"), "probability",
    "measurement_value_dist.csv"
  )
  list(data = dt, warnings = warningsOut)
}


# ---- Aggregate loader ----

#' Load an entire parameter pack
#'
#' Reads all CSV parameter files from a directory, validates columns and
#' probability normalization, and checks referential integrity.
#'
#' @param paramDir Path to directory containing CSV parameter files
#' @param cohortId Optional cohort ID to filter to
#' @return A list with class "parameterPack" containing all parameter tables
#'   and a `warnings` vector
#' @export
loadParameterPack <- function(paramDir, cohortId = NULL) {
  pack <- list()
  pack$paramDir <- paramDir
  pack$warnings <- character(0)

  pack$manifest <- loadManifest(paramDir)

  res <- loadCohortPopulation(paramDir)
  pack$cohortPopulation <- res$data
  pack$warnings <- c(pack$warnings, res$warnings)

  res <- loadTreatmentProfileMix(paramDir)
  pack$treatmentProfileMix <- res$data
  pack$warnings <- c(pack$warnings, res$warnings)

  pack$drugFeatureCatalog <- loadDrugFeatureCatalog(paramDir)
  pack$regimenCatalog <- loadRegimenCatalog(paramDir)
  pack$regimenDrugMap <- loadRegimenDrugMap(paramDir)

  res <- loadProfileRegimenStart(paramDir)
  pack$profileRegimenStart <- res$data
  pack$warnings <- c(pack$warnings, res$warnings)

  res <- loadRegimenTransition(paramDir)
  pack$regimenTransition <- res$data
  pack$warnings <- c(pack$warnings, res$warnings)

  res <- loadRegimenDuration(paramDir)
  pack$regimenDuration <- res$data
  pack$warnings <- c(pack$warnings, res$warnings)

  pack$drugSchedule <- loadDrugSchedule(paramDir)
  pack$drugDropAddRules <- loadDrugDropAddRules(paramDir)
  pack$supportiveDrugRules <- loadSupportiveDrugRules(paramDir)
  pack$conditionSignal <- loadConditionSignal(paramDir)
  pack$procedureSignal <- loadProcedureSignal(paramDir)
  pack$measurementSignal <- loadMeasurementSignal(paramDir)

  res <- loadMeasurementValueDist(paramDir)
  pack$measurementValueDist <- res$data
  pack$warnings <- c(pack$warnings, res$warnings)

  if (!is.null(cohortId)) {
    pack <- filterPackByCohort(pack, cohortId)
  }

  pack$warnings <- c(pack$warnings, validateForeignKeys(pack))

  class(pack) <- "parameterPack"
  pack
}


#' Filter a parameter pack to a single cohort
#' @keywords internal
filterPackByCohort <- function(pack, cohortId) {
  filterCol <- function(dt, col = "cohort_id") {
    if (nrow(dt) == 0 || !col %in% names(dt)) return(dt)
    dt[get(col) == cohortId]
  }
  pack$cohortPopulation <- filterCol(pack$cohortPopulation)
  pack$treatmentProfileMix <- filterCol(pack$treatmentProfileMix)
  pack$regimenCatalog <- filterCol(pack$regimenCatalog)
  pack$profileRegimenStart <- filterCol(pack$profileRegimenStart)
  pack$regimenTransition <- filterCol(pack$regimenTransition)
  pack$regimenDuration <- filterCol(pack$regimenDuration)
  pack$drugDropAddRules <- filterCol(pack$drugDropAddRules)
  pack$conditionSignal <- filterCol(pack$conditionSignal)
  pack$procedureSignal <- filterCol(pack$procedureSignal)
  pack$measurementSignal <- filterCol(pack$measurementSignal)
  pack$measurementValueDist <- filterCol(pack$measurementValueDist)
  pack
}


#' Validate foreign key references across tables
#' @return Character vector of warnings
#' @keywords internal
validateForeignKeys <- function(pack) {
  warningsOut <- character(0)

  drugIds <- pack$drugFeatureCatalog$drug_feature_id
  regimenIds <- pack$regimenCatalog$regimen_id
  profileIds <- unique(pack$treatmentProfileMix$profile_id)

  if (nrow(pack$regimenDrugMap) > 0) {
    badReg <- setdiff(pack$regimenDrugMap$regimen_id, regimenIds)
    for (r in badReg) {
      warningsOut <- c(warningsOut,
        sprintf("regimen_drug_map: unknown regimen_id '%s'", r))
    }
    badDrug <- setdiff(pack$regimenDrugMap$drug_feature_id, drugIds)
    for (d in badDrug) {
      warningsOut <- c(warningsOut,
        sprintf("regimen_drug_map: unknown drug_feature_id '%s'", d))
    }
  }

  if (nrow(pack$drugSchedule) > 0) {
    badReg <- setdiff(pack$drugSchedule$regimen_id, regimenIds)
    for (r in badReg) {
      warningsOut <- c(warningsOut,
        sprintf("drug_schedule: unknown regimen_id '%s'", r))
    }
    badDrug <- setdiff(pack$drugSchedule$drug_feature_id, drugIds)
    for (d in badDrug) {
      warningsOut <- c(warningsOut,
        sprintf("drug_schedule: unknown drug_feature_id '%s'", d))
    }
  }

  if (nrow(pack$profileRegimenStart) > 0) {
    badProf <- setdiff(pack$profileRegimenStart$profile_id, profileIds)
    for (p in badProf) {
      warningsOut <- c(warningsOut,
        sprintf("profile_regimen_start: unknown profile_id '%s'", p))
    }
    badReg <- setdiff(pack$profileRegimenStart$regimen_id, regimenIds)
    for (r in badReg) {
      warningsOut <- c(warningsOut,
        sprintf("profile_regimen_start: unknown regimen_id '%s'", r))
    }
  }

  warningsOut
}


#' Validate a loaded parameter pack
#'
#' Runs all validation checks and returns a summary.
#'
#' @param pack A parameterPack object
#' @return A list with `valid` (logical), `warnings` (character vector),
#'   and `summary` (character vector)
#' @export
validateParameterPack <- function(pack) {
  summaryLines <- character(0)

  summaryLines <- c(summaryLines,
    sprintf("Manifest rows: %d", nrow(pack$manifest)),
    sprintf("Cohort population rows: %d", nrow(pack$cohortPopulation)),
    sprintf("Treatment profiles: %d", nrow(pack$treatmentProfileMix)),
    sprintf("Drug catalog entries: %d", nrow(pack$drugFeatureCatalog)),
    sprintf("Regimen catalog entries: %d", nrow(pack$regimenCatalog)),
    sprintf("Regimen-drug mappings: %d", nrow(pack$regimenDrugMap)),
    sprintf("Profile regimen starts: %d", nrow(pack$profileRegimenStart)),
    sprintf("Regimen transitions: %d", nrow(pack$regimenTransition)),
    sprintf("Regimen durations: %d", nrow(pack$regimenDuration)),
    sprintf("Drug schedules: %d", nrow(pack$drugSchedule))
  )

  list(
    valid = length(pack$warnings) == 0,
    warnings = pack$warnings,
    summary = summaryLines
  )
}
