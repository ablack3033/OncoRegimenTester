#' Parameter loader: loads and validates CSV parameter files
#'
#' @description
#' Each loader reads one CSV, checks required columns, validates data
#' constraints (e.g. probability normalization), and returns a data.table.
#' The aggregate loader [load_parameter_pack()] loads everything at once.

#' @import data.table

# Tolerance for probability sums
.PROB_SUM_TOL <- 0.05


#' Read a CSV parameter file
#' @param path File path
#' @return data.table
#' @keywords internal
read_param_csv <- function(path) {
  if (!file.exists(path)) {
    stop("Parameter file not found: ", path, call. = FALSE)
  }
  dt <- data.table::fread(path, strip.white = TRUE)
  # Trim column names
  setnames(dt, trimws(names(dt)))
  dt
}


#' Check that required columns exist in a data.table
#' @keywords internal
check_columns <- function(dt, required, filename) {
  missing_cols <- setdiff(required, names(dt))
  if (length(missing_cols) > 0) {
    stop(
      filename, " missing required columns: ",
      paste(sort(missing_cols), collapse = ", "),
      call. = FALSE
    )
  }
}


#' Check that probabilities sum to ~1 within groups
#' @return Character vector of warnings
#' @keywords internal
check_probability_groups <- function(dt, group_cols, prob_col, filename) {
  warnings_out <- character(0)
  groups <- dt[, .(total = sum(get(prob_col))), by = group_cols]
  bad <- groups[abs(total - 1.0) > .PROB_SUM_TOL]
  if (nrow(bad) > 0) {
    for (i in seq_len(nrow(bad))) {
      key_str <- paste(
        paste0(group_cols, "=", as.character(bad[i, ..group_cols])),
        collapse = ", "
      )
      warnings_out <- c(
        warnings_out,
        sprintf("%s: group (%s) has %s sum = %.4f",
                filename, key_str, prob_col, bad$total[i])
      )
    }
  }
  warnings_out
}


#' Parse boolean-like values from CSV
#' @keywords internal
parse_bool <- function(x) {
  tolower(trimws(as.character(x))) %in% c("true", "1", "yes", "t")
}


# ---- Individual loaders ----

#' @keywords internal
load_manifest <- function(param_dir) {
  path <- file.path(param_dir, "manifest.csv")
  dt <- read_param_csv(path)
  required <- c("simulator_version", "cohort_id", "cohort_name",
                 "window_version", "min_cell_count", "generation_noise_mode", "notes")
  check_columns(dt, required, "manifest.csv")
  dt
}


#' @keywords internal
load_cohort_population <- function(param_dir) {
  path <- file.path(param_dir, "cohort_population.csv")
  dt <- read_param_csv(path)
  required <- c("cohort_id", "age_group", "sex", "obs_length_bin", "proportion")
  check_columns(dt, required, "cohort_population.csv")
  dt[, proportion := as.numeric(proportion)]
  warnings_out <- check_probability_groups(dt, "cohort_id", "proportion", "cohort_population.csv")
  list(data = dt, warnings = warnings_out)
}


#' @keywords internal
load_treatment_profile_mix <- function(param_dir) {
  path <- file.path(param_dir, "treatment_profile_mix.csv")
  dt <- read_param_csv(path)
  required <- c("cohort_id", "age_group", "sex", "profile_id", "profile_label", "proportion")
  check_columns(dt, required, "treatment_profile_mix.csv")
  dt[, proportion := as.numeric(proportion)]
  warnings_out <- check_probability_groups(
    dt, c("cohort_id", "age_group", "sex"), "proportion", "treatment_profile_mix.csv"
  )
  list(data = dt, warnings = warnings_out)
}


#' @keywords internal
load_drug_feature_catalog <- function(param_dir) {
  path <- file.path(param_dir, "drug_feature_catalog.csv")
  dt <- read_param_csv(path)
  required <- c("drug_feature_id", "drug_feature_name", "drug_class",
                 "is_antineoplastic", "is_supportive", "route_group", "concept_set_name")
  check_columns(dt, required, "drug_feature_catalog.csv")
  dt[, is_antineoplastic := parse_bool(is_antineoplastic)]
  dt[, is_supportive := parse_bool(is_supportive)]
  # Check uniqueness
  dupes <- dt[duplicated(drug_feature_id), drug_feature_id]
  if (length(dupes) > 0) {
    stop("drug_feature_catalog.csv has duplicate drug_feature_id: ",
         paste(dupes, collapse = ", "), call. = FALSE)
  }
  dt
}


#' @keywords internal
load_regimen_catalog <- function(param_dir) {
  path <- file.path(param_dir, "regimen_catalog.csv")
  dt <- read_param_csv(path)
  required <- c("cohort_id", "regimen_id", "regimen_name", "line_preference",
                 "regimen_type", "is_maintenance", "notes")
  check_columns(dt, required, "regimen_catalog.csv")
  dt[, is_maintenance := parse_bool(is_maintenance)]
  dupes <- dt[duplicated(regimen_id), regimen_id]
  if (length(dupes) > 0) {
    stop("regimen_catalog.csv has duplicate regimen_id: ",
         paste(dupes, collapse = ", "), call. = FALSE)
  }
  dt
}


#' @keywords internal
load_regimen_drug_map <- function(param_dir) {
  path <- file.path(param_dir, "regimen_drug_map.csv")
  dt <- read_param_csv(path)
  required <- c("regimen_id", "drug_feature_id", "drug_role", "required_flag")
  check_columns(dt, required, "regimen_drug_map.csv")
  dt[, required_flag := parse_bool(required_flag)]
  dt
}


#' @keywords internal
load_profile_regimen_start <- function(param_dir) {
  path <- file.path(param_dir, "profile_regimen_start.csv")
  dt <- read_param_csv(path)
  required <- c("cohort_id", "profile_id", "line_number", "regimen_id", "probability")
  check_columns(dt, required, "profile_regimen_start.csv")
  dt[, line_number := as.integer(line_number)]
  dt[, probability := as.numeric(probability)]
  warnings_out <- check_probability_groups(
    dt, c("cohort_id", "profile_id", "line_number"), "probability",
    "profile_regimen_start.csv"
  )
  list(data = dt, warnings = warnings_out)
}


#' @keywords internal
load_regimen_transition <- function(param_dir) {
  path <- file.path(param_dir, "regimen_transition.csv")
  dt <- read_param_csv(path)
  required <- c("cohort_id", "profile_id", "from_regimen_id",
                 "to_regimen_id", "transition_type", "probability")
  check_columns(dt, required, "regimen_transition.csv")
  dt[, probability := as.numeric(probability)]
  warnings_out <- check_probability_groups(
    dt, c("cohort_id", "profile_id", "from_regimen_id"), "probability",
    "regimen_transition.csv"
  )
  list(data = dt, warnings = warnings_out)
}


#' @keywords internal
load_regimen_duration <- function(param_dir) {
  path <- file.path(param_dir, "regimen_duration.csv")
  dt <- read_param_csv(path)
  required <- c("cohort_id", "profile_id", "regimen_id", "duration_bin", "probability")
  check_columns(dt, required, "regimen_duration.csv")
  dt[, probability := as.numeric(probability)]
  warnings_out <- check_probability_groups(
    dt, c("cohort_id", "profile_id", "regimen_id"), "probability",
    "regimen_duration.csv"
  )
  list(data = dt, warnings = warnings_out)
}


#' @keywords internal
load_drug_schedule <- function(param_dir) {
  path <- file.path(param_dir, "drug_schedule.csv")
  dt <- read_param_csv(path)
  required <- c("regimen_id", "drug_feature_id", "cycle_length_days",
                 "admin_day_pattern", "exposure_type", "repeat_within_cycle_flag")
  check_columns(dt, required, "drug_schedule.csv")
  dt[, cycle_length_days := as.integer(cycle_length_days)]
  dt[, repeat_within_cycle_flag := parse_bool(repeat_within_cycle_flag)]
  dt
}


#' @keywords internal
load_drug_drop_add_rules <- function(param_dir) {
  path <- file.path(param_dir, "drug_drop_add_rules.csv")
  if (!file.exists(path)) return(data.table())
  dt <- read_param_csv(path)
  required <- c("cohort_id", "profile_id", "regimen_id",
                 "modification_type", "drug_feature_id", "probability")
  check_columns(dt, required, "drug_drop_add_rules.csv")
  dt[, probability := as.numeric(probability)]
  dt
}


#' @keywords internal
load_supportive_drug_rules <- function(param_dir) {
  path <- file.path(param_dir, "supportive_drug_rules.csv")
  if (!file.exists(path)) return(data.table())
  dt <- read_param_csv(path)
  required <- c("regimen_id", "drug_feature_id", "probability",
                 "ignore_for_regimen_detection_flag")
  check_columns(dt, required, "supportive_drug_rules.csv")
  dt[, probability := as.numeric(probability)]
  dt[, ignore_for_regimen_detection_flag := parse_bool(ignore_for_regimen_detection_flag)]
  dt
}


#' @keywords internal
load_condition_signal <- function(param_dir) {
  path <- file.path(param_dir, "condition_signal.csv")
  if (!file.exists(path)) return(data.table())
  dt <- read_param_csv(path)
  required <- c("cohort_id", "profile_id", "window_id", "condition_feature_id", "prevalence")
  check_columns(dt, required, "condition_signal.csv")
  dt[, prevalence := as.numeric(prevalence)]
  dt
}


#' @keywords internal
load_procedure_signal <- function(param_dir) {
  path <- file.path(param_dir, "procedure_signal.csv")
  if (!file.exists(path)) return(data.table())
  dt <- read_param_csv(path)
  required <- c("cohort_id", "profile_id", "window_id", "procedure_feature_id", "prevalence")
  check_columns(dt, required, "procedure_signal.csv")
  dt[, prevalence := as.numeric(prevalence)]
  dt
}


#' @keywords internal
load_measurement_signal <- function(param_dir) {
  path <- file.path(param_dir, "measurement_signal.csv")
  if (!file.exists(path)) return(data.table())
  dt <- read_param_csv(path)
  required <- c("cohort_id", "profile_id", "window_id", "measurement_feature_id", "prevalence")
  check_columns(dt, required, "measurement_signal.csv")
  dt[, prevalence := as.numeric(prevalence)]
  dt
}


#' @keywords internal
load_measurement_value_dist <- function(param_dir) {
  path <- file.path(param_dir, "measurement_value_dist.csv")
  if (!file.exists(path)) return(list(data = data.table(), warnings = character(0)))
  dt <- read_param_csv(path)
  required <- c("cohort_id", "profile_id", "measurement_feature_id",
                 "value_bin", "probability")
  check_columns(dt, required, "measurement_value_dist.csv")
  dt[, probability := as.numeric(probability)]
  warnings_out <- check_probability_groups(
    dt, c("cohort_id", "profile_id", "measurement_feature_id"), "probability",
    "measurement_value_dist.csv"
  )
  list(data = dt, warnings = warnings_out)
}


# ---- Aggregate loader ----

#' Load an entire parameter pack
#'
#' Reads all CSV parameter files from a directory, validates columns and
#' probability normalization, and checks referential integrity.
#'
#' @param param_dir Path to directory containing CSV parameter files
#' @param cohort_id Optional cohort ID to filter to
#' @return A list with class "parameter_pack" containing all parameter tables
#'   and a `warnings` vector
#' @export
load_parameter_pack <- function(param_dir, cohort_id = NULL) {
  pack <- list()
  pack$param_dir <- param_dir
  pack$warnings <- character(0)

  pack$manifest <- load_manifest(param_dir)

  res <- load_cohort_population(param_dir)
  pack$cohort_population <- res$data
  pack$warnings <- c(pack$warnings, res$warnings)

  res <- load_treatment_profile_mix(param_dir)
  pack$treatment_profile_mix <- res$data
  pack$warnings <- c(pack$warnings, res$warnings)

  pack$drug_feature_catalog <- load_drug_feature_catalog(param_dir)
  pack$regimen_catalog <- load_regimen_catalog(param_dir)
  pack$regimen_drug_map <- load_regimen_drug_map(param_dir)

  res <- load_profile_regimen_start(param_dir)
  pack$profile_regimen_start <- res$data
  pack$warnings <- c(pack$warnings, res$warnings)

  res <- load_regimen_transition(param_dir)
  pack$regimen_transition <- res$data
  pack$warnings <- c(pack$warnings, res$warnings)

  res <- load_regimen_duration(param_dir)
  pack$regimen_duration <- res$data
  pack$warnings <- c(pack$warnings, res$warnings)

  pack$drug_schedule <- load_drug_schedule(param_dir)
  pack$drug_drop_add_rules <- load_drug_drop_add_rules(param_dir)
  pack$supportive_drug_rules <- load_supportive_drug_rules(param_dir)
  pack$condition_signal <- load_condition_signal(param_dir)
  pack$procedure_signal <- load_procedure_signal(param_dir)
  pack$measurement_signal <- load_measurement_signal(param_dir)

  res <- load_measurement_value_dist(param_dir)
  pack$measurement_value_dist <- res$data
  pack$warnings <- c(pack$warnings, res$warnings)

  # Optional cohort filter

  if (!is.null(cohort_id)) {
    pack <- filter_pack_by_cohort(pack, cohort_id)
  }

  # Cross-reference validation
  pack$warnings <- c(pack$warnings, validate_foreign_keys(pack))

  class(pack) <- "parameter_pack"
  pack
}


#' Filter a parameter pack to a single cohort
#' @keywords internal
filter_pack_by_cohort <- function(pack, cohort_id) {
  filter_col <- function(dt, col = "cohort_id") {
    if (nrow(dt) == 0 || !col %in% names(dt)) return(dt)
    dt[get(col) == cohort_id]
  }
  pack$cohort_population <- filter_col(pack$cohort_population)
  pack$treatment_profile_mix <- filter_col(pack$treatment_profile_mix)
  pack$regimen_catalog <- filter_col(pack$regimen_catalog)
  pack$profile_regimen_start <- filter_col(pack$profile_regimen_start)
  pack$regimen_transition <- filter_col(pack$regimen_transition)
  pack$regimen_duration <- filter_col(pack$regimen_duration)
  pack$drug_drop_add_rules <- filter_col(pack$drug_drop_add_rules)
  pack$condition_signal <- filter_col(pack$condition_signal)
  pack$procedure_signal <- filter_col(pack$procedure_signal)
  pack$measurement_signal <- filter_col(pack$measurement_signal)
  pack$measurement_value_dist <- filter_col(pack$measurement_value_dist)
  pack
}


#' Validate foreign key references across tables
#' @return Character vector of warnings
#' @keywords internal
validate_foreign_keys <- function(pack) {
  warnings_out <- character(0)

  drug_ids <- pack$drug_feature_catalog$drug_feature_id
  regimen_ids <- pack$regimen_catalog$regimen_id
  profile_ids <- unique(pack$treatment_profile_mix$profile_id)

  # regimen_drug_map -> drugs and regimens
  if (nrow(pack$regimen_drug_map) > 0) {
    bad_reg <- setdiff(pack$regimen_drug_map$regimen_id, regimen_ids)
    for (r in bad_reg) {
      warnings_out <- c(warnings_out,
        sprintf("regimen_drug_map: unknown regimen_id '%s'", r))
    }
    bad_drug <- setdiff(pack$regimen_drug_map$drug_feature_id, drug_ids)
    for (d in bad_drug) {
      warnings_out <- c(warnings_out,
        sprintf("regimen_drug_map: unknown drug_feature_id '%s'", d))
    }
  }

  # drug_schedule -> regimens and drugs
  if (nrow(pack$drug_schedule) > 0) {
    bad_reg <- setdiff(pack$drug_schedule$regimen_id, regimen_ids)
    for (r in bad_reg) {
      warnings_out <- c(warnings_out,
        sprintf("drug_schedule: unknown regimen_id '%s'", r))
    }
    bad_drug <- setdiff(pack$drug_schedule$drug_feature_id, drug_ids)
    for (d in bad_drug) {
      warnings_out <- c(warnings_out,
        sprintf("drug_schedule: unknown drug_feature_id '%s'", d))
    }
  }

  # profile_regimen_start -> profiles and regimens
  if (nrow(pack$profile_regimen_start) > 0) {
    bad_prof <- setdiff(pack$profile_regimen_start$profile_id, profile_ids)
    for (p in bad_prof) {
      warnings_out <- c(warnings_out,
        sprintf("profile_regimen_start: unknown profile_id '%s'", p))
    }
    bad_reg <- setdiff(pack$profile_regimen_start$regimen_id, regimen_ids)
    for (r in bad_reg) {
      warnings_out <- c(warnings_out,
        sprintf("profile_regimen_start: unknown regimen_id '%s'", r))
    }
  }

  warnings_out
}


#' Validate a loaded parameter pack
#'
#' Runs all validation checks and returns a summary.
#'
#' @param pack A parameter_pack object
#' @return A list with `valid` (logical), `warnings` (character vector),
#'   and `summary` (character vector)
#' @export
validate_parameter_pack <- function(pack) {
  summary_lines <- character(0)

  summary_lines <- c(summary_lines,
    sprintf("Manifest rows: %d", nrow(pack$manifest)),
    sprintf("Cohort population rows: %d", nrow(pack$cohort_population)),
    sprintf("Treatment profiles: %d", nrow(pack$treatment_profile_mix)),
    sprintf("Drug catalog entries: %d", nrow(pack$drug_feature_catalog)),
    sprintf("Regimen catalog entries: %d", nrow(pack$regimen_catalog)),
    sprintf("Regimen-drug mappings: %d", nrow(pack$regimen_drug_map)),
    sprintf("Profile regimen starts: %d", nrow(pack$profile_regimen_start)),
    sprintf("Regimen transitions: %d", nrow(pack$regimen_transition)),
    sprintf("Regimen durations: %d", nrow(pack$regimen_duration)),
    sprintf("Drug schedules: %d", nrow(pack$drug_schedule))
  )

  list(
    valid = length(pack$warnings) == 0,
    warnings = pack$warnings,
    summary = summary_lines
  )
}
