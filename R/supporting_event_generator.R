#' Supporting event generator
#'
#' @description
#' Generates optional events: conditions, procedures, measurements.
#' These are secondary to drug exposures -- they align with treatment
#' windows to add plausibility but are not the primary testing signal.

#' @import data.table


#' Map a window_id to a day range relative to index
#'
#' Supported formats:
#' - "pre_index": -90 to 0
#' - "line_N": start/end of line N
#' - "post_treatment": after last line end
#' @keywords internal
window_to_day_range <- function(window_id, line_plans) {
  if (window_id == "pre_index") {
    return(c(-90L, 0L))
  }
  if (grepl("^line_", window_id)) {
    line_num <- as.integer(sub("^line_", "", window_id))
    matching <- line_plans[line_plans$line_number == line_num, ]
    if (nrow(matching) > 0) {
      return(c(matching$line_start_day[1], matching$line_end_day[1]))
    }
    return(c(0L, 90L))
  }
  if (window_id == "post_treatment") {
    if (nrow(line_plans) > 0) {
      last_end <- max(line_plans$line_end_day)
      return(c(last_end, last_end + 180L))
    }
    return(c(0L, 180L))
  }
  c(0L, 90L)
}


#' Generate all supporting events for all patients
#'
#' @param patients data.table from generate_patients
#' @param line_plans data.table from generate_therapy
#' @param pack parameter_pack
#' @param config simulator_config
#' @return A list with elements: conditions, procedures, measurements
#' @export
generate_supporting_events <- function(patients, line_plans, pack, config) {
  set.seed(config$seed + 2L)

  all_conditions <- list()
  all_procedures <- list()
  all_measurements <- list()
  ci <- 0L; pi <- 0L; mi <- 0L

  for (i in seq_len(nrow(patients))) {
    pat <- patients[i]
    pat_lines <- line_plans[line_plans$patient_id == pat$patient_id, ]

    # Conditions
    if (config$generate_conditions && nrow(pack$condition_signal) > 0) {
      signals <- pack$condition_signal[
        pack$condition_signal$cohort_id == pat$cohort_id &
        pack$condition_signal$profile_id == pat$profile_id,
      ]
      for (j in seq_len(nrow(signals))) {
        sig <- signals[j]
        if (runif(1) < sig$prevalence) {
          day_range <- window_to_day_range(sig$window_id, pat_lines)
          day <- sample(day_range[1]:max(day_range[1], day_range[2] - 1L), 1)
          ci <- ci + 1L
          all_conditions[[ci]] <- data.table(
            patient_id = pat$patient_id,
            condition_feature_id = sig$condition_feature_id,
            condition_start_day = day,
            condition_end_day = NA_integer_
          )
        }
      }
    }

    # Procedures
    if (config$generate_procedures && nrow(pack$procedure_signal) > 0) {
      signals <- pack$procedure_signal[
        pack$procedure_signal$cohort_id == pat$cohort_id &
        pack$procedure_signal$profile_id == pat$profile_id,
      ]
      for (j in seq_len(nrow(signals))) {
        sig <- signals[j]
        if (runif(1) < sig$prevalence) {
          day_range <- window_to_day_range(sig$window_id, pat_lines)
          day <- sample(day_range[1]:max(day_range[1], day_range[2] - 1L), 1)
          pi <- pi + 1L
          all_procedures[[pi]] <- data.table(
            patient_id = pat$patient_id,
            procedure_feature_id = sig$procedure_feature_id,
            procedure_day = day
          )
        }
      }
    }

    # Measurements
    if (config$generate_measurements && nrow(pack$measurement_signal) > 0) {
      signals <- pack$measurement_signal[
        pack$measurement_signal$cohort_id == pat$cohort_id &
        pack$measurement_signal$profile_id == pat$profile_id,
      ]
      for (j in seq_len(nrow(signals))) {
        sig <- signals[j]
        if (runif(1) < sig$prevalence) {
          day_range <- window_to_day_range(sig$window_id, pat_lines)
          day <- sample(day_range[1]:max(day_range[1], day_range[2] - 1L), 1)

          # Sample value bin if available
          value_bin <- NA_character_
          if (nrow(pack$measurement_value_dist) > 0) {
            vd <- pack$measurement_value_dist[
              pack$measurement_value_dist$cohort_id == pat$cohort_id &
              pack$measurement_value_dist$profile_id == pat$profile_id &
              pack$measurement_value_dist$measurement_feature_id == sig$measurement_feature_id,
            ]
            if (nrow(vd) > 0) {
              p <- as.numeric(vd$probability)
              p <- p / sum(p)
              value_bin <- sample(vd$value_bin, 1, prob = p)
            }
          }

          mi <- mi + 1L
          all_measurements[[mi]] <- data.table(
            patient_id = pat$patient_id,
            measurement_feature_id = sig$measurement_feature_id,
            measurement_day = day,
            value_bin = value_bin
          )
        }
      }
    }
  }

  empty_cond <- data.table(
    patient_id = integer(0), condition_feature_id = character(0),
    condition_start_day = integer(0), condition_end_day = integer(0)
  )
  empty_proc <- data.table(
    patient_id = integer(0), procedure_feature_id = character(0),
    procedure_day = integer(0)
  )
  empty_meas <- data.table(
    patient_id = integer(0), measurement_feature_id = character(0),
    measurement_day = integer(0), value_bin = character(0)
  )

  list(
    conditions = if (length(all_conditions) > 0) rbindlist(all_conditions) else empty_cond,
    procedures = if (length(all_procedures) > 0) rbindlist(all_procedures) else empty_proc,
    measurements = if (length(all_measurements) > 0) rbindlist(all_measurements) else empty_meas
  )
}
