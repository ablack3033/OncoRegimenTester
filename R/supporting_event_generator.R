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
windowToDayRange <- function(windowId, linePlans) {
  if (windowId == "pre_index") {
    return(c(-90L, 0L))
  }
  if (grepl("^line_", windowId)) {
    lineNum <- as.integer(sub("^line_", "", windowId))
    matching <- linePlans[linePlans$line_number == lineNum, ]
    if (nrow(matching) > 0) {
      return(c(matching$line_start_day[1], matching$line_end_day[1]))
    }
    return(c(0L, 90L))
  }
  if (windowId == "post_treatment") {
    if (nrow(linePlans) > 0) {
      lastEnd <- max(linePlans$line_end_day)
      return(c(lastEnd, lastEnd + 180L))
    }
    return(c(0L, 180L))
  }
  c(0L, 90L)
}


#' Generate all supporting events for all patients
#'
#' @param patients data.table from generatePatients
#' @param linePlans data.table from generateTherapy
#' @param pack parameterPack
#' @param config simulatorConfig
#' @return A list with elements: conditions, procedures, measurements
#' @export
generateSupportingEvents <- function(patients, linePlans, pack, config) {
  set.seed(config$seed + 2L)

  allConditions <- list()
  allProcedures <- list()
  allMeasurements <- list()
  ci <- 0L; pi <- 0L; mi <- 0L

  for (i in seq_len(nrow(patients))) {
    pat <- patients[i]
    patLines <- linePlans[linePlans$patient_id == pat$patient_id, ]

    if (config$generateConditions && nrow(pack$conditionSignal) > 0) {
      signals <- pack$conditionSignal[
        pack$conditionSignal$cohort_id == pat$cohort_id &
        pack$conditionSignal$profile_id == pat$profile_id,
      ]
      for (j in seq_len(nrow(signals))) {
        sig <- signals[j]
        if (runif(1) < sig$prevalence) {
          dayRange <- windowToDayRange(sig$window_id, patLines)
          day <- sample(dayRange[1]:max(dayRange[1], dayRange[2] - 1L), 1)
          ci <- ci + 1L
          allConditions[[ci]] <- data.table(
            patient_id = pat$patient_id,
            condition_feature_id = sig$condition_feature_id,
            condition_start_day = day,
            condition_end_day = NA_integer_
          )
        }
      }
    }

    if (config$generateProcedures && nrow(pack$procedureSignal) > 0) {
      signals <- pack$procedureSignal[
        pack$procedureSignal$cohort_id == pat$cohort_id &
        pack$procedureSignal$profile_id == pat$profile_id,
      ]
      for (j in seq_len(nrow(signals))) {
        sig <- signals[j]
        if (runif(1) < sig$prevalence) {
          dayRange <- windowToDayRange(sig$window_id, patLines)
          day <- sample(dayRange[1]:max(dayRange[1], dayRange[2] - 1L), 1)
          pi <- pi + 1L
          allProcedures[[pi]] <- data.table(
            patient_id = pat$patient_id,
            procedure_feature_id = sig$procedure_feature_id,
            procedure_day = day
          )
        }
      }
    }

    if (config$generateMeasurements && nrow(pack$measurementSignal) > 0) {
      signals <- pack$measurementSignal[
        pack$measurementSignal$cohort_id == pat$cohort_id &
        pack$measurementSignal$profile_id == pat$profile_id,
      ]
      for (j in seq_len(nrow(signals))) {
        sig <- signals[j]
        if (runif(1) < sig$prevalence) {
          dayRange <- windowToDayRange(sig$window_id, patLines)
          day <- sample(dayRange[1]:max(dayRange[1], dayRange[2] - 1L), 1)

          valueBin <- NA_character_
          if (nrow(pack$measurementValueDist) > 0) {
            vd <- pack$measurementValueDist[
              pack$measurementValueDist$cohort_id == pat$cohort_id &
              pack$measurementValueDist$profile_id == pat$profile_id &
              pack$measurementValueDist$measurement_feature_id == sig$measurement_feature_id,
            ]
            if (nrow(vd) > 0) {
              p <- as.numeric(vd$probability)
              p <- p / sum(p)
              valueBin <- sample(vd$value_bin, 1, prob = p)
            }
          }

          mi <- mi + 1L
          allMeasurements[[mi]] <- data.table(
            patient_id = pat$patient_id,
            measurement_feature_id = sig$measurement_feature_id,
            measurement_day = day,
            value_bin = valueBin
          )
        }
      }
    }
  }

  emptyCond <- data.table(
    patient_id = integer(0), condition_feature_id = character(0),
    condition_start_day = integer(0), condition_end_day = integer(0)
  )
  emptyProc <- data.table(
    patient_id = integer(0), procedure_feature_id = character(0),
    procedure_day = integer(0)
  )
  emptyMeas <- data.table(
    patient_id = integer(0), measurement_feature_id = character(0),
    measurement_day = integer(0), value_bin = character(0)
  )

  list(
    conditions = if (length(allConditions) > 0) rbindlist(allConditions) else emptyCond,
    procedures = if (length(allProcedures) > 0) rbindlist(allProcedures) else emptyProc,
    measurements = if (length(allMeasurements) > 0) rbindlist(allMeasurements) else emptyMeas
  )
}
