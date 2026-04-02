#' Patient generator
#'
#' @description
#' Generates synthetic patients from demographic and profile parameters.
#' Each patient gets demographic attributes sampled from cohort_population,
#' a treatment profile sampled from treatment_profile_mix, and an index date.

#' @import data.table

#' Sample one item from a weighted vector with optional jitter
#' @keywords internal
weightedSample <- function(items, weights, noise) {
  w <- as.numeric(weights)
  if (sum(w) <= 0) {
    return(sample(items, 1))
  }
  w <- w / sum(w)
  w <- jitterProbabilities(w, noise)
  sample(items, size = 1, prob = w)
}


#' Parse an obs_length_bin like "365_730" to a concrete day count
#' @keywords internal
parseObsLengthBin <- function(binLabel) {
  if (grepl("_plus$", binLabel)) {
    start <- as.integer(sub("_plus$", "", binLabel))
    return(start + sample.int(366, 1) - 1L)
  }
  parts <- as.integer(strsplit(binLabel, "_")[[1]])
  if (length(parts) != 2 || anyNA(parts)) {
    warning("Invalid obs_length_bin: ", binLabel, ". Using default 365.")
    return(365L)
  }
  sample(parts[1]:parts[2], 1)
}


#' Generate synthetic patients
#'
#' @param pack A parameterPack object
#' @param config A simulatorConfig object
#' @return A data.table of synthetic patients
#' @export
generatePatients <- function(pack, config) {
  set.seed(config$seed)

  pop <- pack$cohortPopulation
  profiles <- pack$treatmentProfileMix

  if (nrow(pop) == 0) stop("No cohort_population rows available", call. = FALSE)
  if (nrow(profiles) == 0) stop("No treatment_profile_mix rows available", call. = FALSE)

  calDays <- as.integer(config$calendarEnd - config$calendarStart)
  if (calDays <= 0) stop("calendarEnd must be after calendarStart", call. = FALSE)

  patients <- vector("list", config$nPatients)

  for (pid in seq_len(config$nPatients)) {
    demoIdx <- as.integer(weightedSample(
      seq_len(nrow(pop)), pop$proportion, config$noise
    ))
    chosenDemo <- pop[demoIdx]

    matching <- profiles[
      age_group == chosenDemo$age_group & sex == chosenDemo$sex
    ]
    if (nrow(matching) == 0) {
      matching <- profiles[cohort_id == chosenDemo$cohort_id]
    }
    if (nrow(matching) == 0) {
      stop(sprintf("No treatment profiles for age_group=%s, sex=%s",
                   chosenDemo$age_group, chosenDemo$sex), call. = FALSE)
    }

    profIdx <- as.integer(weightedSample(
      seq_len(nrow(matching)), matching$proportion, config$noise
    ))
    chosenProf <- matching[profIdx]

    indexDay <- sample.int(calDays + 1L, 1) - 1L
    obsLength <- parseObsLengthBin(chosenDemo$obs_length_bin)

    patients[[pid]] <- data.table(
      patient_id = pid,
      cohort_id = chosenDemo$cohort_id,
      age_group = chosenDemo$age_group,
      sex = chosenDemo$sex,
      obs_length_bin = chosenDemo$obs_length_bin,
      obs_length_days = obsLength,
      profile_id = chosenProf$profile_id,
      profile_label = chosenProf$profile_label,
      index_day = indexDay
    )
  }

  rbindlist(patients)
}
