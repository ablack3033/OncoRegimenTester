#' Patient generator
#'
#' @description
#' Generates synthetic patients from demographic and profile parameters.
#' Each patient gets demographic attributes sampled from cohort_population,
#' a treatment profile sampled from treatment_profile_mix, and an index date.

#' @import data.table

#' Sample one item from a weighted vector with optional jitter
#' @keywords internal
weighted_sample <- function(items, weights, noise) {
  w <- as.numeric(weights)
  if (sum(w) <= 0) {
    return(sample(items, 1))
  }
  w <- w / sum(w)
  w <- jitter_probabilities(w, noise)
  sample(items, size = 1, prob = w)
}


#' Parse an obs_length_bin like "365_730" to a concrete day count
#' @keywords internal
parse_obs_length_bin <- function(bin_label) {
  if (grepl("_plus$", bin_label)) {
    start <- as.integer(sub("_plus$", "", bin_label))
    return(start + sample.int(366, 1) - 1L)
  }
  parts <- as.integer(strsplit(bin_label, "_")[[1]])
  sample(parts[1]:parts[2], 1)
}


#' Generate synthetic patients
#'
#' @param pack A parameter_pack object
#' @param config A simulator_config object
#' @return A data.table of synthetic patients
#' @export
generate_patients <- function(pack, config) {
  set.seed(config$seed)

  pop <- pack$cohort_population
  profiles <- pack$treatment_profile_mix

  if (nrow(pop) == 0) stop("No cohort_population rows available", call. = FALSE)
  if (nrow(profiles) == 0) stop("No treatment_profile_mix rows available", call. = FALSE)

  cal_days <- as.integer(config$calendar_end - config$calendar_start)
  if (cal_days <= 0) stop("calendar_end must be after calendar_start", call. = FALSE)

  patients <- vector("list", config$n_patients)

  for (pid in seq_len(config$n_patients)) {
    # Sample demographic stratum
    demo_idx <- as.integer(weighted_sample(
      seq_len(nrow(pop)), pop$proportion, config$noise
    ))
    chosen_demo <- pop[demo_idx]

    # Sample treatment profile conditional on age_group and sex
    matching <- profiles[
      age_group == chosen_demo$age_group & sex == chosen_demo$sex
    ]
    if (nrow(matching) == 0) {
      # Fallback to all profiles for this cohort
      matching <- profiles[cohort_id == chosen_demo$cohort_id]
    }
    if (nrow(matching) == 0) {
      stop(sprintf("No treatment profiles for age_group=%s, sex=%s",
                   chosen_demo$age_group, chosen_demo$sex), call. = FALSE)
    }

    prof_idx <- as.integer(weighted_sample(
      seq_len(nrow(matching)), matching$proportion, config$noise
    ))
    chosen_prof <- matching[prof_idx]

    index_day <- sample.int(cal_days + 1L, 1) - 1L
    obs_length <- parse_obs_length_bin(chosen_demo$obs_length_bin)

    patients[[pid]] <- data.table(
      patient_id = pid,
      cohort_id = chosen_demo$cohort_id,
      age_group = chosen_demo$age_group,
      sex = chosen_demo$sex,
      obs_length_bin = chosen_demo$obs_length_bin,
      obs_length_days = obs_length,
      profile_id = chosen_prof$profile_id,
      profile_label = chosen_prof$profile_label,
      index_day = index_day
    )
  }

  rbindlist(patients)
}
