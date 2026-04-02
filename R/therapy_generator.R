#' Therapy generator: lines, regimens, and drug exposures
#'
#' @description
#' The core generative engine. For each patient it:
#' 1. Plans lines of therapy (which regimens, in what order)
#' 2. Applies regimen modifications (drug add/drop/substitute)
#' 3. Generates individual drug exposure records from schedules
#'
#' The ground truth (what regimen was intended) is preserved separately
#' from the drug exposures (what would appear in the data). This separation
#' is the whole point: regimen detection algorithms only see the exposures.

#' @import data.table


#' Parse a duration_bin like "29_84" to a concrete day count
#' @keywords internal
parse_duration_bin <- function(bin_label, noise) {
  if (grepl("_plus$", bin_label)) {
    start <- as.integer(sub("_plus$", "", bin_label))
    base <- start + sample.int(180, 1) - 1L
  } else {
    parts <- as.integer(strsplit(bin_label, "_")[[1]])
    base <- sample(parts[1]:parts[2], 1)
  }
  jitter_duration(base, noise)
}


#' Parse comma-separated admin_day_pattern into 0-indexed offsets
#' @keywords internal
parse_admin_day_pattern <- function(pattern) {
  # Days are 1-indexed in CSV, convert to 0-indexed
  as.integer(trimws(strsplit(pattern, ",")[[1]])) - 1L
}


#' Sample a starting regimen for a given line number
#' @keywords internal
sample_regimen_for_line <- function(pack, cohort_id, profile_id, line_number, noise) {
  cid <- cohort_id
  pid <- profile_id
  ln <- line_number
  candidates <- pack$profile_regimen_start[
    pack$profile_regimen_start$cohort_id == cid &
    pack$profile_regimen_start$profile_id == pid &
    pack$profile_regimen_start$line_number == ln,
  ]
  if (nrow(candidates) == 0) return(NULL)
  weighted_sample(candidates$regimen_id, candidates$probability, noise)
}


#' Sample a regimen duration in days
#' @keywords internal
sample_duration <- function(pack, cohort_id, profile_id, regimen_id, noise) {
  cid <- cohort_id; pid <- profile_id; rid <- regimen_id
  candidates <- pack$regimen_duration[
    pack$regimen_duration$cohort_id == cid &
    pack$regimen_duration$profile_id == pid &
    pack$regimen_duration$regimen_id == rid,
  ]
  if (nrow(candidates) == 0) {
    # Fallback: any duration row for this regimen
    candidates <- pack$regimen_duration[pack$regimen_duration$regimen_id == rid, ]
  }
  if (nrow(candidates) == 0) {
    return(sample(28:168, 1))  # default fallback
  }
  chosen_bin <- weighted_sample(candidates$duration_bin, candidates$probability, noise)
  parse_duration_bin(chosen_bin, noise)
}


#' Sample a transition from current regimen
#' @return Named list with transition_type and to_regimen_id
#' @keywords internal
sample_transition <- function(pack, cohort_id, profile_id, from_regimen_id, noise) {
  cid <- cohort_id; pid <- profile_id; frid <- from_regimen_id
  candidates <- pack$regimen_transition[
    pack$regimen_transition$cohort_id == cid &
    pack$regimen_transition$profile_id == pid &
    pack$regimen_transition$from_regimen_id == frid,
  ]
  if (nrow(candidates) == 0) {
    return(list(transition_type = "discontinue", to_regimen_id = ""))
  }
  w <- as.numeric(candidates$probability)
  w <- w / sum(w)
  w <- jitter_probabilities(w, noise)
  idx <- sample.int(nrow(candidates), 1, prob = w)
  list(
    transition_type = candidates$transition_type[idx],
    to_regimen_id = candidates$to_regimen_id[idx]
  )
}


#' Generate therapy plans and drug exposures for all patients
#'
#' @param patients data.table from generate_patients
#' @param pack parameter_pack
#' @param config simulator_config
#' @return A list with elements: line_plans, drug_exposures, ground_truth_episodes
#' @export
generate_therapy <- function(patients, pack, config) {
  set.seed(config$seed + 1L)  # offset seed for therapy generation

  all_lines <- vector("list", nrow(patients))
  all_exposures <- vector("list", nrow(patients))
  all_episodes <- vector("list", nrow(patients))

  for (i in seq_len(nrow(patients))) {
    pat <- patients[i]

    lines <- generate_line_plans_for_patient(pat, pack, config)
    exposures <- generate_drug_exposures_for_patient(pat, lines, pack, config)
    episodes <- build_ground_truth_episodes(pat, lines, pack)

    all_lines[[i]] <- lines
    all_exposures[[i]] <- exposures
    all_episodes[[i]] <- episodes
  }

  list(
    line_plans = rbindlist(all_lines),
    drug_exposures = rbindlist(all_exposures),
    ground_truth_episodes = rbindlist(all_episodes)
  )
}


#' Generate line plans for one patient
#' @keywords internal
generate_line_plans_for_patient <- function(pat, pack, config) {
  lines <- list()
  current_day <- 0L
  max_lines <- 4L

  for (line_num in seq_len(max_lines)) {
    regimen_id <- sample_regimen_for_line(
      pack, pat$cohort_id, pat$profile_id, line_num, config$noise
    )
    if (is.null(regimen_id)) break

    duration <- sample_duration(
      pack, pat$cohort_id, pat$profile_id, regimen_id, config$noise
    )

    line_start <- current_day
    line_end <- current_day + duration

    transition <- sample_transition(
      pack, pat$cohort_id, pat$profile_id, regimen_id, config$noise
    )

    lines[[line_num]] <- data.table(
      patient_id = pat$patient_id,
      line_number = line_num,
      regimen_id = regimen_id,
      line_start_day = line_start,
      line_end_day = line_end,
      transition_type = transition$transition_type
    )

    if (transition$transition_type == "discontinue") break

    # Gap between lines (7-28 days)
    gap <- sample(7L:28L, 1)
    current_day <- line_end + gap
  }

  if (length(lines) == 0) {
    return(data.table(
      patient_id = integer(0), line_number = integer(0),
      regimen_id = character(0), line_start_day = integer(0),
      line_end_day = integer(0), transition_type = character(0)
    ))
  }

  rbindlist(lines)
}


#' Get active drugs for a regimen after applying modifications
#' @keywords internal
get_active_drugs <- function(pack, regimen_id, pat) {
  rid <- regimen_id
  base_drugs <- pack$regimen_drug_map[pack$regimen_drug_map$regimen_id == rid, ]
  active_ids <- base_drugs$drug_feature_id

  # Apply modification rules
  cid <- pat$cohort_id; pid <- pat$profile_id
  mods <- pack$drug_drop_add_rules[
    pack$drug_drop_add_rules$regimen_id == rid &
    pack$drug_drop_add_rules$cohort_id == cid &
    pack$drug_drop_add_rules$profile_id == pid,
  ]

  if (nrow(mods) > 0) {
    for (j in seq_len(nrow(mods))) {
      mod <- mods[j]
      if (runif(1) < mod$probability) {
        if (mod$modification_type == "drop") {
          active_ids <- setdiff(active_ids, mod$drug_feature_id)
        } else if (mod$modification_type == "add") {
          active_ids <- union(active_ids, mod$drug_feature_id)
        } else if (mod$modification_type == "substitute") {
          # Drop first non-anchor and add substitute
          non_anchors <- base_drugs[
            drug_role != "anchor" & drug_feature_id %in% active_ids,
            drug_feature_id
          ]
          if (length(non_anchors) > 0) {
            active_ids <- setdiff(active_ids, non_anchors[1])
          }
          active_ids <- union(active_ids, mod$drug_feature_id)
        }
      }
    }
  }

  base_drugs[base_drugs$drug_feature_id %in% active_ids, ]
}


#' Generate drug exposures for one patient from their line plans
#' @export
generate_drug_exposures <- function(pat, line_plans, pack, config) {
  generate_drug_exposures_for_patient(pat, line_plans, pack, config)
}

#' @keywords internal
generate_drug_exposures_for_patient <- function(pat, line_plans, pack, config) {
  if (nrow(line_plans) == 0) {
    return(data.table(
      patient_id = integer(0), drug_feature_id = character(0),
      drug_exposure_start_day = integer(0), drug_exposure_end_day = integer(0),
      regimen_id = character(0), line_number = integer(0),
      is_supportive = logical(0)
    ))
  }

  exposures <- list()
  exp_idx <- 0L

  for (li in seq_len(nrow(line_plans))) {
    line <- line_plans[li]
    active_drugs <- get_active_drugs(pack, line$regimen_id, pat)

    # Build schedule lookup for this regimen
    rid <- line$regimen_id
    schedules <- pack$drug_schedule[pack$drug_schedule$regimen_id == rid, ]
    sched_map <- setNames(
      split(schedules, seq_len(nrow(schedules))),
      schedules$drug_feature_id
    )

    line_duration <- line$line_end_day - line$line_start_day

    for (di in seq_len(nrow(active_drugs))) {
      drug_id <- active_drugs$drug_feature_id[di]
      sched <- sched_map[[drug_id]]
      if (is.null(sched) || nrow(sched) == 0) next
      sched <- sched[1, ]  # take first if duplicates

      if (sched$exposure_type == "continuous") {
        # Continuous oral: one long exposure spanning the regimen
        if (!should_drop_exposure(config$noise)) {
          exp_idx <- exp_idx + 1L
          exposures[[exp_idx]] <- data.table(
            patient_id = pat$patient_id,
            drug_feature_id = drug_id,
            drug_exposure_start_day = line$line_start_day,
            drug_exposure_end_day = line$line_end_day,
            regimen_id = line$regimen_id,
            line_number = line$line_number,
            is_supportive = FALSE
          )
        }
      } else {
        # Cyclic administration
        admin_days <- parse_admin_day_pattern(sched$admin_day_pattern)
        cycle_len <- sched$cycle_length_days

        cycle_start <- 0L
        while (cycle_start < line_duration) {
          for (admin_day in admin_days) {
            actual_day <- cycle_start + admin_day
            actual_day <- jitter_cycle_day(actual_day, config$noise)

            if (actual_day >= line_duration) next

            abs_start <- line$line_start_day + actual_day

            # IV drugs ~1 day, infusions ~2 days
            exp_len <- if (sched$exposure_type == "infusion") 2L else 1L

            if (should_drop_exposure(config$noise)) next

            exp_idx <- exp_idx + 1L
            exposures[[exp_idx]] <- data.table(
              patient_id = pat$patient_id,
              drug_feature_id = drug_id,
              drug_exposure_start_day = abs_start,
              drug_exposure_end_day = abs_start + exp_len,
              regimen_id = line$regimen_id,
              line_number = line$line_number,
              is_supportive = FALSE
            )
          }
          cycle_start <- cycle_start + cycle_len
        }
      }
    }

    # Supportive drugs
    supp_rules <- pack$supportive_drug_rules[
      pack$supportive_drug_rules$regimen_id == rid,
    ]
    if (nrow(supp_rules) > 0) {
      for (si in seq_len(nrow(supp_rules))) {
        rule <- supp_rules[si]
        if (runif(1) < rule$probability) {
          n_admin <- max(1L, sample.int(max(2L, line_duration %/% 14L + 1L), 1))
          for (ai in seq_len(n_admin)) {
            day_offset <- sample.int(max(1L, line_duration), 1) - 1L
            abs_day <- line$line_start_day + day_offset

            if (should_drop_exposure(config$noise)) next

            exp_idx <- exp_idx + 1L
            exposures[[exp_idx]] <- data.table(
              patient_id = pat$patient_id,
              drug_feature_id = rule$drug_feature_id,
              drug_exposure_start_day = abs_day,
              drug_exposure_end_day = abs_day + 1L,
              regimen_id = line$regimen_id,
              line_number = line$line_number,
              is_supportive = TRUE
            )
          }
        }
      }
    }
  }

  if (length(exposures) == 0) {
    return(data.table(
      patient_id = integer(0), drug_feature_id = character(0),
      drug_exposure_start_day = integer(0), drug_exposure_end_day = integer(0),
      regimen_id = character(0), line_number = integer(0),
      is_supportive = logical(0)
    ))
  }

  rbindlist(exposures)
}


#' Build ground truth regimen episodes from line plans
#' @keywords internal
build_ground_truth_episodes <- function(pat, line_plans, pack) {
  if (nrow(line_plans) == 0) {
    return(data.table(
      patient_id = integer(0), episode_id = integer(0),
      regimen_id = character(0), regimen_name = character(0),
      episode_start_day = integer(0), episode_end_day = integer(0),
      line_number = integer(0)
    ))
  }

  regimen_map <- setNames(pack$regimen_catalog$regimen_name,
                          pack$regimen_catalog$regimen_id)

  data.table(
    patient_id = pat$patient_id,
    episode_id = seq_len(nrow(line_plans)),
    regimen_id = line_plans$regimen_id,
    regimen_name = ifelse(
      line_plans$regimen_id %in% names(regimen_map),
      regimen_map[line_plans$regimen_id],
      line_plans$regimen_id
    ),
    episode_start_day = line_plans$line_start_day,
    episode_end_day = line_plans$line_end_day,
    line_number = line_plans$line_number
  )
}
