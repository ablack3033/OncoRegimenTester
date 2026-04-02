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
parseDurationBin <- function(binLabel, noise) {
  if (grepl("_plus$", binLabel)) {
    start <- as.integer(sub("_plus$", "", binLabel))
    base <- start + sample.int(180, 1) - 1L
  } else {
    parts <- as.integer(strsplit(binLabel, "_")[[1]])
    if (length(parts) != 2 || anyNA(parts)) {
      warning("Invalid duration_bin: ", binLabel, ". Using default 84.")
      return(jitterDuration(84L, noise))
    }
    base <- sample(parts[1]:parts[2], 1)
  }
  jitterDuration(base, noise)
}


#' Parse comma-separated admin_day_pattern into 0-indexed offsets
#' @keywords internal
parseAdminDayPattern <- function(pattern) {
  # Days are 1-indexed in CSV, convert to 0-indexed
  as.integer(trimws(strsplit(pattern, ",")[[1]])) - 1L
}


#' Sample a starting regimen for a given line number
#' @keywords internal
sampleRegimenForLine <- function(pack, cohortId, profileId, lineNumber, noise) {
  cid <- cohortId
  pid <- profileId
  ln <- lineNumber
  candidates <- pack$profileRegimenStart[
    pack$profileRegimenStart$cohort_id == cid &
    pack$profileRegimenStart$profile_id == pid &
    pack$profileRegimenStart$line_number == ln,
  ]
  if (nrow(candidates) == 0) return(NULL)
  weightedSample(candidates$regimen_id, candidates$probability, noise)
}


#' Sample a regimen duration in days
#' @keywords internal
sampleDuration <- function(pack, cohortId, profileId, regimenId, noise) {
  cid <- cohortId; pid <- profileId; rid <- regimenId
  candidates <- pack$regimenDuration[
    pack$regimenDuration$cohort_id == cid &
    pack$regimenDuration$profile_id == pid &
    pack$regimenDuration$regimen_id == rid,
  ]
  if (nrow(candidates) == 0) {
    candidates <- pack$regimenDuration[pack$regimenDuration$regimen_id == rid, ]
  }
  if (nrow(candidates) == 0) {
    return(sample(28:168, 1))
  }
  chosenBin <- weightedSample(candidates$duration_bin, candidates$probability, noise)
  parseDurationBin(chosenBin, noise)
}


#' Sample a transition from current regimen
#' @return Named list with transitionType and toRegimenId
#' @keywords internal
sampleTransition <- function(pack, cohortId, profileId, fromRegimenId, noise) {
  cid <- cohortId; pid <- profileId; frid <- fromRegimenId
  candidates <- pack$regimenTransition[
    pack$regimenTransition$cohort_id == cid &
    pack$regimenTransition$profile_id == pid &
    pack$regimenTransition$from_regimen_id == frid,
  ]
  if (nrow(candidates) == 0) {
    return(list(transitionType = "discontinue", toRegimenId = ""))
  }
  idx <- as.integer(weightedSample(
    seq_len(nrow(candidates)), candidates$probability, noise
  ))
  list(
    transitionType = candidates$transition_type[idx],
    toRegimenId = candidates$to_regimen_id[idx]
  )
}


#' Generate therapy plans and drug exposures for all patients
#'
#' @param patients data.table from generatePatients
#' @param pack parameterPack
#' @param config simulatorConfig
#' @return A list with elements: linePlans, drugExposures, groundTruthEpisodes
#' @export
generateTherapy <- function(patients, pack, config) {
  set.seed(config$seed + 1L)

  allLines <- vector("list", nrow(patients))
  allExposures <- vector("list", nrow(patients))
  allEpisodes <- vector("list", nrow(patients))

  for (i in seq_len(nrow(patients))) {
    pat <- patients[i]

    lines <- generateLinePlansForPatient(pat, pack, config)
    exposures <- generateDrugExposuresForPatient(pat, lines, pack, config)
    episodes <- buildGroundTruthEpisodes(pat, lines, pack)

    allLines[[i]] <- lines
    allExposures[[i]] <- exposures
    allEpisodes[[i]] <- episodes
  }

  list(
    linePlans = rbindlist(allLines),
    drugExposures = rbindlist(allExposures),
    groundTruthEpisodes = rbindlist(allEpisodes)
  )
}


#' Generate line plans for one patient
#' @keywords internal
generateLinePlansForPatient <- function(pat, pack, config) {
  lines <- list()
  currentDay <- 0L
  maxLines <- 4L

  for (lineNum in seq_len(maxLines)) {
    regimenId <- sampleRegimenForLine(
      pack, pat$cohort_id, pat$profile_id, lineNum, config$noise
    )
    if (is.null(regimenId)) break

    duration <- sampleDuration(
      pack, pat$cohort_id, pat$profile_id, regimenId, config$noise
    )

    lineStart <- currentDay
    lineEnd <- currentDay + duration

    transition <- sampleTransition(
      pack, pat$cohort_id, pat$profile_id, regimenId, config$noise
    )

    lines[[lineNum]] <- data.table(
      patient_id = pat$patient_id,
      line_number = lineNum,
      regimen_id = regimenId,
      line_start_day = lineStart,
      line_end_day = lineEnd,
      transition_type = transition$transitionType
    )

    if (transition$transitionType == "discontinue") break

    gap <- sample(7L:28L, 1)
    currentDay <- lineEnd + gap
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
getActiveDrugs <- function(pack, regimenId, pat) {
  rid <- regimenId
  baseDrugs <- pack$regimenDrugMap[pack$regimenDrugMap$regimen_id == rid, ]
  activeIds <- baseDrugs$drug_feature_id

  cid <- pat$cohort_id; pid <- pat$profile_id
  mods <- pack$drugDropAddRules[
    pack$drugDropAddRules$regimen_id == rid &
    pack$drugDropAddRules$cohort_id == cid &
    pack$drugDropAddRules$profile_id == pid,
  ]

  if (nrow(mods) > 0) {
    for (j in seq_len(nrow(mods))) {
      mod <- mods[j]
      if (runif(1) < mod$probability) {
        if (mod$modification_type == "drop") {
          activeIds <- setdiff(activeIds, mod$drug_feature_id)
        } else if (mod$modification_type == "add") {
          activeIds <- union(activeIds, mod$drug_feature_id)
        } else if (mod$modification_type == "substitute") {
          nonAnchors <- baseDrugs[
            drug_role != "anchor" & drug_feature_id %in% activeIds,
            drug_feature_id
          ]
          if (length(nonAnchors) > 0) {
            activeIds <- setdiff(activeIds, nonAnchors[1])
          }
          activeIds <- union(activeIds, mod$drug_feature_id)
        }
      }
    }
  }

  baseDrugs[baseDrugs$drug_feature_id %in% activeIds, ]
}


#' Generate drug exposures for one patient from their line plans
#' @param pat Single-row data.table for one patient
#' @param linePlans data.table of line plans for this patient
#' @param pack parameterPack
#' @param config simulatorConfig
#' @return data.table of drug exposures
#' @export
generateDrugExposures <- function(pat, linePlans, pack, config) {
  generateDrugExposuresForPatient(pat, linePlans, pack, config)
}

#' @keywords internal
generateDrugExposuresForPatient <- function(pat, linePlans, pack, config) {
  if (nrow(linePlans) == 0) {
    return(data.table(
      patient_id = integer(0), drug_feature_id = character(0),
      drug_exposure_start_day = integer(0), drug_exposure_end_day = integer(0),
      regimen_id = character(0), line_number = integer(0),
      is_supportive = logical(0)
    ))
  }

  # Pre-build lookups by regimen_id
  schedByRegimen <- split(pack$drugSchedule, pack$drugSchedule$regimen_id)
  suppByRegimen <- split(pack$supportiveDrugRules, pack$supportiveDrugRules$regimen_id)

  exposures <- list()
  expIdx <- 0L

  for (li in seq_len(nrow(linePlans))) {
    line <- linePlans[li]
    activeDrugs <- getActiveDrugs(pack, line$regimen_id, pat)

    rid <- line$regimen_id
    schedules <- schedByRegimen[[rid]]
    if (is.null(schedules)) schedules <- pack$drugSchedule[0, ]
    schedMap <- setNames(
      split(schedules, seq_len(nrow(schedules))),
      schedules$drug_feature_id
    )

    lineDuration <- line$line_end_day - line$line_start_day
    if (lineDuration <= 0L) next

    for (di in seq_len(nrow(activeDrugs))) {
      drugId <- activeDrugs$drug_feature_id[di]
      sched <- schedMap[[drugId]]
      if (is.null(sched) || nrow(sched) == 0) next
      sched <- sched[1, ]

      if (sched$exposure_type == "continuous") {
        if (!shouldDropExposure(config$noise)) {
          expIdx <- expIdx + 1L
          exposures[[expIdx]] <- data.table(
            patient_id = pat$patient_id,
            drug_feature_id = drugId,
            drug_exposure_start_day = line$line_start_day,
            drug_exposure_end_day = line$line_end_day,
            regimen_id = line$regimen_id,
            line_number = line$line_number,
            is_supportive = FALSE
          )
        }
      } else {
        adminDays <- parseAdminDayPattern(sched$admin_day_pattern)
        cycleLen <- sched$cycle_length_days

        cycleStart <- 0L
        while (cycleStart < lineDuration) {
          for (adminDay in adminDays) {
            actualDay <- cycleStart + adminDay
            actualDay <- jitterCycleDay(actualDay, config$noise)

            if (actualDay >= lineDuration) next

            absStart <- line$line_start_day + actualDay
            expLen <- if (sched$exposure_type == "infusion") 2L else 1L

            if (shouldDropExposure(config$noise)) next

            expIdx <- expIdx + 1L
            exposures[[expIdx]] <- data.table(
              patient_id = pat$patient_id,
              drug_feature_id = drugId,
              drug_exposure_start_day = absStart,
              drug_exposure_end_day = absStart + expLen,
              regimen_id = line$regimen_id,
              line_number = line$line_number,
              is_supportive = FALSE
            )
          }
          cycleStart <- cycleStart + cycleLen
        }
      }
    }

    suppRules <- suppByRegimen[[rid]]
    if (is.null(suppRules)) suppRules <- pack$supportiveDrugRules[0, ]
    if (nrow(suppRules) > 0) {
      for (si in seq_len(nrow(suppRules))) {
        rule <- suppRules[si]
        if (runif(1) < rule$probability) {
          nAdmin <- max(1L, sample.int(max(2L, lineDuration %/% 14L + 1L), 1))
          for (ai in seq_len(nAdmin)) {
            dayOffset <- sample.int(max(1L, lineDuration), 1) - 1L
            absDay <- line$line_start_day + dayOffset

            if (shouldDropExposure(config$noise)) next

            expIdx <- expIdx + 1L
            exposures[[expIdx]] <- data.table(
              patient_id = pat$patient_id,
              drug_feature_id = rule$drug_feature_id,
              drug_exposure_start_day = absDay,
              drug_exposure_end_day = absDay + 1L,
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
buildGroundTruthEpisodes <- function(pat, linePlans, pack) {
  if (nrow(linePlans) == 0) {
    return(data.table(
      patient_id = integer(0), episode_id = integer(0),
      regimen_id = character(0), regimen_name = character(0),
      episode_start_day = integer(0), episode_end_day = integer(0),
      line_number = integer(0)
    ))
  }

  regimenMap <- setNames(pack$regimenCatalog$regimen_name,
                         pack$regimenCatalog$regimen_id)

  data.table(
    patient_id = pat$patient_id,
    episode_id = seq_len(nrow(linePlans)),
    regimen_id = linePlans$regimen_id,
    regimen_name = ifelse(
      linePlans$regimen_id %in% names(regimenMap),
      regimenMap[linePlans$regimen_id],
      linePlans$regimen_id
    ),
    episode_start_day = linePlans$line_start_day,
    episode_end_day = linePlans$line_end_day,
    line_number = linePlans$line_number
  )
}
