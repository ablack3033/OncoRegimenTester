#' Minimal reference regimen detector for self-validation
#'
#' @description
#' A SIMPLE, OPTIONAL reference detector so the simulator can sanity-check
#' its own output. It is NOT the production regimen detection algorithm --
#' that belongs in a separate package.
#'
#' Groups nearby drug exposures, ignores supportive drugs, matches drug
#' combinations against known regimen templates.

#' @import data.table

DEFAULT_OVERLAP_WINDOW <- 28L


#' Build regimen signatures from the parameter pack
#'
#' Each regimen's required antineoplastic drugs form its signature.
#' @keywords internal
buildRegimenSignatures <- function(pack) {
  drugCat <- setNames(
    pack$drugFeatureCatalog$is_antineoplastic,
    pack$drugFeatureCatalog$drug_feature_id
  )

  sigs <- list()
  for (i in seq_len(nrow(pack$regimenCatalog))) {
    regId <- pack$regimenCatalog$regimen_id[i]
    drugs <- pack$regimenDrugMap[pack$regimenDrugMap$regimen_id == regId, ]
    sigDrugs <- drugs[
      drugs$required_flag == TRUE &
      drugs$drug_feature_id %in% names(drugCat) &
      drugCat[drugs$drug_feature_id] == TRUE,
      "drug_feature_id"
    ][[1]]
    if (length(sigDrugs) > 0) {
      sigs[[regId]] <- sort(sigDrugs)
    }
  }
  sigs
}


#' Get drug IDs that should be ignored for detection
#' @keywords internal
getIgnorableDrugs <- function(pack) {
  ignorable <- character(0)
  if (nrow(pack$supportiveDrugRules) > 0) {
    ignorable <- union(
      ignorable,
      pack$supportiveDrugRules[
        pack$supportiveDrugRules$ignore_for_regimen_detection_flag == TRUE,
        drug_feature_id
      ]
    )
  }
  ignorable <- union(
    ignorable,
    pack$drugFeatureCatalog[
      pack$drugFeatureCatalog$is_supportive == TRUE,
      drug_feature_id
    ]
  )
  ignorable
}


#' Detect regimen episodes from raw drug exposures
#'
#' @param exposures data.table of drug exposures
#' @param pack parameterPack
#' @param overlapWindow Days gap that starts a new episode
#' @return A list with `episodes` and `lines` data.tables
#' @export
detectRegimens <- function(exposures, pack,
                           overlapWindow = DEFAULT_OVERLAP_WINDOW) {
  emptyEpisodes <- data.table(
    patient_id = integer(0), episode_id = integer(0),
    episode_start_day = integer(0), episode_end_day = integer(0),
    detected_regimen_id = character(0), detected_regimen_name = character(0),
    line_number = integer(0), drug_set_signature = character(0)
  )
  emptyLines <- data.table(
    patient_id = integer(0), line_number = integer(0),
    line_start_day = integer(0), line_end_day = integer(0),
    regimen_sequence_summary = character(0)
  )

  if (nrow(exposures) == 0) {
    return(list(episodes = emptyEpisodes, lines = emptyLines))
  }

  ignorable <- getIgnorableDrugs(pack)
  regimenSigs <- buildRegimenSignatures(pack)
  regimenNames <- setNames(
    pack$regimenCatalog$regimen_name,
    pack$regimenCatalog$regimen_id
  )

  patientIds <- unique(exposures$patient_id)
  allEpisodes <- list()
  allLines <- list()

  for (pid in patientIds) {
    patExp <- exposures[
      exposures$patient_id == pid &
      !(exposures$drug_feature_id %in% ignorable),
    ]
    patExp <- patExp[order(patExp$drug_exposure_start_day), ]

    if (nrow(patExp) == 0) next

    # Group into windows by gaps
    windows <- list()
    currentWindow <- list(patExp[1, ])
    lastEnd <- patExp$drug_exposure_end_day[1]
    for (k in seq_len(nrow(patExp))[-1]) {
      if (patExp$drug_exposure_start_day[k] - lastEnd > overlapWindow) {
        windows[[length(windows) + 1L]] <- rbindlist(currentWindow)
        currentWindow <- list(patExp[k, ])
        lastEnd <- patExp$drug_exposure_end_day[k]
      } else {
        currentWindow[[length(currentWindow) + 1L]] <- patExp[k, ]
        lastEnd <- max(lastEnd, patExp$drug_exposure_end_day[k])
      }
    }
    windows[[length(windows) + 1L]] <- rbindlist(currentWindow)

    episodeId <- 0L
    lineNumber <- 1L
    prevDrugSet <- character(0)
    epList <- list()

    for (wi in seq_along(windows)) {
      w <- windows[[wi]]
      drugSet <- sort(unique(w$drug_feature_id))
      startDay <- min(w$drug_exposure_start_day)
      endDay <- max(w$drug_exposure_end_day)

      detectedId <- "unknown"
      detectedName <- "Unknown"
      for (regId in names(regimenSigs)) {
        sig <- regimenSigs[[regId]]
        if (all(sig %in% drugSet)) {
          detectedId <- regId
          detectedName <- ifelse(
            regId %in% names(regimenNames),
            regimenNames[regId], regId
          )
          break
        }
      }

      if (wi > 1 && length(prevDrugSet) > 0) {
        overlapCount <- length(intersect(drugSet, prevDrugSet))
        if (overlapCount < length(prevDrugSet) * 0.5) {
          lineNumber <- lineNumber + 1L
        }
      }
      prevDrugSet <- drugSet

      episodeId <- episodeId + 1L
      epList[[episodeId]] <- data.table(
        patient_id = pid,
        episode_id = episodeId,
        episode_start_day = startDay,
        episode_end_day = endDay,
        detected_regimen_id = detectedId,
        detected_regimen_name = detectedName,
        line_number = lineNumber,
        drug_set_signature = paste(drugSet, collapse = "+")
      )
    }

    if (length(epList) > 0) {
      episodesDt <- rbindlist(epList)
      allEpisodes[[length(allEpisodes) + 1L]] <- episodesDt

      for (ln in unique(episodesDt$line_number)) {
        lnEps <- episodesDt[episodesDt$line_number == ln, ]
        allLines[[length(allLines) + 1L]] <- data.table(
          patient_id = pid,
          line_number = ln,
          line_start_day = min(lnEps$episode_start_day),
          line_end_day = max(lnEps$episode_end_day),
          regimen_sequence_summary = paste(lnEps$detected_regimen_name, collapse = " -> ")
        )
      }
    }
  }

  list(
    episodes = if (length(allEpisodes) > 0) rbindlist(allEpisodes) else emptyEpisodes,
    lines = if (length(allLines) > 0) rbindlist(allLines) else emptyLines
  )
}
