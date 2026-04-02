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

# Default gap (days) before starting a new episode
DEFAULT_OVERLAP_WINDOW <- 28L


#' Build regimen signatures from the parameter pack
#'
#' Each regimen's required antineoplastic drugs form its signature.
#' @keywords internal
build_regimen_signatures <- function(pack) {
  drug_cat <- setNames(
    pack$drug_feature_catalog$is_antineoplastic,
    pack$drug_feature_catalog$drug_feature_id
  )

  sigs <- list()
  for (i in seq_len(nrow(pack$regimen_catalog))) {
    reg_id <- pack$regimen_catalog$regimen_id[i]
    drugs <- pack$regimen_drug_map[pack$regimen_drug_map$regimen_id == reg_id, ]
    sig_drugs <- drugs[
      drugs$required_flag == TRUE &
      drugs$drug_feature_id %in% names(drug_cat) &
      drug_cat[drugs$drug_feature_id] == TRUE,
      "drug_feature_id"
    ][[1]]
    if (length(sig_drugs) > 0) {
      sigs[[reg_id]] <- sort(sig_drugs)
    }
  }
  sigs
}


#' Get drug IDs that should be ignored for detection
#' @keywords internal
get_ignorable_drugs <- function(pack) {
  ignorable <- character(0)
  if (nrow(pack$supportive_drug_rules) > 0) {
    ignorable <- union(
      ignorable,
      pack$supportive_drug_rules[
        pack$supportive_drug_rules$ignore_for_regimen_detection_flag == TRUE,
        drug_feature_id
      ]
    )
  }
  ignorable <- union(
    ignorable,
    pack$drug_feature_catalog[
      pack$drug_feature_catalog$is_supportive == TRUE,
      drug_feature_id
    ]
  )
  ignorable
}


#' Detect regimen episodes from raw drug exposures
#'
#' @param exposures data.table of drug exposures
#' @param pack parameter_pack
#' @param overlap_window Days gap that starts a new episode
#' @return A list with `episodes` and `lines` data.tables
#' @export
detect_regimens <- function(exposures, pack,
                            overlap_window = DEFAULT_OVERLAP_WINDOW) {
  empty_episodes <- data.table(
    patient_id = integer(0), episode_id = integer(0),
    episode_start_day = integer(0), episode_end_day = integer(0),
    detected_regimen_id = character(0), detected_regimen_name = character(0),
    line_number = integer(0), drug_set_signature = character(0)
  )
  empty_lines <- data.table(
    patient_id = integer(0), line_number = integer(0),
    line_start_day = integer(0), line_end_day = integer(0),
    regimen_sequence_summary = character(0)
  )

  if (nrow(exposures) == 0) {
    return(list(episodes = empty_episodes, lines = empty_lines))
  }

  ignorable <- get_ignorable_drugs(pack)
  regimen_sigs <- build_regimen_signatures(pack)
  regimen_names <- setNames(
    pack$regimen_catalog$regimen_name,
    pack$regimen_catalog$regimen_id
  )

  # Process per patient
  patient_ids <- unique(exposures$patient_id)
  all_episodes <- list()
  all_lines <- list()

  for (pid in patient_ids) {
    pat_exp <- exposures[
      exposures$patient_id == pid &
      !(exposures$drug_feature_id %in% ignorable),
    ]
    pat_exp <- pat_exp[order(pat_exp$drug_exposure_start_day), ]

    if (nrow(pat_exp) == 0) next

    # Group into windows by gaps
    windows <- list()
    current_window <- list(pat_exp[1, ])
    for (k in seq_len(nrow(pat_exp))[-1]) {
      last_end <- max(sapply(current_window, function(x) x$drug_exposure_end_day))
      if (pat_exp$drug_exposure_start_day[k] - last_end > overlap_window) {
        windows[[length(windows) + 1L]] <- rbindlist(current_window)
        current_window <- list(pat_exp[k, ])
      } else {
        current_window[[length(current_window) + 1L]] <- pat_exp[k, ]
      }
    }
    windows[[length(windows) + 1L]] <- rbindlist(current_window)

    episode_id <- 0L
    line_number <- 1L
    prev_drug_set <- character(0)
    ep_list <- list()

    for (wi in seq_along(windows)) {
      w <- windows[[wi]]
      drug_set <- sort(unique(w$drug_feature_id))
      start_day <- min(w$drug_exposure_start_day)
      end_day <- max(w$drug_exposure_end_day)

      # Match against regimen signatures
      detected_id <- "unknown"
      detected_name <- "Unknown"
      for (reg_id in names(regimen_sigs)) {
        sig <- regimen_sigs[[reg_id]]
        if (all(sig %in% drug_set)) {
          detected_id <- reg_id
          detected_name <- ifelse(
            reg_id %in% names(regimen_names),
            regimen_names[reg_id], reg_id
          )
          break
        }
      }

      # Line advancement: if drug set changes significantly
      if (wi > 1 && length(prev_drug_set) > 0) {
        overlap_count <- length(intersect(drug_set, prev_drug_set))
        if (overlap_count < length(prev_drug_set) * 0.5) {
          line_number <- line_number + 1L
        }
      }
      prev_drug_set <- drug_set

      episode_id <- episode_id + 1L
      ep_list[[episode_id]] <- data.table(
        patient_id = pid,
        episode_id = episode_id,
        episode_start_day = start_day,
        episode_end_day = end_day,
        detected_regimen_id = detected_id,
        detected_regimen_name = detected_name,
        line_number = line_number,
        drug_set_signature = paste(drug_set, collapse = "+")
      )
    }

    if (length(ep_list) > 0) {
      episodes_dt <- rbindlist(ep_list)
      all_episodes[[length(all_episodes) + 1L]] <- episodes_dt

      # Derive lines
      for (ln in unique(episodes_dt$line_number)) {
        ln_eps <- episodes_dt[episodes_dt$line_number == ln, ]
        all_lines[[length(all_lines) + 1L]] <- data.table(
          patient_id = pid,
          line_number = ln,
          line_start_day = min(ln_eps$episode_start_day),
          line_end_day = max(ln_eps$episode_end_day),
          regimen_sequence_summary = paste(ln_eps$detected_regimen_name, collapse = " -> ")
        )
      }
    }
  }

  list(
    episodes = if (length(all_episodes) > 0) rbindlist(all_episodes) else empty_episodes,
    lines = if (length(all_lines) > 0) rbindlist(all_lines) else empty_lines
  )
}
