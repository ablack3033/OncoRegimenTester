#' Validation utilities
#'
#' @description
#' Compares generated output back to parameter targets.
#' Produces a human-readable validation report.

#' @import data.table


#' Validate simulation output against parameter targets
#'
#' @param patients data.table of generated patients
#' @param therapy list from generateTherapy
#' @param pack parameterPack
#' @return A list with `report` (character vector) and `checks` (list of check results)
#' @export
validateOutput <- function(patients, therapy, pack) {
  report <- character(0)
  checks <- list()

  report <- c(report, "=== OncoRegimenTester Validation Report ===", "")

  report <- c(report, "--- Row Counts ---")
  report <- c(report,
    sprintf("  Patients: %d", nrow(patients)),
    sprintf("  Line plans: %d", nrow(therapy$linePlans)),
    sprintf("  Drug exposures: %d", nrow(therapy$drugExposures)),
    sprintf("  Ground truth episodes: %d", nrow(therapy$groundTruthEpisodes)),
    ""
  )

  report <- c(report, "--- Profile Mix (observed vs target) ---")
  if (nrow(patients) > 0 && nrow(pack$treatmentProfileMix) > 0) {
    obsProfiles <- patients[, .N, by = profile_id]
    obsProfiles[, obs_proportion := N / sum(N)]

    target <- pack$treatmentProfileMix[,
      .(target_proportion = mean(proportion)),
      by = profile_id
    ]
    target[, target_proportion := target_proportion / sum(target_proportion)]

    merged <- merge(obsProfiles, target, by = "profile_id", all = TRUE)
    merged[is.na(obs_proportion), obs_proportion := 0]
    merged[is.na(target_proportion), target_proportion := 0]

    for (i in seq_len(nrow(merged))) {
      report <- c(report, sprintf(
        "  %s: observed=%.3f, target=%.3f",
        merged$profile_id[i], merged$obs_proportion[i], merged$target_proportion[i]
      ))
    }
    checks$profileMix <- merged
  }
  report <- c(report, "")

  report <- c(report, "--- Regimen Start Frequency (Line 1) ---")
  if (nrow(therapy$linePlans) > 0) {
    l1 <- therapy$linePlans[therapy$linePlans$line_number == 1, ]
    if (nrow(l1) > 0) {
      obsReg <- l1[, .N, by = regimen_id]
      obsReg[, obs_proportion := N / sum(N)]

      for (i in seq_len(nrow(obsReg))) {
        report <- c(report, sprintf(
          "  %s: observed=%.3f (n=%d)",
          obsReg$regimen_id[i], obsReg$obs_proportion[i], obsReg$N[i]
        ))
      }
      checks$l1Regimens <- obsReg
    }
  }
  report <- c(report, "")

  report <- c(report, "--- Regimen Duration Summary (days) ---")
  if (nrow(therapy$linePlans) > 0) {
    lp <- copy(therapy$linePlans)
    lp[, duration := line_end_day - line_start_day]
    durSummary <- lp[, .(
      mean_dur = mean(duration),
      median_dur = median(duration),
      min_dur = min(duration),
      max_dur = max(duration)
    ), by = regimen_id]

    for (i in seq_len(nrow(durSummary))) {
      report <- c(report, sprintf(
        "  %s: mean=%.0f, median=%.0f, range=[%d, %d]",
        durSummary$regimen_id[i], durSummary$mean_dur[i],
        durSummary$median_dur[i], durSummary$min_dur[i], durSummary$max_dur[i]
      ))
    }
    checks$durationSummary <- durSummary
  }
  report <- c(report, "")

  report <- c(report, "--- Drug Exposures per Regimen ---")
  if (nrow(therapy$drugExposures) > 0) {
    de <- therapy$drugExposures[therapy$drugExposures$is_supportive == FALSE, ]
    if (nrow(de) > 0) {
      drugCounts <- de[, .N, by = .(regimen_id, drug_feature_id)]
      drugCounts <- drugCounts[order(regimen_id, -N)]
      for (i in seq_len(min(20, nrow(drugCounts)))) {
        report <- c(report, sprintf(
          "  %s / %s: %d exposures",
          drugCounts$regimen_id[i], drugCounts$drug_feature_id[i],
          drugCounts$N[i]
        ))
      }
      checks$drugCounts <- drugCounts
    }
  }
  report <- c(report, "")

  if (length(pack$warnings) > 0) {
    report <- c(report, "--- Parameter Pack Warnings ---")
    for (w in pack$warnings) {
      report <- c(report, sprintf("  WARNING: %s", w))
    }
    report <- c(report, "")
  }

  report <- c(report, "=== End of Validation Report ===")

  list(report = report, checks = checks)
}
