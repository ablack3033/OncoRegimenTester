#' Validation utilities
#'
#' @description
#' Compares generated output back to parameter targets.
#' Produces a human-readable validation report.

#' @import data.table


#' Validate simulation output against parameter targets
#'
#' @param patients data.table of generated patients
#' @param therapy list from generate_therapy
#' @param pack parameter_pack
#' @return A list with `report` (character vector) and `checks` (list of check results)
#' @export
validate_output <- function(patients, therapy, pack) {
  report <- character(0)
  checks <- list()

  report <- c(report, "=== OncoRegimenTester Validation Report ===", "")

  # Row counts
  report <- c(report, "--- Row Counts ---")
  report <- c(report,
    sprintf("  Patients: %d", nrow(patients)),
    sprintf("  Line plans: %d", nrow(therapy$line_plans)),
    sprintf("  Drug exposures: %d", nrow(therapy$drug_exposures)),
    sprintf("  Ground truth episodes: %d", nrow(therapy$ground_truth_episodes)),
    ""
  )

  # Profile mix comparison
  report <- c(report, "--- Profile Mix (observed vs target) ---")
  if (nrow(patients) > 0 && nrow(pack$treatment_profile_mix) > 0) {
    obs_profiles <- patients[, .N, by = profile_id]
    obs_profiles[, obs_proportion := N / sum(N)]

    # Get target proportions (averaged across strata)
    target <- pack$treatment_profile_mix[,
      .(target_proportion = mean(proportion)),
      by = profile_id
    ]
    target[, target_proportion := target_proportion / sum(target_proportion)]

    merged <- merge(obs_profiles, target, by = "profile_id", all = TRUE)
    merged[is.na(obs_proportion), obs_proportion := 0]
    merged[is.na(target_proportion), target_proportion := 0]

    for (i in seq_len(nrow(merged))) {
      report <- c(report, sprintf(
        "  %s: observed=%.3f, target=%.3f",
        merged$profile_id[i], merged$obs_proportion[i], merged$target_proportion[i]
      ))
    }
    checks$profile_mix <- merged
  }
  report <- c(report, "")

  # Regimen start frequency
  report <- c(report, "--- Regimen Start Frequency (Line 1) ---")
  if (nrow(therapy$line_plans) > 0) {
    l1 <- therapy$line_plans[therapy$line_plans$line_number == 1, ]
    if (nrow(l1) > 0) {
      obs_reg <- l1[, .N, by = regimen_id]
      obs_reg[, obs_proportion := N / sum(N)]

      for (i in seq_len(nrow(obs_reg))) {
        report <- c(report, sprintf(
          "  %s: observed=%.3f (n=%d)",
          obs_reg$regimen_id[i], obs_reg$obs_proportion[i], obs_reg$N[i]
        ))
      }
      checks$l1_regimens <- obs_reg
    }
  }
  report <- c(report, "")

  # Duration summary
  report <- c(report, "--- Regimen Duration Summary (days) ---")
  if (nrow(therapy$line_plans) > 0) {
    lp <- therapy$line_plans
    lp[, duration := line_end_day - line_start_day]
    dur_summary <- lp[, .(
      mean_dur = mean(duration),
      median_dur = median(duration),
      min_dur = min(duration),
      max_dur = max(duration)
    ), by = regimen_id]

    for (i in seq_len(nrow(dur_summary))) {
      report <- c(report, sprintf(
        "  %s: mean=%.0f, median=%.0f, range=[%d, %d]",
        dur_summary$regimen_id[i], dur_summary$mean_dur[i],
        dur_summary$median_dur[i], dur_summary$min_dur[i], dur_summary$max_dur[i]
      ))
    }
    checks$duration_summary <- dur_summary
  }
  report <- c(report, "")

  # Drug exposure counts per regimen
  report <- c(report, "--- Drug Exposures per Regimen ---")
  if (nrow(therapy$drug_exposures) > 0) {
    de <- therapy$drug_exposures[therapy$drug_exposures$is_supportive == FALSE, ]
    if (nrow(de) > 0) {
      drug_counts <- de[, .N, by = .(regimen_id, drug_feature_id)]
      drug_counts <- drug_counts[order(regimen_id, -N)]
      for (i in seq_len(min(20, nrow(drug_counts)))) {
        report <- c(report, sprintf(
          "  %s / %s: %d exposures",
          drug_counts$regimen_id[i], drug_counts$drug_feature_id[i],
          drug_counts$N[i]
        ))
      }
      checks$drug_counts <- drug_counts
    }
  }
  report <- c(report, "")

  # Parameter pack warnings
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
