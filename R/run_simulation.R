#' Run a complete simulation
#'
#' @description
#' Top-level function that orchestrates the full simulation pipeline:
#' load parameters -> generate patients -> generate therapy -> detect regimens
#' -> render OMOP -> validate.

#' @import data.table


#' Run the full simulation pipeline
#'
#' @param config A simulator_config object. If NULL, uses defaults.
#' @return A list with all generated data and validation report
#' @export
run_simulation <- function(config = NULL) {
  if (is.null(config)) {
    config <- simulator_config()
  }

  message("Loading parameter pack from: ", config$param_dir)
  pack <- load_parameter_pack(config$param_dir, cohort_id = config$cohort_id)

  if (length(pack$warnings) > 0) {
    message(sprintf("Parameter pack loaded with %d warnings", length(pack$warnings)))
    for (w in pack$warnings) {
      message("  WARNING: ", w)
    }
  } else {
    message("Parameter pack loaded successfully")
  }

  message(sprintf("Generating %d patients...", config$n_patients))
  patients <- generate_patients(pack, config)

  message("Generating therapy (lines, regimens, drug exposures)...")
  therapy <- generate_therapy(patients, pack, config)

  message(sprintf("  Generated %d line plans", nrow(therapy$line_plans)))
  message(sprintf("  Generated %d drug exposures", nrow(therapy$drug_exposures)))

  message("Generating supporting events...")
  supporting <- generate_supporting_events(
    patients, therapy$line_plans, pack, config
  )
  message(sprintf("  Conditions: %d, Procedures: %d, Measurements: %d",
                  nrow(supporting$conditions),
                  nrow(supporting$procedures),
                  nrow(supporting$measurements)))

  message("Running reference regimen detector...")
  detected <- detect_regimens(therapy$drug_exposures, pack)
  message(sprintf("  Detected %d episodes, %d lines",
                  nrow(detected$episodes), nrow(detected$lines)))

  if (config$omop_output) {
    message("Rendering OMOP-like output tables...")
    render_omop(patients, therapy, supporting, detected, config)
  }

  message("Running validation...")
  validation <- validate_output(patients, therapy, pack)

  message("")
  cat(paste(validation$report, collapse = "\n"), "\n")

  invisible(list(
    patients = patients,
    therapy = therapy,
    supporting = supporting,
    detected = detected,
    validation = validation,
    pack = pack,
    config = config
  ))
}
