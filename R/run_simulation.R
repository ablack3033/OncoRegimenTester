#' Run a complete simulation
#'
#' @description
#' Top-level function that orchestrates the full simulation pipeline:
#' load parameters -> generate patients -> generate therapy -> detect regimens
#' -> render OMOP -> validate.

#' @import data.table


#' Run the full simulation pipeline
#'
#' @param config A simulatorConfig object. If NULL, uses defaults.
#' @return A list with all generated data and validation report
#' @export
runSimulation <- function(config = NULL) {
  if (is.null(config)) {
    config <- simulatorConfig()
  }

  message("Loading parameter pack from: ", config$paramDir)
  pack <- loadParameterPack(config$paramDir, cohortId = config$cohortId)

  if (length(pack$warnings) > 0) {
    message(sprintf("Parameter pack loaded with %d warnings", length(pack$warnings)))
    for (w in pack$warnings) {
      message("  WARNING: ", w)
    }
  } else {
    message("Parameter pack loaded successfully")
  }

  message(sprintf("Generating %d patients...", config$nPatients))
  patients <- generatePatients(pack, config)

  message("Generating therapy (lines, regimens, drug exposures)...")
  therapy <- generateTherapy(patients, pack, config)

  message(sprintf("  Generated %d line plans", nrow(therapy$linePlans)))
  message(sprintf("  Generated %d drug exposures", nrow(therapy$drugExposures)))

  message("Generating supporting events...")
  supporting <- generateSupportingEvents(
    patients, therapy$linePlans, pack, config
  )
  message(sprintf("  Conditions: %d, Procedures: %d, Measurements: %d",
                  nrow(supporting$conditions),
                  nrow(supporting$procedures),
                  nrow(supporting$measurements)))

  message("Running reference regimen detector...")
  detected <- detectRegimens(therapy$drugExposures, pack)
  message(sprintf("  Detected %d episodes, %d lines",
                  nrow(detected$episodes), nrow(detected$lines)))

  if (config$omopOutput) {
    message("Rendering OMOP-like output tables...")
    renderOmop(patients, therapy, supporting, detected, config)
  }

  message("Running validation...")
  validation <- validateOutput(patients, therapy, pack)

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
