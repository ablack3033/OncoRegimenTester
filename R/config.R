#' Configuration for the OncoRegimenTester simulator
#'
#' @description
#' Creates configuration objects that control simulation behavior.
#' All settings are explicit and interpretable.

#' Create a noise configuration
#'
#' Controls generation-time stochastic jitter for realism.
#' This is NOT differential privacy -- it is simulation variability
#' so that outputs are not perfectly deterministic copies of the parameter tables.
#'
#' @param probabilityJitterSd SD on logit scale applied to sampling probabilities
#' @param durationJitterDays Max +/- days added to sampled regimen durations
#' @param cycleJitterDays Max +/- days shift on admin days within a cycle
#' @param transitionJitterSd SD on logit scale for transition probabilities
#' @param missingnessRate Probability of dropping an individual drug exposure
#' @return A list with class "noiseConfig"
#' @export
noiseConfig <- function(probabilityJitterSd = 0,
                        durationJitterDays = 0L,
                        cycleJitterDays = 0L,
                        transitionJitterSd = 0,
                        missingnessRate = 0) {
  cfg <- list(
    probabilityJitterSd = probabilityJitterSd,
    durationJitterDays = as.integer(durationJitterDays),
    cycleJitterDays = as.integer(cycleJitterDays),
    transitionJitterSd = transitionJitterSd,
    missingnessRate = missingnessRate
  )
  class(cfg) <- "noiseConfig"
  cfg
}


#' Create a simulator configuration
#'
#' @param paramDir Directory containing the CSV parameter pack
#' @param outputDir Directory for generated outputs
#' @param seed Master random seed for reproducibility
#' @param nPatients Number of synthetic patients to generate
#' @param calendarStart Earliest possible index date (ISO format)
#' @param calendarEnd Latest possible index date (ISO format)
#' @param noise A noiseConfig object
#' @param generateConditions Whether to generate condition records
#' @param generateProcedures Whether to generate procedure records
#' @param generateMeasurements Whether to generate measurement records
#' @param omopOutput Whether to render OMOP-like output tables
#' @param cohortId If set, restrict to this cohort only
#' @return A list with class "simulatorConfig"
#' @export
simulatorConfig <- function(paramDir = "inst/example_params/mCRC",
                            outputDir = "output",
                            seed = 42L,
                            nPatients = 100L,
                            calendarStart = "2018-01-01",
                            calendarEnd = "2023-12-31",
                            noise = noiseConfig(),
                            generateConditions = TRUE,
                            generateProcedures = TRUE,
                            generateMeasurements = TRUE,
                            omopOutput = TRUE,
                            cohortId = NULL) {
  cfg <- list(
    paramDir = paramDir,
    outputDir = outputDir,
    seed = as.integer(seed),
    nPatients = as.integer(nPatients),
    calendarStart = as.Date(calendarStart),
    calendarEnd = as.Date(calendarEnd),
    noise = noise,
    generateConditions = generateConditions,
    generateProcedures = generateProcedures,
    generateMeasurements = generateMeasurements,
    omopOutput = omopOutput,
    cohortId = cohortId
  )
  class(cfg) <- "simulatorConfig"
  cfg
}
