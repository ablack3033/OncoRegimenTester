test_that("end-to-end simulation runs without error", {
  paramDir <- system.file("example_params", "mCRC", package = "OncoRegimenTester")
  if (paramDir == "") {
    paramDir <- file.path("inst", "example_params", "mCRC")
  }
  if (!dir.exists(paramDir)) skip("Example params not found")

  outputDir <- tempdir()

  config <- simulatorConfig(
    paramDir = paramDir,
    outputDir = file.path(outputDir, "test_output"),
    nPatients = 20L,
    seed = 42L,
    omopOutput = TRUE,
    noise = noiseConfig(
      probabilityJitterSd = 0.1,
      durationJitterDays = 3L,
      cycleJitterDays = 1L,
      missingnessRate = 0.05
    )
  )

  result <- runSimulation(config)

  expect_true(nrow(result$patients) == 20)
  expect_true(nrow(result$therapy$linePlans) > 0)
  expect_true(nrow(result$therapy$drugExposures) > 0)
  expect_true(nrow(result$therapy$groundTruthEpisodes) > 0)
  expect_true(nrow(result$detected$episodes) > 0)

  expect_true(file.exists(file.path(config$outputDir, "person.csv")))
  expect_true(file.exists(file.path(config$outputDir, "drug_exposure.csv")))
  expect_true(file.exists(file.path(config$outputDir, "ground_truth_regimen_episodes.csv")))

  expect_true(length(result$validation$report) > 0)

  unlink(file.path(outputDir, "test_output"), recursive = TRUE)
})


test_that("simulation is reproducible with same seed", {
  paramDir <- system.file("example_params", "mCRC", package = "OncoRegimenTester")
  if (paramDir == "") {
    paramDir <- file.path("inst", "example_params", "mCRC")
  }
  if (!dir.exists(paramDir)) skip("Example params not found")

  config <- simulatorConfig(
    paramDir = paramDir,
    outputDir = tempfile(),
    nPatients = 5L,
    seed = 777L,
    omopOutput = FALSE
  )

  r1 <- runSimulation(config)
  r2 <- runSimulation(config)

  expect_identical(r1$patients$profile_id, r2$patients$profile_id)
  expect_identical(
    r1$therapy$drugExposures$drug_feature_id,
    r2$therapy$drugExposures$drug_feature_id
  )
})
