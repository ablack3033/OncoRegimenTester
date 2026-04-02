test_that("end-to-end simulation runs without error", {
  param_dir <- system.file("example_params", "mCRC", package = "OncoRegimenTester")
  if (param_dir == "") {
    param_dir <- file.path("inst", "example_params", "mCRC")
  }
  if (!dir.exists(param_dir)) skip("Example params not found")

  output_dir <- tempdir()

  config <- simulator_config(
    param_dir = param_dir,
    output_dir = file.path(output_dir, "test_output"),
    n_patients = 20L,
    seed = 42L,
    omop_output = TRUE,
    noise = noise_config(
      probability_jitter_sd = 0.1,
      duration_jitter_days = 3L,
      cycle_jitter_days = 1L,
      missingness_rate = 0.05
    )
  )

  result <- run_simulation(config)

  # Check that all outputs exist
  expect_true(nrow(result$patients) == 20)
  expect_true(nrow(result$therapy$line_plans) > 0)
  expect_true(nrow(result$therapy$drug_exposures) > 0)
  expect_true(nrow(result$therapy$ground_truth_episodes) > 0)
  expect_true(nrow(result$detected$episodes) > 0)

  # Check output files exist
  expect_true(file.exists(file.path(config$output_dir, "person.csv")))
  expect_true(file.exists(file.path(config$output_dir, "drug_exposure.csv")))
  expect_true(file.exists(file.path(config$output_dir, "ground_truth_regimen_episodes.csv")))

  # Validation report should exist
  expect_true(length(result$validation$report) > 0)

  # Clean up
  unlink(file.path(output_dir, "test_output"), recursive = TRUE)
})


test_that("simulation is reproducible with same seed", {
  param_dir <- system.file("example_params", "mCRC", package = "OncoRegimenTester")
  if (param_dir == "") {
    param_dir <- file.path("inst", "example_params", "mCRC")
  }
  if (!dir.exists(param_dir)) skip("Example params not found")

  config <- simulator_config(
    param_dir = param_dir,
    output_dir = tempfile(),
    n_patients = 5L,
    seed = 777L,
    omop_output = FALSE
  )

  r1 <- run_simulation(config)
  r2 <- run_simulation(config)

  expect_identical(r1$patients$profile_id, r2$patients$profile_id)
  expect_identical(
    r1$therapy$drug_exposures$drug_feature_id,
    r2$therapy$drug_exposures$drug_feature_id
  )
})
