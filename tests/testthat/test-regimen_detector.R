test_that("detect_regimens finds episodes from generated exposures", {
  param_dir <- system.file("example_params", "mCRC", package = "OncoRegimenTester")
  if (param_dir == "") {
    param_dir <- file.path("inst", "example_params", "mCRC")
  }
  if (!dir.exists(param_dir)) skip("Example params not found")

  pack <- load_parameter_pack(param_dir)
  config <- simulator_config(param_dir = param_dir, n_patients = 10L, seed = 42L)
  patients <- generate_patients(pack, config)
  therapy <- generate_therapy(patients, pack, config)

  detected <- detect_regimens(therapy$drug_exposures, pack)

  expect_true(nrow(detected$episodes) > 0)
  expect_true(nrow(detected$lines) > 0)

  # Every detected episode should have a patient_id from our patient set
  expect_true(all(detected$episodes$patient_id %in% patients$patient_id))
  expect_true(all(detected$lines$patient_id %in% patients$patient_id))
})


test_that("detect_regimens handles empty input", {
  param_dir <- system.file("example_params", "mCRC", package = "OncoRegimenTester")
  if (param_dir == "") {
    param_dir <- file.path("inst", "example_params", "mCRC")
  }
  if (!dir.exists(param_dir)) skip("Example params not found")

  pack <- load_parameter_pack(param_dir)
  empty_exp <- data.table::data.table(
    patient_id = integer(0), drug_feature_id = character(0),
    drug_exposure_start_day = integer(0), drug_exposure_end_day = integer(0),
    regimen_id = character(0), line_number = integer(0),
    is_supportive = logical(0)
  )

  detected <- detect_regimens(empty_exp, pack)
  expect_equal(nrow(detected$episodes), 0)
  expect_equal(nrow(detected$lines), 0)
})
