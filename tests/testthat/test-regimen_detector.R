test_that("detectRegimens finds episodes from generated exposures", {
  paramDir <- system.file("example_params", "mCRC", package = "OncoRegimenTester")
  if (paramDir == "") {
    paramDir <- file.path("inst", "example_params", "mCRC")
  }
  if (!dir.exists(paramDir)) skip("Example params not found")

  pack <- loadParameterPack(paramDir)
  config <- simulatorConfig(paramDir = paramDir, nPatients = 10L, seed = 42L)
  patients <- generatePatients(pack, config)
  therapy <- generateTherapy(patients, pack, config)

  detected <- detectRegimens(therapy$drugExposures, pack)

  expect_true(nrow(detected$episodes) > 0)
  expect_true(nrow(detected$lines) > 0)

  expect_true(all(detected$episodes$patient_id %in% patients$patient_id))
  expect_true(all(detected$lines$patient_id %in% patients$patient_id))
})


test_that("detectRegimens handles empty input", {
  paramDir <- system.file("example_params", "mCRC", package = "OncoRegimenTester")
  if (paramDir == "") {
    paramDir <- file.path("inst", "example_params", "mCRC")
  }
  if (!dir.exists(paramDir)) skip("Example params not found")

  pack <- loadParameterPack(paramDir)
  emptyExp <- data.table::data.table(
    patient_id = integer(0), drug_feature_id = character(0),
    drug_exposure_start_day = integer(0), drug_exposure_end_day = integer(0),
    regimen_id = character(0), line_number = integer(0),
    is_supportive = logical(0)
  )

  detected <- detectRegimens(emptyExp, pack)
  expect_equal(nrow(detected$episodes), 0)
  expect_equal(nrow(detected$lines), 0)
})
