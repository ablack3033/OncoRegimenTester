test_that("generateTherapy produces line plans and drug exposures", {
  paramDir <- system.file("example_params", "mCRC", package = "OncoRegimenTester")
  if (paramDir == "") {
    paramDir <- file.path("inst", "example_params", "mCRC")
  }
  if (!dir.exists(paramDir)) skip("Example params not found")

  pack <- loadParameterPack(paramDir)
  config <- simulatorConfig(paramDir = paramDir, nPatients = 10L, seed = 42L)
  patients <- generatePatients(pack, config)

  therapy <- generateTherapy(patients, pack, config)

  expect_true(nrow(therapy$linePlans) > 0)
  expect_true(nrow(therapy$drugExposures) > 0)
  expect_true(nrow(therapy$groundTruthEpisodes) > 0)

  # Every patient should have at least one line plan
  expect_true(all(patients$patient_id %in% therapy$linePlans$patient_id))

  # Drug exposures should reference valid patient IDs
  expect_true(all(therapy$drugExposures$patient_id %in% patients$patient_id))

  # Line numbers should be positive
  expect_true(all(therapy$linePlans$line_number >= 1))
})


test_that("parseAdminDayPattern works correctly", {
  result <- OncoRegimenTester:::parseAdminDayPattern("1")
  expect_equal(result, 0L)

  result <- OncoRegimenTester:::parseAdminDayPattern("1, 8")
  expect_equal(result, c(0L, 7L))

  result <- OncoRegimenTester:::parseAdminDayPattern("1,8,15")
  expect_equal(result, c(0L, 7L, 14L))
})
