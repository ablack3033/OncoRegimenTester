test_that("generatePatients creates correct number of patients", {
  paramDir <- system.file("example_params", "mCRC", package = "OncoRegimenTester")
  if (paramDir == "") {
    paramDir <- file.path("inst", "example_params", "mCRC")
  }
  if (!dir.exists(paramDir)) skip("Example params not found")

  pack <- loadParameterPack(paramDir)
  config <- simulatorConfig(paramDir = paramDir, nPatients = 20L, seed = 123L)
  patients <- generatePatients(pack, config)

  expect_equal(nrow(patients), 20)
  expect_true("patient_id" %in% names(patients))
  expect_true("profile_id" %in% names(patients))
  expect_true("cohort_id" %in% names(patients))
  expect_true(all(patients$patient_id == 1:20))
})


test_that("generatePatients is deterministic with same seed", {
  paramDir <- system.file("example_params", "mCRC", package = "OncoRegimenTester")
  if (paramDir == "") {
    paramDir <- file.path("inst", "example_params", "mCRC")
  }
  if (!dir.exists(paramDir)) skip("Example params not found")

  pack <- loadParameterPack(paramDir)
  config <- simulatorConfig(paramDir = paramDir, nPatients = 10L, seed = 99L)

  p1 <- generatePatients(pack, config)
  p2 <- generatePatients(pack, config)

  expect_identical(p1$profile_id, p2$profile_id)
  expect_identical(p1$age_group, p2$age_group)
  expect_identical(p1$sex, p2$sex)
})
