test_that("generate_patients creates correct number of patients", {
  param_dir <- system.file("example_params", "mCRC", package = "OncoRegimenTester")
  if (param_dir == "") {
    param_dir <- file.path("inst", "example_params", "mCRC")
  }
  if (!dir.exists(param_dir)) skip("Example params not found")

  pack <- load_parameter_pack(param_dir)
  config <- simulator_config(param_dir = param_dir, n_patients = 20L, seed = 123L)
  patients <- generate_patients(pack, config)

  expect_equal(nrow(patients), 20)
  expect_true("patient_id" %in% names(patients))
  expect_true("profile_id" %in% names(patients))
  expect_true("cohort_id" %in% names(patients))
  expect_true(all(patients$patient_id == 1:20))
})


test_that("generate_patients is deterministic with same seed", {
  param_dir <- system.file("example_params", "mCRC", package = "OncoRegimenTester")
  if (param_dir == "") {
    param_dir <- file.path("inst", "example_params", "mCRC")
  }
  if (!dir.exists(param_dir)) skip("Example params not found")

  pack <- load_parameter_pack(param_dir)
  config <- simulator_config(param_dir = param_dir, n_patients = 10L, seed = 99L)

  p1 <- generate_patients(pack, config)
  p2 <- generate_patients(pack, config)

  expect_identical(p1$profile_id, p2$profile_id)
  expect_identical(p1$age_group, p2$age_group)
  expect_identical(p1$sex, p2$sex)
})
