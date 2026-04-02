test_that("generate_therapy produces line plans and drug exposures", {
  param_dir <- system.file("example_params", "mCRC", package = "OncoRegimenTester")
  if (param_dir == "") {
    param_dir <- file.path("inst", "example_params", "mCRC")
  }
  if (!dir.exists(param_dir)) skip("Example params not found")

  pack <- load_parameter_pack(param_dir)
  config <- simulator_config(param_dir = param_dir, n_patients = 10L, seed = 42L)
  patients <- generate_patients(pack, config)

  therapy <- generate_therapy(patients, pack, config)

  expect_true(nrow(therapy$line_plans) > 0)
  expect_true(nrow(therapy$drug_exposures) > 0)
  expect_true(nrow(therapy$ground_truth_episodes) > 0)

  # Every patient should have at least one line plan
  expect_true(all(patients$patient_id %in% therapy$line_plans$patient_id))

  # Drug exposures should reference valid patient IDs
  expect_true(all(therapy$drug_exposures$patient_id %in% patients$patient_id))

  # Line numbers should be positive
  expect_true(all(therapy$line_plans$line_number >= 1))
})


test_that("parse_admin_day_pattern works correctly", {
  result <- OncoRegimenTester:::parse_admin_day_pattern("1")
  expect_equal(result, 0L)

  result <- OncoRegimenTester:::parse_admin_day_pattern("1, 8")
  expect_equal(result, c(0L, 7L))

  result <- OncoRegimenTester:::parse_admin_day_pattern("1,8,15")
  expect_equal(result, c(0L, 7L, 14L))
})
