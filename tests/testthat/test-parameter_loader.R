test_that("load_parameter_pack loads example mCRC data", {
  param_dir <- system.file("example_params", "mCRC", package = "OncoRegimenTester")
  if (param_dir == "") {
    param_dir <- file.path("inst", "example_params", "mCRC")
  }
  if (!dir.exists(param_dir)) skip("Example params not found")

  pack <- load_parameter_pack(param_dir)

  expect_s3_class(pack, "parameter_pack")
  expect_true(nrow(pack$manifest) > 0)
  expect_true(nrow(pack$cohort_population) > 0)
  expect_true(nrow(pack$treatment_profile_mix) > 0)
  expect_true(nrow(pack$drug_feature_catalog) > 0)
  expect_true(nrow(pack$regimen_catalog) > 0)
  expect_true(nrow(pack$regimen_drug_map) > 0)
  expect_true(nrow(pack$profile_regimen_start) > 0)
  expect_true(nrow(pack$regimen_transition) > 0)
  expect_true(nrow(pack$regimen_duration) > 0)
  expect_true(nrow(pack$drug_schedule) > 0)
})


test_that("check_columns catches missing columns", {
  dt <- data.table::data.table(a = 1, b = 2)
  expect_error(
    OncoRegimenTester:::check_columns(dt, c("a", "b", "c"), "test.csv"),
    "missing required columns"
  )
})


test_that("check_probability_groups detects bad sums", {
  dt <- data.table::data.table(
    group = c("A", "A", "B", "B"),
    prob = c(0.3, 0.3, 0.5, 0.5)
  )
  warnings <- OncoRegimenTester:::check_probability_groups(
    dt, "group", "prob", "test.csv"
  )
  # Group A sums to 0.6, should warn

expect_true(length(warnings) > 0)
  expect_true(any(grepl("group.*A", warnings, ignore.case = TRUE)))
})


test_that("parse_bool handles various inputs", {
  expect_true(OncoRegimenTester:::parse_bool("TRUE"))
  expect_true(OncoRegimenTester:::parse_bool("1"))
  expect_true(OncoRegimenTester:::parse_bool("yes"))
  expect_false(OncoRegimenTester:::parse_bool("FALSE"))
  expect_false(OncoRegimenTester:::parse_bool("0"))
  expect_false(OncoRegimenTester:::parse_bool("no"))
})


test_that("validate_parameter_pack returns summary", {
  param_dir <- system.file("example_params", "mCRC", package = "OncoRegimenTester")
  if (param_dir == "") {
    param_dir <- file.path("inst", "example_params", "mCRC")
  }
  if (!dir.exists(param_dir)) skip("Example params not found")

  pack <- load_parameter_pack(param_dir)
  result <- validate_parameter_pack(pack)

  expect_true(is.list(result))
  expect_true("valid" %in% names(result))
  expect_true("summary" %in% names(result))
  expect_true(length(result$summary) > 0)
})
