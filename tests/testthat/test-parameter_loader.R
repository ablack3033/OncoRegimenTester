test_that("loadParameterPack loads example mCRC data", {
  paramDir <- system.file("example_params", "mCRC", package = "OncoRegimenTester")
  if (paramDir == "") {
    paramDir <- file.path("inst", "example_params", "mCRC")
  }
  if (!dir.exists(paramDir)) skip("Example params not found")

  pack <- loadParameterPack(paramDir)

  expect_s3_class(pack, "parameterPack")
  expect_true(nrow(pack$manifest) > 0)
  expect_true(nrow(pack$cohortPopulation) > 0)
  expect_true(nrow(pack$treatmentProfileMix) > 0)
  expect_true(nrow(pack$drugFeatureCatalog) > 0)
  expect_true(nrow(pack$regimenCatalog) > 0)
  expect_true(nrow(pack$regimenDrugMap) > 0)
  expect_true(nrow(pack$profileRegimenStart) > 0)
  expect_true(nrow(pack$regimenTransition) > 0)
  expect_true(nrow(pack$regimenDuration) > 0)
  expect_true(nrow(pack$drugSchedule) > 0)
})


test_that("checkColumns catches missing columns", {
  dt <- data.table::data.table(a = 1, b = 2)
  expect_error(
    OncoRegimenTester:::checkColumns(dt, c("a", "b", "c"), "test.csv"),
    "missing required columns"
  )
})


test_that("checkProbabilityGroups detects bad sums", {
  dt <- data.table::data.table(
    group = c("A", "A", "B", "B"),
    prob = c(0.3, 0.3, 0.5, 0.5)
  )
  warnings <- OncoRegimenTester:::checkProbabilityGroups(
    dt, "group", "prob", "test.csv"
  )
  # Group A sums to 0.6, should warn
  expect_true(length(warnings) > 0)
  expect_true(any(grepl("group.*A", warnings, ignore.case = TRUE)))
})


test_that("parseBool handles various inputs", {
  expect_true(OncoRegimenTester:::parseBool("TRUE"))
  expect_true(OncoRegimenTester:::parseBool("1"))
  expect_true(OncoRegimenTester:::parseBool("yes"))
  expect_false(OncoRegimenTester:::parseBool("FALSE"))
  expect_false(OncoRegimenTester:::parseBool("0"))
  expect_false(OncoRegimenTester:::parseBool("no"))
})


test_that("validateParameterPack returns summary", {
  paramDir <- system.file("example_params", "mCRC", package = "OncoRegimenTester")
  if (paramDir == "") {
    paramDir <- file.path("inst", "example_params", "mCRC")
  }
  if (!dir.exists(paramDir)) skip("Example params not found")

  pack <- loadParameterPack(paramDir)
  result <- validateParameterPack(pack)

  expect_true(is.list(result))
  expect_true("valid" %in% names(result))
  expect_true("summary" %in% names(result))
  expect_true(length(result$summary) > 0)
})
