ex <- system.file("extdata", "example_markers.csv", package = "omg")

test_that("bundled reference data is present", {
  expect_true(file.exists(omg_default_markers()))
  expect_true(file.exists(omg_default_orthogroups()))
})

test_that("omg() runs end-to-end and returns expected structure", {
  res <- omg(ex, write_files = FALSE, verbose = FALSE)
  expect_type(res, "list")
  expect_named(res, c("predictions", "pairwise", "combined", "extract", "out_dir"))
  expect_s3_class(res$predictions, "data.frame")
  expect_true(all(c("query_cluster", "consolidated_cell_type",
                    "prediction_confidence") %in% names(res$predictions)))
})

test_that("known root clusters get the expected cell types", {
  res <- omg(ex, write_files = FALSE, verbose = FALSE)
  pred <- setNames(res$predictions$consolidated_cell_type,
                   res$predictions$query_cluster)
  # truth: 0 = Lateral root cap, 1 = Root cortex, 3 = Xylem, 4 = Root endodermis
  expect_equal(pred[["1"]], "Cortex")
  expect_equal(pred[["3"]], "Xylem")
  expect_equal(pred[["4"]], "Endodermis")
})

test_that("omg_check_reference reports a high mapping rate for the example", {
  chk <- omg_check_reference(ex)
  expect_gt(chk$rate, 0.8)
})

test_that("omg_list_species and helpers work", {
  sp <- omg_list_species()
  expect_true("arabidopsis_thaliana" %in% sp$marker_species)
  expect_s3_class(omg_paper_15species(), "data.frame")
  expect_type(omg_cell_type_groups(), "list")
})
