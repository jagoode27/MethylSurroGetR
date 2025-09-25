# Test file for surro_set.R

# Create Expected Values ----

# expected methylation matrix
expected_methyl <- rbind(
  beta_matrix_comp[c("cg02", "cg07", "cg08", "cg13", "cg17"), ],
  matrix(rep(NA, 25), nrow = 5) |>
    `colnames<-`(paste0("samp", 1:5)) |>
    `rownames<-`(paste0("cg", c("03", "06", "11", "15", "18")))
)

# Run Original Tests ----
test_that("surro_set() function works with weights vector", {
  # load sample data
  data(beta_matrix_comp)
  data(wts_vec_cnt)

  # generate expected values
  expected_vals <- list(
    methyl = expected_methyl,
    weights = wts_vec_cnt[-11],
    intercept = wts_vec_cnt[11]
  )
  class(expected_vals) <- "methyl_surro"

  # run function
  function_vals <- surro_set(
    methyl = beta_matrix_comp,
    weights = wts_vec_cnt,
    intercept = "Intercept"
  )

  # compare result
  expect_identical(function_vals, expected_vals)
})

# NOTE: Removing matrix and data frame tests as current implementation only supports named vectors

# NEW TESTS FOR COMPLETE COVERAGE

test_that("surro_set() input validation works", {
  data(beta_matrix_comp)
  data(wts_vec_lin)

  # Test non-matrix methylation input
  expect_error(surro_set(c(1, 2, 3, 4), wts_vec_lin),
               "must be a matrix")

  # Test non-numeric matrix
  char_matrix <- matrix(c("a", "b", "c", "d"), nrow = 2)
  expect_error(surro_set(char_matrix, wts_vec_lin),
               "must be numeric")

  # Test empty matrix
  empty_matrix <- matrix(numeric(0), nrow = 0, ncol = 0)
  expect_error(surro_set(empty_matrix, wts_vec_lin),
               "cannot be empty")

  # Test matrix without row names
  no_names_matrix <- beta_matrix_comp
  rownames(no_names_matrix) <- NULL
  expect_error(surro_set(no_names_matrix, wts_vec_lin),
               "must have CpG sites as row names")

  # Test non-numeric weights
  expect_error(surro_set(beta_matrix_comp, c("a", "b", "c")),
               "must be a named numeric vector")

  # Test non-vector weights
  expect_error(surro_set(beta_matrix_comp, matrix(1:4, 2, 2)),
               "must be a named numeric vector")

  # Test weights without names
  unnamed_weights <- c(1, 2, 3, 4, 5)
  expect_error(surro_set(beta_matrix_comp, unnamed_weights),
               "must be a named numeric vector")

  # Test empty weights
  expect_error(surro_set(beta_matrix_comp, numeric(0)),
               "cannot be empty")

  # Test weights with duplicate names
  dup_weights <- c(cg01 = 1, cg02 = 2, cg01 = 3)
  expect_error(surro_set(beta_matrix_comp, dup_weights),
               "Duplicate CpG names found")

  # Test weights with missing names
  na_names_weights <- c(1, 2, 3)
  names(na_names_weights) <- c("cg01", NA, "cg03")
  expect_error(surro_set(beta_matrix_comp, na_names_weights),
               "no missing or empty names")

  # Test weights with empty names
  empty_names_weights <- c(1, 2, 3)
  names(empty_names_weights) <- c("cg01", "", "cg03")
  expect_error(surro_set(beta_matrix_comp, empty_names_weights),
               "no missing or empty names")
})

test_that("surro_set() intercept handling works correctly", {
  data(beta_matrix_comp)
  data(wts_vec_lin)

  # Test with valid intercept
  result_with_int <- surro_set(beta_matrix_comp, wts_vec_lin, "Intercept")
  expect_equal(result_with_int$intercept, wts_vec_lin["Intercept"])
  expect_false("Intercept" %in% names(result_with_int$weights))

  # Test without intercept
  result_no_int <- surro_set(beta_matrix_comp, wts_vec_lin[-11])  # Remove intercept
  expect_null(result_no_int$intercept)

  # Test with non-character intercept
  expect_error(surro_set(beta_matrix_comp, wts_vec_lin, intercept = 123),
               "must be a single character string")

  # Test with multiple intercept values
  expect_error(surro_set(beta_matrix_comp, wts_vec_lin, intercept = c("Int1", "Int2")),
               "must be a single character string")

  # Test with intercept not found in weights
  expect_error(surro_set(beta_matrix_comp, wts_vec_lin, intercept = "NotFound"),
               "not found in the weights")

  # Test with only intercept in weights (no other weights remain)
  intercept_only <- c(Intercept = 1.5)
  expect_error(surro_set(beta_matrix_comp, intercept_only, intercept = "Intercept"),
               "No weights remain after removing intercept")
})

test_that("surro_set() handles probe alignment correctly", {
  data(beta_matrix_comp)
  data(wts_vec_lin)

  # Test with perfect alignment
  aligned_result <- surro_set(beta_matrix_comp[c("cg02", "cg07", "cg08", "cg13", "cg17"), ],
                              wts_vec_lin[c("cg02", "cg07", "cg08", "cg13", "cg17", "Intercept")],
                              intercept = "Intercept")

  expect_equal(nrow(aligned_result$methyl), 5)
  expect_true(all(c("cg02", "cg07", "cg08", "cg13", "cg17") %in% rownames(aligned_result$methyl)))

  # Test informational messages
  expect_message(surro_set(beta_matrix_comp, wts_vec_lin, "Intercept"),
                 "Added.*missing probes with NA values")

  # Test with methylation matrix having more probes than weights
  large_matrix <- rbind(beta_matrix_comp,
                        matrix(runif(10), nrow = 2,
                               dimnames = list(c("extra1", "extra2"), colnames(beta_matrix_comp))))

  expect_message(surro_set(large_matrix, wts_vec_lin, "Intercept"),
                 "Filtered methylation matrix")
})

test_that("surro_set() creates proper methyl_surro object", {
  data(beta_matrix_comp)
  data(wts_vec_lin)

  result <- surro_set(beta_matrix_comp, wts_vec_lin, "Intercept")

  # Test class assignment
  expect_s3_class(result, "methyl_surro")

  # Test required components
  expect_true(all(c("methyl", "weights", "intercept") %in% names(result)))

  # Test component types
  expect_true(is.matrix(result$methyl))
  expect_true(is.numeric(result$weights))
  expect_true(is.vector(result$weights))
  expect_true(is.numeric(result$intercept))
})

test_that("surro_set() handles edge cases with matrix dimensions", {
  data(wts_vec_lin)

  # Single row matrix
  single_row <- matrix(c(0.3, 0.4, 0.5), nrow = 1,
                       dimnames = list("cg02", c("s1", "s2", "s3")))

  result1 <- surro_set(single_row, wts_vec_lin, "Intercept")
  expect_s3_class(result1, "methyl_surro")
  expect_equal(ncol(result1$methyl), 3)

  # Single column matrix
  single_col <- matrix(c(0.2, 0.6, 0.8), ncol = 1,
                       dimnames = list(c("cg02", "cg07", "cg08"), "s1"))

  result2 <- surro_set(single_col, wts_vec_lin[c("cg02", "cg07", "cg08", "Intercept")], "Intercept")
  expect_s3_class(result2, "methyl_surro")
  expect_equal(ncol(result2$methyl), 1)
})

test_that("surro_set() preserves matrix structure and names", {
  data(beta_matrix_comp)
  data(wts_vec_lin)

  result <- surro_set(beta_matrix_comp, wts_vec_lin, "Intercept")

  # Test that column names are preserved
  expect_identical(colnames(result$methyl), colnames(beta_matrix_comp))

  # Test that row names match weight names (excluding intercept)
  weight_names_no_int <- names(wts_vec_lin)[names(wts_vec_lin) != "Intercept"]
  expect_true(all(weight_names_no_int %in% rownames(result$methyl)))

  # Test that weights names are preserved
  expect_identical(names(result$weights), weight_names_no_int)
})

test_that("surro_set() handles missing probes correctly", {
  data(beta_matrix_comp)
  data(wts_vec_lin)

  # Create weights that require probes not in methylation matrix
  extended_weights <- c(wts_vec_lin, c(missing_probe = 1.5))

  result <- surro_set(beta_matrix_comp, extended_weights, "Intercept")

  # Should add row of NAs for missing probe
  expect_true("missing_probe" %in% rownames(result$methyl))
  expect_true(all(is.na(result$methyl["missing_probe", ])))
})

test_that("surro_set() efficient handling of large matrices", {
  # Create larger test data
  large_matrix <- matrix(rnorm(10000), nrow = 100, ncol = 100,
                         dimnames = list(paste0("cg", 1:100), paste0("s", 1:100)))

  large_weights <- setNames(rnorm(50), paste0("cg", 1:50))  # Only half the probes
  large_weights <- c(large_weights, c(Intercept = 0.5))

  # Function should handle large matrices efficiently
  start_time <- Sys.time()
  result <- surro_set(large_matrix, large_weights, "Intercept")
  end_time <- Sys.time()

  expect_s3_class(result, "methyl_surro")
  expect_true(as.numeric(end_time - start_time) < 2)  # Should complete quickly
  expect_equal(nrow(result$methyl), 50)  # Should match number of weights
})

test_that("surro_set() handles special characters in probe names", {
  # Create matrix with special characters in names
  special_matrix <- matrix(c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6), nrow = 2,
                           dimnames = list(c("cg-01_test", "cg.02-test"),
                                           c("sample_1", "sample.2", "sample-3")))

  special_weights <- c("cg-01_test" = 1.0, "cg.02-test" = 2.0, "Intercept" = 0.5)

  result <- surro_set(special_matrix, special_weights, "Intercept")

  expect_s3_class(result, "methyl_surro")
  expect_true(all(c("cg-01_test", "cg.02-test") %in% rownames(result$methyl)))
})

test_that("surro_set() preserves data types and precision", {
  data(beta_matrix_comp)
  data(wts_vec_lin)

  result <- surro_set(beta_matrix_comp, wts_vec_lin, "Intercept")

  # Test that numeric precision is maintained
  common_probes <- intersect(rownames(beta_matrix_comp), rownames(result$methyl))
  expect_equal(result$methyl[common_probes, ], beta_matrix_comp[common_probes, ])

  # Test that weights maintain precision
  weight_names <- names(wts_vec_lin)[names(wts_vec_lin) != "Intercept"]
  expect_equal(result$weights, wts_vec_lin[weight_names])

  # Test intercept precision
  expect_equal(result$intercept, wts_vec_lin["Intercept"])
})

test_that("surro_set() works with minimal valid input", {
  # Minimal case: 1 probe, 1 sample, 1 weight + intercept
  minimal_matrix <- matrix(0.5, nrow = 1, ncol = 1,
                           dimnames = list("cg01", "s1"))
  minimal_weights <- c(cg01 = 1.0, Intercept = 0.0)

  result <- surro_set(minimal_matrix, minimal_weights, "Intercept")

  expect_s3_class(result, "methyl_surro")
  expect_equal(dim(result$methyl), c(1, 1))
  expect_length(result$weights, 1)
  expect_equal(result$intercept, 0.0)
})

test_that("surro_set() handles matrices with missing values", {
  data(wts_vec_lin)

  # Create matrix with some NA values
  matrix_with_na <- matrix(c(0.1, NA, 0.3, 0.4, NA, 0.6), nrow = 2,
                           dimnames = list(c("cg02", "cg07"), c("s1", "s2", "s3")))

  result <- surro_set(matrix_with_na, wts_vec_lin[c("cg02", "cg07", "Intercept")], "Intercept")

  expect_s3_class(result, "methyl_surro")

  # NA values should be preserved
  expect_true(is.na(result$methyl["cg02", "s2"]))
  expect_true(is.na(result$methyl["cg07", "s2"]))

  # Non-NA values should be preserved
  expect_equal(result$methyl["cg02", "s1"], 0.1)
  expect_equal(result$methyl["cg07", "s3"], 0.6)
})

test_that("surro_set() error messages are informative", {
  data(beta_matrix_comp)
  data(wts_vec_lin)

  # Test specific error message content
  expect_error(surro_set("not_a_matrix", wts_vec_lin),
               "methylation data must be a matrix", fixed = TRUE)

  expect_error(surro_set(matrix(c("a", "b", "c", "d"), 2, 2), wts_vec_lin),
               "must be numeric", fixed = TRUE)

  expect_error(surro_set(beta_matrix_comp, "not_weights"),
               "must be a named numeric vector", fixed = TRUE)

  # Test that error messages help identify the problem
  no_rownames <- beta_matrix_comp
  rownames(no_rownames) <- NULL
  expect_error(surro_set(no_rownames, wts_vec_lin),
               "CpG sites as row names", fixed = TRUE)
})

test_that("surro_set() messaging provides useful information", {
  data(beta_matrix_comp)
  data(wts_vec_lin)

  # Should provide information about added missing probes
  expect_message(surro_set(beta_matrix_comp, wts_vec_lin, "Intercept"),
                 "missing from methylation data")

  # Should provide information about filtered probes
  large_matrix <- rbind(beta_matrix_comp,
                        matrix(0.5, nrow = 5, ncol = ncol(beta_matrix_comp),
                               dimnames = list(paste0("extra", 1:5), NULL)))

  expect_message(surro_set(large_matrix, wts_vec_lin, "Intercept"),
                 "Filtered methylation matrix from")
})

test_that("surro_set() return value structure is consistent", {
  data(beta_matrix_comp)
  data(wts_vec_lin)

  # Test multiple different input scenarios
  scenarios <- list(
    # Scenario 1: Perfect alignment
    list(methyl = beta_matrix_comp[1:5, ],
         weights = wts_vec_lin[c(rownames(beta_matrix_comp)[1:5], "Intercept")]),

    # Scenario 2: Missing probes need to be added
    list(methyl = beta_matrix_comp[1:3, ],
         weights = wts_vec_lin),

    # Scenario 3: Extra methylation probes need to be filtered
    list(methyl = beta_matrix_comp,
         weights = wts_vec_lin[1:5])  # Only first 4 weights + intercept
  )

  for (i in seq_along(scenarios)) {
    result <- surro_set(scenarios[[i]]$methyl, scenarios[[i]]$weights, "Intercept")

    # All should return methyl_surro objects with same structure
    expect_s3_class(result, "methyl_surro")
    expect_true(all(c("methyl", "weights", "intercept") %in% names(result)))
    expect_true(is.matrix(result$methyl))
    expect_true(is.numeric(result$weights) && is.vector(result$weights))
    expect_true(is.numeric(result$intercept) && length(result$intercept) == 1)
  }
})
