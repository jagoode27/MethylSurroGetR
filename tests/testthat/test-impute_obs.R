# Test file for impute_obs.R

test_that("impute_obs() function correctly handles 'mean' method", {
  # load sample data
  data(methyl_surro_miss)

  # generate expected values
  expected_vals <- methyl_surro_miss
  expected_vals$methyl[1, c(2, 3, 5)] <- mean(expected_vals$methyl[1, ], na.rm = TRUE)
  expected_vals$methyl[2, 3] <- mean(expected_vals$methyl[2, ], na.rm = TRUE)

  # run function
  function_vals <- impute_obs(methyl_surro_miss, method = "mean", min_nonmiss_prop = 0)

  # compare result
  expect_equal(function_vals, expected_vals)
})

test_that("impute_obs() function correctly handles 'mean' method with min_nonmiss_prop", {
  # load sample data
  data(methyl_surro_miss)

  # generate expected values
  expected_vals <- methyl_surro_miss
  expected_vals$methyl[2, 3] <- mean(expected_vals$methyl[2, ], na.rm = TRUE)

  # run function
  function_vals <- impute_obs(methyl_surro_miss, method = "mean", min_nonmiss_prop = 0.5)

  # compare result
  expect_equal(function_vals, expected_vals)
})

test_that("impute_obs() function correctly handles 'median' method", {
  # load sample data
  data(methyl_surro_miss)

  # generate expected values
  expected_vals <- methyl_surro_miss
  expected_vals$methyl[1, c(2, 3, 5)] <- median(expected_vals$methyl[1, ], na.rm = TRUE)
  expected_vals$methyl[2, 3] <- median(expected_vals$methyl[2, ], na.rm = TRUE)

  # run function
  function_vals <- impute_obs(methyl_surro_miss, method = "median", min_nonmiss_prop = 0)

  # compare result
  expect_equal(function_vals, expected_vals)
})

test_that("impute_obs() function correctly handles 'median' method with min_nonmiss_prop", {
  # load sample data
  data(methyl_surro_miss)

  # generate expected values
  expected_vals <- methyl_surro_miss
  expected_vals$methyl[2, 3] <- median(expected_vals$methyl[2, ], na.rm = TRUE)

  # run function
  function_vals <- impute_obs(methyl_surro_miss, method = "median", min_nonmiss_prop = 0.5)

  # compare result
  expect_equal(function_vals, expected_vals)
})

# NEW TESTS FOR COMPLETE COVERAGE

test_that("impute_obs() input validation works", {
  data(methyl_surro_miss)

  # Test non-methyl_surro object
  expect_error(impute_obs(list(methyl = matrix(1:4, 2, 2))),
               "Input must be an object of class 'methyl_surro'")

  # Test missing methyl component
  bad_surro <- methyl_surro_miss
  bad_surro$methyl <- NULL
  expect_error(impute_obs(bad_surro), "must contain a 'methyl' matrix component")

  # Test non-matrix methyl component
  bad_surro2 <- methyl_surro_miss
  bad_surro2$methyl <- c(1, 2, 3, 4)
  expect_error(impute_obs(bad_surro2), "must contain a 'methyl' matrix component")

  # Test invalid min_nonmiss_prop
  expect_error(impute_obs(methyl_surro_miss, min_nonmiss_prop = -0.1),
               "must be a numeric value between 0 and 1")
  expect_error(impute_obs(methyl_surro_miss, min_nonmiss_prop = 1.1),
               "must be a numeric value between 0 and 1")
  expect_error(impute_obs(methyl_surro_miss, min_nonmiss_prop = "0.5"),
               "must be a numeric value between 0 and 1")

  # Test invalid method
  expect_error(impute_obs(methyl_surro_miss, method = "invalid_method"),
               "'arg' should be one of")
})

test_that("impute_obs() handles edge cases", {
  # Test empty matrix
  empty_surro <- list(
    methyl = matrix(numeric(0), nrow = 0, ncol = 0),
    weights = c(cg1 = 1.0),
    intercept = 0.5
  )
  class(empty_surro) <- "methyl_surro"

  expect_error(impute_obs(empty_surro), "cannot be empty")

  # Test matrix with no missing values
  data(beta_matrix_comp)
  complete_surro <- list(
    methyl = beta_matrix_comp[1:5, ],
    weights = setNames(rep(1, 5), rownames(beta_matrix_comp)[1:5]),
    intercept = 0
  )
  class(complete_surro) <- "methyl_surro"

  expect_message(result <- impute_obs(complete_surro, verbose = TRUE),
                 "No missing values found")
  expect_identical(result, complete_surro)

  # Test completely missing probes only
  all_miss_matrix <- matrix(NA, nrow = 3, ncol = 4,
                            dimnames = list(c("cg1", "cg2", "cg3"),
                                            c("s1", "s2", "s3", "s4")))
  all_miss_surro <- list(
    methyl = all_miss_matrix,
    weights = setNames(c(1, 2, 3), c("cg1", "cg2", "cg3")),
    intercept = 0
  )
  class(all_miss_surro) <- "methyl_surro"

  expect_message(result <- impute_obs(all_miss_surro, verbose = TRUE),
                 "No probes meet the imputation threshold")
  expect_identical(result$methyl, all_miss_matrix)
})

test_that("impute_obs() verbose messaging works", {
  data(methyl_surro_miss)

  expect_message(
    impute_obs(methyl_surro_miss, method = "mean", verbose = TRUE),
    "Starting imputation analysis"
  )

  expect_message(
    impute_obs(methyl_surro_miss, method = "mean", verbose = TRUE),
    "Probe analysis:"
  )

  # Test non-verbose mode
  expect_message(
    result <- impute_obs(methyl_surro_miss, method = "mean", verbose = FALSE),
    "Imputed.*values.*probes"
  )
})

test_that("impute_obs() return_stats functionality works", {
  data(methyl_surro_miss)

  result <- impute_obs(methyl_surro_miss, method = "mean",
                       min_nonmiss_prop = 0, return_stats = TRUE)

  expect_true("imputation_stats" %in% names(result))
  expect_s3_class(result$imputation_stats, "imputation_stats")

  stats <- result$imputation_stats
  expect_true(all(c("method", "min_nonmiss_prop", "timestamp", "n_total_probes",
                    "n_complete_probes", "n_completely_missing_probes",
                    "n_partially_missing_probes", "n_probes_imputed",
                    "n_probes_skipped", "n_values_imputed",
                    "n_missing_before_imputation", "n_missing_after_imputation",
                    "imputation_rate", "probes_imputed", "probes_skipped") %in% names(stats)))

  expect_equal(stats$method, "mean")
  expect_equal(stats$min_nonmiss_prop, 0)
  expect_true(is.numeric(stats$n_total_probes))
  expect_true(stats$n_values_imputed > 0)
})

test_that("imputation_stats print method works", {
  data(methyl_surro_miss)

  result <- impute_obs(methyl_surro_miss, method = "median",
                       min_nonmiss_prop = 0.3, return_stats = TRUE)

  expect_output(print(result$imputation_stats), "Imputation Statistics")
  expect_output(print(result$imputation_stats), "Method: median")
  expect_output(print(result$imputation_stats), "Threshold:")
  expect_output(print(result$imputation_stats), "Probe Summary:")
  expect_output(print(result$imputation_stats), "Imputation Results:")
})

test_that("impute_obs() handles different min_nonmiss_prop thresholds", {
  data(methyl_surro_miss)

  # Test with high threshold (should impute fewer probes)
  result_high <- impute_obs(methyl_surro_miss, method = "mean",
                            min_nonmiss_prop = 0.8, return_stats = TRUE)

  # Test with low threshold (should impute more probes)
  result_low <- impute_obs(methyl_surro_miss, method = "mean",
                           min_nonmiss_prop = 0.2, return_stats = TRUE)

  expect_true(result_low$imputation_stats$n_probes_imputed >=
                result_high$imputation_stats$n_probes_imputed)
})

test_that("impute_obs() preserves object structure", {
  data(methyl_surro_miss)

  result <- impute_obs(methyl_surro_miss, method = "mean")

  expect_s3_class(result, "methyl_surro")
  expect_true(all(c("methyl", "weights", "intercept") %in% names(result)))
  expect_identical(result$weights, methyl_surro_miss$weights)
  expect_identical(result$intercept, methyl_surro_miss$intercept)
})

test_that("impute_obs() creates proper copy of input object", {
  data(methyl_surro_miss)
  original <- methyl_surro_miss

  result <- impute_obs(methyl_surro_miss, method = "mean")

  # Original object should be unchanged
  expect_identical(methyl_surro_miss, original)

  # Result should be different (imputed values)
  expect_false(identical(result$methyl, methyl_surro_miss$methyl))
})

test_that("impute_obs() handles single row/column edge cases", {
  # Single row with missing values
  single_row_matrix <- matrix(c(0.3, NA, 0.5), nrow = 1,
                              dimnames = list("cg1", c("s1", "s2", "s3")))
  single_row_surro <- list(
    methyl = single_row_matrix,
    weights = c(cg1 = 1.0),
    intercept = 0
  )
  class(single_row_surro) <- "methyl_surro"

  result <- impute_obs(single_row_surro, method = "mean")
  expected_value <- mean(c(0.3, 0.5))  # mean of non-missing values
  expect_equal(result$methyl[1, 2], expected_value)

  # Single column with missing values
  single_col_matrix <- matrix(c(0.2, NA, 0.8), ncol = 1,
                              dimnames = list(c("cg1", "cg2", "cg3"), "s1"))
  single_col_surro <- list(
    methyl = single_col_matrix,
    weights = setNames(c(1, 2, 3), c("cg1", "cg2", "cg3")),
    intercept = 0
  )
  class(single_col_surro) <- "methyl_surro"

  result2 <- impute_obs(single_col_surro, method = "median")
  # Row 2 has only 1 sample with NA, so can't be imputed with default threshold
  expect_true(is.na(result2$methyl[2, 1]))
})

test_that("impute_obs() skipped probes recommendations work", {
  data(methyl_surro_miss)

  # Test with very high threshold to create skipped probes
  expect_message(
    impute_obs(methyl_surro_miss, method = "mean",
               min_nonmiss_prop = 0.9, verbose = TRUE),
    "Recommendation:"
  )

  expect_message(
    impute_obs(methyl_surro_miss, method = "mean",
               min_nonmiss_prop = 0.9, verbose = TRUE),
    "reference_fill"
  )
})

test_that("create_imputation_stats helper function works correctly", {
  # This tests the internal helper function indirectly
  data(methyl_surro_miss)

  result <- impute_obs(methyl_surro_miss, method = "mean",
                       return_stats = TRUE, verbose = TRUE)

  stats <- result$imputation_stats

  # Test that all required fields are present
  required_fields <- c("method", "min_nonmiss_prop", "timestamp", "n_total_probes",
                       "imputation_rate", "probes_imputed", "probes_skipped")
  expect_true(all(required_fields %in% names(stats)))

  # Test that values make sense
  expect_true(stats$imputation_rate >= 0 && stats$imputation_rate <= 1)
  expect_true(stats$n_total_probes > 0)
  expect_true(length(stats$probes_imputed) == stats$n_probes_imputed)
  expect_true(length(stats$probes_skipped) == stats$n_probes_skipped)
})

test_that("impute_obs() handles extreme missing patterns", {
  # Create matrix with mixed patterns
  extreme_matrix <- matrix(c(
    0.1, NA, 0.3, 0.4, 0.5,  # 20% missing
    NA, NA, NA, NA, NA,      # 100% missing
    0.2, 0.3, NA, NA, NA,    # 60% missing
    0.4, 0.5, 0.6, 0.7, 0.8  # 0% missing
  ), nrow = 4, byrow = TRUE,
  dimnames = list(c("cg1", "cg2", "cg3", "cg4"), paste0("s", 1:5)))

  extreme_surro <- list(
    methyl = extreme_matrix,
    weights = setNames(c(1, 2, 3, 4), c("cg1", "cg2", "cg3", "cg4")),
    intercept = 0
  )
  class(extreme_surro) <- "methyl_surro"

  # Test with 50% threshold
  result <- impute_obs(extreme_surro, method = "mean",
                       min_nonmiss_prop = 0.5, return_stats = TRUE)

  # Only cg1 should be imputed (80% complete > 50% threshold)
  # cg2 is completely missing, cg3 is 40% complete < 50% threshold
  # cg4 is complete
  expect_equal(result$imputation_stats$n_probes_imputed, 1)
  expect_equal(result$imputation_stats$probes_imputed, "cg1")
  expect_equal(result$imputation_stats$n_probes_skipped, 1)  # cg3
  expect_equal(result$imputation_stats$n_completely_missing_probes, 1)  # cg2
})

test_that("impute_obs() handles large matrices efficiently", {
  # Create larger matrix for performance testing
  large_matrix <- matrix(rnorm(5000), nrow = 100, ncol = 50,
                         dimnames = list(paste0("cg", 1:100), paste0("s", 1:50)))

  # Add some missing values
  missing_indices <- sample(5000, 500)  # 10% missing
  large_matrix[missing_indices] <- NA

  large_surro <- list(
    methyl = large_matrix,
    weights = setNames(rep(1, 100), paste0("cg", 1:100)),
    intercept = 0
  )
  class(large_surro) <- "methyl_surro"

  # Should complete reasonably quickly
  start_time <- Sys.time()
  result <- impute_obs(large_surro, method = "mean", min_nonmiss_prop = 0.8)
  end_time <- Sys.time()

  expect_s3_class(result, "methyl_surro")
  expect_true(as.numeric(end_time - start_time) < 5)  # Should complete in under 5 seconds
})

test_that("impute_obs() methods produce different results", {
  data(methyl_surro_miss)

  # Create data where mean and median differ significantly
  skewed_matrix <- matrix(c(
    0.1, 0.1, 0.9, 0.9, NA,  # mean = 0.5, median = 0.5
    0.1, 0.2, 0.8, NA, NA   # mean = 0.37, median = 0.15
  ), nrow = 2, byrow = TRUE,
  dimnames = list(c("cg1", "cg2"), paste0("s", 1:5)))

  skewed_surro <- list(
    methyl = skewed_matrix,
    weights = setNames(c(1, 2), c("cg1", "cg2")),
    intercept = 0
  )
  class(skewed_surro) <- "methyl_surro"

  result_mean <- impute_obs(skewed_surro, method = "mean", min_nonmiss_prop = 0)
  result_median <- impute_obs(skewed_surro, method = "median", min_nonmiss_prop = 0)

  # For the second row, mean and median should be different
  expect_false(result_mean$methyl["cg2", "s4"] == result_median$methyl["cg2", "s4"])
})

test_that("impute_obs() handles numeric precision correctly", {
  # Test with very small numbers
  precision_matrix <- matrix(c(
    1e-10, NA, 1e-9, 1e-8,
    1e-15, 1e-14, NA, 1e-12
  ), nrow = 2, byrow = TRUE,
  dimnames = list(c("cg1", "cg2"), paste0("s", 1:4)))

  precision_surro <- list(
    methyl = precision_matrix,
    weights = setNames(c(1, 2), c("cg1", "cg2")),
    intercept = 0
  )
  class(precision_surro) <- "methyl_surro"

  result <- impute_obs(precision_surro, method = "mean")

  # Should handle small numbers without underflow
  expect_true(all(is.finite(result$methyl[!is.na(result$methyl)])))
  expect_false(any(result$methyl[!is.na(result$methyl)] == 0))
})
