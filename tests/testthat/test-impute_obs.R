# Test file for impute_obs.R

test_that("impute_obs() function correctly handles 'mean' method", {
  # load sample data
  data(beta_matrix_miss)
  data(wts_df)

  # create vectors from data frames
  wts_vec_lin <- setNames(wts_df$wt_lin, rownames(wts_df))

  # create methyl_surro object
  methyl_surro_miss <- surro_set(
    methyl = beta_matrix_miss,
    weights = wts_vec_lin,
    intercept = "Intercept"
  )

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
  data(beta_matrix_miss)
  data(wts_df)

  # create vectors from data frames
  wts_vec_lin <- setNames(wts_df$wt_lin, rownames(wts_df))

  # create methyl_surro object
  methyl_surro_miss <- surro_set(
    methyl = beta_matrix_miss,
    weights = wts_vec_lin,
    intercept = "Intercept"
  )

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
  data(beta_matrix_miss)
  data(wts_df)

  # create vectors from data frames
  wts_vec_lin <- setNames(wts_df$wt_lin, rownames(wts_df))

  # create methyl_surro object
  methyl_surro_miss <- surro_set(
    methyl = beta_matrix_miss,
    weights = wts_vec_lin,
    intercept = "Intercept"
  )

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
  data(beta_matrix_miss)
  data(wts_df)

  # create vectors from data frames
  wts_vec_lin <- setNames(wts_df$wt_lin, rownames(wts_df))

  # create methyl_surro object
  methyl_surro_miss <- surro_set(
    methyl = beta_matrix_miss,
    weights = wts_vec_lin,
    intercept = "Intercept"
  )

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
  data(beta_matrix_miss)
  data(wts_df)

  # create vectors from data frames
  wts_vec_lin <- setNames(wts_df$wt_lin, rownames(wts_df))

  # create methyl_surro object
  methyl_surro_miss <- surro_set(
    methyl = beta_matrix_miss,
    weights = wts_vec_lin,
    intercept = "Intercept"
  )

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
  expect_identical(result$methyl, complete_surro$methyl)
})

test_that("impute_obs() verbose messaging works", {
  data(beta_matrix_miss)
  data(wts_df)

  # create vectors from data frames
  wts_vec_lin <- setNames(wts_df$wt_lin, rownames(wts_df))

  # create methyl_surro object
  methyl_surro_miss <- surro_set(
    methyl = beta_matrix_miss,
    weights = wts_vec_lin,
    intercept = "Intercept"
  )

  expect_message(impute_obs(methyl_surro_miss, verbose = TRUE),
                 "Starting imputation analysis")

  expect_message(impute_obs(methyl_surro_miss, verbose = TRUE),
                 "Found.*missing values")

  expect_message(impute_obs(methyl_surro_miss, verbose = TRUE),
                 "Probe analysis:")

  # Test non-verbose mode
  expect_message(impute_obs(methyl_surro_miss, verbose = FALSE),
                 "Imputed.*values.*probes using.*method")
})

test_that("impute_obs() return_stats functionality works", {
  data(beta_matrix_miss)
  data(wts_df)

  # create vectors from data frames
  wts_vec_lin <- setNames(wts_df$wt_lin, rownames(wts_df))

  # create methyl_surro object
  methyl_surro_miss <- surro_set(
    methyl = beta_matrix_miss,
    weights = wts_vec_lin,
    intercept = "Intercept"
  )

  result <- impute_obs(methyl_surro_miss, return_stats = TRUE)

  expect_true("imputation_stats" %in% names(result))
  expect_s3_class(result$imputation_stats, "imputation_stats")

  stats <- result$imputation_stats
  required_fields <- c("method", "min_nonmiss_prop", "timestamp", "n_total_probes",
                       "n_complete_probes", "n_completely_missing_probes",
                       "n_partially_missing_probes", "n_probes_imputed",
                       "n_probes_skipped", "n_values_imputed",
                       "n_missing_before_imputation", "n_missing_after_imputation",
                       "imputation_rate", "probes_imputed", "probes_skipped")

  expect_true(all(required_fields %in% names(stats)))
  expect_equal(stats$method, "mean")
  expect_true(stats$imputation_rate >= 0 && stats$imputation_rate <= 1)
})

test_that("impute_obs() print method for imputation_stats works", {
  data(beta_matrix_miss)
  data(wts_df)

  # create vectors from data frames
  wts_vec_lin <- setNames(wts_df$wt_lin, rownames(wts_df))

  # create methyl_surro object
  methyl_surro_miss <- surro_set(
    methyl = beta_matrix_miss,
    weights = wts_vec_lin,
    intercept = "Intercept"
  )

  result <- impute_obs(methyl_surro_miss, return_stats = TRUE)

  expect_output(print(result$imputation_stats), "Imputation Statistics")
  expect_output(print(result$imputation_stats), "Method:")
  expect_output(print(result$imputation_stats), "Probe Summary:")
  expect_output(print(result$imputation_stats), "Imputation Results:")
})

test_that("impute_obs() preserves object structure", {
  data(beta_matrix_miss)
  data(wts_df)

  # create vectors from data frames
  wts_vec_lin <- setNames(wts_df$wt_lin, rownames(wts_df))

  # create methyl_surro object
  methyl_surro_miss <- surro_set(
    methyl = beta_matrix_miss,
    weights = wts_vec_lin,
    intercept = "Intercept"
  )

  result <- impute_obs(methyl_surro_miss)

  expect_s3_class(result, "methyl_surro")
  expect_true(all(c("methyl", "weights", "intercept") %in% names(result)))
  expect_identical(result$weights, methyl_surro_miss$weights)
  expect_identical(result$intercept, methyl_surro_miss$intercept)
})

test_that("impute_obs() creates proper copy of input object", {
  data(beta_matrix_miss)
  data(wts_df)

  # create vectors from data frames
  wts_vec_lin <- setNames(wts_df$wt_lin, rownames(wts_df))

  # create methyl_surro object
  methyl_surro_miss <- surro_set(
    methyl = beta_matrix_miss,
    weights = wts_vec_lin,
    intercept = "Intercept"
  )

  original <- methyl_surro_miss

  result <- impute_obs(methyl_surro_miss)

  # Original object should be unchanged
  expect_identical(methyl_surro_miss, original)

  # Result should be different (imputed values)
  expect_false(identical(result$methyl, methyl_surro_miss$methyl))
})

test_that("impute_obs() handles different missing patterns", {
  data(wts_df)

  # create vectors from data frames
  wts_vec_lin <- setNames(wts_df$wt_lin, rownames(wts_df))

  # Create matrix with different missing patterns
  mixed_missing_matrix <- matrix(c(
    0.1, 0.2, 0.3, 0.4, 0.5,  # complete
    NA, NA, NA, NA, NA,       # completely missing
    0.2, NA, 0.4, 0.5, 0.6,   # partially missing (20% missing)
    0.3, NA, NA, NA, 0.7      # partially missing (60% missing)
  ), nrow = 4, byrow = TRUE,
  dimnames = list(c("cg02", "cg03", "cg07", "cg08"), paste0("samp", 1:5)))

  mixed_surro <- list(
    methyl = mixed_missing_matrix,
    weights = wts_vec_lin[c("cg02", "cg03", "cg07", "cg08")],
    intercept = wts_vec_lin["Intercept"]
  )
  class(mixed_surro) <- "methyl_surro"

  # Test with threshold that should exclude the 60% missing probe
  result <- impute_obs(mixed_surro, min_nonmiss_prop = 0.5)

  # Should impute the 20% missing probe but not the 60% missing one
  expect_false(is.na(result$methyl["cg07", "samp2"]))  # Should be imputed
  expect_true(is.na(result$methyl["cg08", "samp2"]))   # Should remain NA
})

test_that("impute_obs() handles edge cases with thresholds", {
  data(wts_df)

  # create vectors from data frames
  wts_vec_lin <- setNames(wts_df$wt_lin, rownames(wts_df))

  # Create matrix where no probes meet the threshold
  all_high_missing <- matrix(c(
    NA, 0.2, NA, NA, NA,  # 80% missing
    NA, NA, 0.4, NA, NA,  # 80% missing
    NA, NA, NA, 0.5, NA   # 80% missing
  ), nrow = 3, byrow = TRUE,
  dimnames = list(c("cg02", "cg07", "cg08"), paste0("samp", 1:5)))

  high_missing_surro <- list(
    methyl = all_high_missing,
    weights = wts_vec_lin[c("cg02", "cg07", "cg08")],
    intercept = wts_vec_lin["Intercept"]
  )
  class(high_missing_surro) <- "methyl_surro"

  # With high threshold, no probes should be imputed
  result <- impute_obs(high_missing_surro, min_nonmiss_prop = 0.9, verbose = TRUE)

  expect_message(impute_obs(high_missing_surro, min_nonmiss_prop = 0.9, verbose = TRUE),
                 "No probes meet the imputation threshold")

  # Matrix should be unchanged
  expect_identical(result$methyl, high_missing_surro$methyl)
})

test_that("impute_obs() numerical precision is maintained", {
  data(wts_df)

  # create vectors from data frames
  wts_vec_lin <- setNames(wts_df$wt_lin, rownames(wts_df))

  # Create matrix with high precision values
  precision_matrix <- matrix(c(
    0.123456789, NA, 0.987654321,
    0.111111111, 0.222222222, 0.333333333
  ), nrow = 2, byrow = TRUE,
  dimnames = list(c("cg02", "cg07"), c("s1", "s2", "s3")))

  precision_surro <- list(
    methyl = precision_matrix,
    weights = wts_vec_lin[c("cg02", "cg07")],
    intercept = wts_vec_lin["Intercept"]
  )
  class(precision_surro) <- "methyl_surro"

  result <- impute_obs(precision_surro, method = "mean")

  # Should maintain high precision for mean calculation
  expected_mean <- mean(c(0.123456789, 0.987654321))
  expect_equal(result$methyl["cg02", "s2"], expected_mean, tolerance = 1e-10)
})

test_that("impute_obs() handles single row/column matrices", {
  data(wts_df)

  # create vectors from data frames
  wts_vec_lin <- setNames(wts_df$wt_lin, rownames(wts_df))

  # Single row matrix
  single_row <- matrix(c(0.3, NA, 0.5), nrow = 1,
                       dimnames = list("cg02", c("s1", "s2", "s3")))
  single_row_surro <- list(
    methyl = single_row,
    weights = wts_vec_lin["cg02"],
    intercept = wts_vec_lin["Intercept"]
  )
  class(single_row_surro) <- "methyl_surro"

  result1 <- impute_obs(single_row_surro)
  expect_equal(result1$methyl[1, 2], mean(c(0.3, 0.5)))

  # Single column matrix
  single_col <- matrix(c(0.2, NA, 0.8), ncol = 1,
                       dimnames = list(c("cg02", "cg07", "cg08"), "s1"))
  single_col_surro <- list(
    methyl = single_col,
    weights = wts_vec_lin[c("cg02", "cg07", "cg08")],
    intercept = wts_vec_lin["Intercept"]
  )
  class(single_col_surro) <- "methyl_surro"

  # Single column means no imputation possible (can't calculate row means with one value)
  expect_message(result2 <- impute_obs(single_col_surro, verbose = TRUE),
                 "No probes meet the imputation threshold")
})

test_that("impute_obs() handles large matrices efficiently", {
  # Create larger test dataset
  large_matrix <- matrix(rnorm(10000), nrow = 100, ncol = 100,
                         dimnames = list(paste0("cg", 1:100), paste0("s", 1:100)))

  # Add random missing patterns
  missing_indices <- sample(10000, 1000)  # 10% missing
  large_matrix[missing_indices] <- NA

  large_weights <- setNames(rep(1, 100), paste0("cg", 1:100))
  large_surro <- list(
    methyl = large_matrix,
    weights = large_weights,
    intercept = 0
  )
  class(large_surro) <- "methyl_surro"

  # Should complete reasonably quickly
  start_time <- Sys.time()
  result <- impute_obs(large_surro, min_nonmiss_prop = 0.5)
  end_time <- Sys.time()

  expect_s3_class(result, "methyl_surro")
  expect_true(as.numeric(end_time - start_time) < 5)  # Should complete in under 5 seconds
})

test_that("impute_obs() recommendations are helpful", {
  data(beta_matrix_miss)
  data(wts_df)

  # create vectors from data frames
  wts_vec_lin <- setNames(wts_df$wt_lin, rownames(wts_df))

  # create methyl_surro object
  methyl_surro_miss <- surro_set(
    methyl = beta_matrix_miss,
    weights = wts_vec_lin,
    intercept = "Intercept"
  )

  # Test that recommendations are provided for skipped probes
  expect_message(
    impute_obs(methyl_surro_miss, min_nonmiss_prop = 0.5, verbose = TRUE),
    "Recommendation:.*probes were skipped"
  )

  # Test recommendations for completely missing probes
  expect_message(
    impute_obs(methyl_surro_miss, min_nonmiss_prop = 0, verbose = TRUE),
    "completely missing probes detected.*reference_fill"
  )
})
