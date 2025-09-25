# Test file for methyl_miss.R

test_that("methyl_miss() function works", {
  # load sample data
  data(beta_matrix_miss)
  data(wts_vec_lin)

  # generate expected values
  expected_vals <- list(
    missing_obs = c(0.6, 0.2) |>
      `names<-`(c("cg02", "cg07")),
    missing_probes = c("cg03", "cg06", "cg11", "cg15", "cg18")
  )

  # run function
  function_vals <- surro_set(
    methyl = beta_matrix_miss,
    weights = wts_vec_lin,  # Fixed: was wts_vec_cnt
    intercept = "Intercept"
  ) |>
    methyl_miss()

  # compare result - only check the main components, not summary
  expect_identical(function_vals$missing_obs, expected_vals$missing_obs)
  expect_identical(function_vals$missing_probes, expected_vals$missing_probes)
})

# NEW TESTS FOR COMPLETE COVERAGE

test_that("methyl_miss() input validation works", {
  # Test non-methyl_surro object
  expect_error(methyl_miss(list(methyl = matrix(1:4, 2, 2))),
               "Input must be an object of class 'methyl_surro'")

  # Test missing methyl component
  bad_surro <- list(weights = c(1, 2), intercept = 0)
  class(bad_surro) <- "methyl_surro"
  expect_error(methyl_miss(bad_surro), "must contain a 'methyl' matrix component")

  # Test non-matrix methyl component
  bad_surro2 <- list(methyl = c(1, 2, 3), weights = c(1, 2), intercept = 0)
  class(bad_surro2) <- "methyl_surro"
  expect_error(methyl_miss(bad_surro2), "must contain a 'methyl' matrix component")

  # Test empty matrix
  empty_surro <- list(
    methyl = matrix(numeric(0), nrow = 0, ncol = 0),
    weights = c(cg1 = 1.0),
    intercept = 0.5
  )
  class(empty_surro) <- "methyl_surro"
  expect_error(methyl_miss(empty_surro), "cannot be empty")
})

test_that("methyl_miss() handles complete data correctly", {
  data(beta_matrix_comp)

  complete_surro <- list(
    methyl = beta_matrix_comp,
    weights = setNames(rep(1, nrow(beta_matrix_comp)), rownames(beta_matrix_comp)),
    intercept = 0
  )
  class(complete_surro) <- "methyl_surro"

  result <- methyl_miss(complete_surro)

  expect_s3_class(result, "methyl_miss")
  expect_length(result$missing_obs, 0)
  expect_length(result$missing_probes, 0)
  expect_equal(result$summary$n_complete_probes, nrow(beta_matrix_comp))
  expect_equal(result$summary$overall_missing_rate, 0)
})

test_that("methyl_miss() handles all missing data correctly", {
  # Create matrix with all NA values
  all_na_matrix <- matrix(NA, nrow = 3, ncol = 4,
                          dimnames = list(c("cg1", "cg2", "cg3"),
                                          c("s1", "s2", "s3", "s4")))

  all_na_surro <- list(
    methyl = all_na_matrix,
    weights = setNames(c(1, 2, 3), c("cg1", "cg2", "cg3")),
    intercept = 0
  )
  class(all_na_surro) <- "methyl_surro"

  result <- methyl_miss(all_na_surro)

  expect_length(result$missing_obs, 0)  # No partial missing
  expect_length(result$missing_probes, 3)  # All probes completely missing
  expect_equal(result$summary$overall_missing_rate, 1.0)
  expect_equal(result$summary$n_completely_missing_probes, 3)
})

test_that("methyl_miss() summary statistics are correct", {
  data(methyl_surro_miss)

  result <- methyl_miss(methyl_surro_miss)

  # Test that summary contains all required fields
  required_fields <- c("total_probes", "total_samples", "n_complete_probes",
                       "n_missing_obs", "n_missing_probes", "overall_missing_rate",
                       "missing_obs_rate", "missing_probes_rate", "complete_probes_rate")

  expect_true(all(required_fields %in% names(result$summary)))

  # Test that rates sum correctly
  expect_equal(
    result$summary$complete_probes_rate +
      result$summary$missing_obs_rate +
      result$summary$missing_probes_rate,
    1.0,
    tolerance = 1e-10
  )

  # Test that counts sum correctly
  expect_equal(
    result$summary$n_complete_probes +
      result$summary$n_missing_obs +
      result$summary$n_missing_probes,
    result$summary$total_probes
  )

  # Test that overall missing rate makes sense
  expect_true(result$summary$overall_missing_rate >= 0)
  expect_true(result$summary$overall_missing_rate <= 1)
})

test_that("methyl_miss() print method works", {
  data(methyl_surro_miss)

  result <- methyl_miss(methyl_surro_miss)

  expect_output(print(result), "Missing Data Summary for methyl_surro Object")
  expect_output(print(result), "Total probes:")
  expect_output(print(result), "Total samples:")
  expect_output(print(result), "Complete probes:")
  expect_output(print(result), "Overall missing rate:")
})

test_that("methyl_miss() print method handles long lists correctly", {
  # Create a surro object with many missing probes to test truncation
  large_matrix <- matrix(rnorm(200), nrow = 20, ncol = 10,
                         dimnames = list(paste0("cg", 1:20), paste0("s", 1:10)))

  # Make many probes completely missing
  large_matrix[11:20, ] <- NA

  # Make some probes partially missing
  large_matrix[1:5, 1:3] <- NA

  large_surro <- list(
    methyl = large_matrix,
    weights = setNames(rep(1, 20), paste0("cg", 1:20)),
    intercept = 0
  )
  class(large_surro) <- "methyl_surro"

  result <- methyl_miss(large_surro)

  # Test that print method truncates long lists
  expect_output(print(result), "... and \\d+ more")
})

test_that("methyl_miss() handles edge cases with missing proportions", {
  # Create matrix with specific missing patterns
  edge_matrix <- matrix(c(
    0.1, NA, 0.3, 0.4,      # 25% missing (1/4)
    NA, NA, 0.5, 0.6,       # 50% missing (2/4)
    0.7, 0.8, 0.9, NA,      # 25% missing (1/4)
    NA, NA, NA, NA,         # 100% missing (4/4)
    0.2, 0.3, 0.4, 0.5      # 0% missing (0/4)
  ), nrow = 5, byrow = TRUE,
  dimnames = list(paste0("cg", 1:5), paste0("s", 1:4)))

  edge_surro <- list(
    methyl = edge_matrix,
    weights = setNames(c(1, 2, 3, 4, 5), paste0("cg", 1:5)),
    intercept = 0
  )
  class(edge_surro) <- "methyl_surro"

  result <- methyl_miss(edge_surro)

  # Check specific missing proportions
  expect_equal(result$missing_obs["cg1"], 0.25)
  expect_equal(result$missing_obs["cg2"], 0.50)
  expect_equal(result$missing_obs["cg3"], 0.25)
  expect_true("cg4" %in% result$missing_probes)
  expect_false("cg5" %in% names(result$missing_obs))  # complete probe
})

test_that("methyl_miss() works with single row/column matrices", {
  # Single row matrix
  single_row <- matrix(c(0.3, NA, 0.5), nrow = 1,
                       dimnames = list("cg1", c("s1", "s2", "s3")))
  single_row_surro <- list(
    methyl = single_row,
    weights = c(cg1 = 1.0),
    intercept = 0
  )
  class(single_row_surro) <- "methyl_surro"

  result1 <- methyl_miss(single_row_surro)
  expect_equal(result1$summary$total_probes, 1)
  expect_equal(result1$summary$total_samples, 3)
  expect_equal(result1$missing_obs["cg1"], 1/3)

  # Single column matrix
  single_col <- matrix(c(0.2, NA, 0.8), ncol = 1,
                       dimnames = list(c("cg1", "cg2", "cg3"), "s1"))
  single_col_surro <- list(
    methyl = single_col,
    weights = setNames(c(1, 2, 3), c("cg1", "cg2", "cg3")),
    intercept = 0
  )
  class(single_col_surro) <- "methyl_surro"

  result2 <- methyl_miss(single_col_surro)
  expect_equal(result2$summary$total_probes, 3)
  expect_equal(result2$summary$total_samples, 1)
  expect_true("cg2" %in% result2$missing_probes)  # completely missing in single column
})

test_that("methyl_miss() class assignment works correctly", {
  data(methyl_surro_miss)

  result <- methyl_miss(methyl_surro_miss)

  expect_s3_class(result, "methyl_miss")
  expect_true("methyl_miss" %in% class(result))
})

test_that("methyl_miss() handles matrices without row names", {
  # Create matrix without row names
  no_names_matrix <- matrix(c(0.1, NA, 0.3, 0.4, NA, 0.6), nrow = 2)
  colnames(no_names_matrix) <- c("s1", "s2", "s3")

  no_names_surro <- list(
    methyl = no_names_matrix,
    weights = c("1" = 1, "2" = 2),  # Use row indices as names
    intercept = 0
  )
  class(no_names_surro) <- "methyl_surro"

  result <- methyl_miss(no_names_surro)

  # Should still work, using row indices as names
  expect_s3_class(result, "methyl_miss")
  expect_true(length(result$missing_obs) > 0 || length(result$missing_probes) > 0)
})

test_that("methyl_miss() efficiently handles large matrices", {
  # Test with larger matrix to ensure efficiency
  large_matrix <- matrix(rnorm(1000), nrow = 100, ncol = 10,
                         dimnames = list(paste0("cg", 1:100), paste0("s", 1:10)))

  # Add some missing data patterns
  large_matrix[sample(1000, 100)] <- NA  # Random 10% missing
  large_matrix[91:100, ] <- NA  # Last 10 probes completely missing

  large_surro <- list(
    methyl = large_matrix,
    weights = setNames(rep(1, 100), paste0("cg", 1:100)),
    intercept = 0
  )
  class(large_surro) <- "methyl_surro"

  # Function should complete quickly
  start_time <- Sys.time()
  result <- methyl_miss(large_surro)
  end_time <- Sys.time()

  expect_s3_class(result, "methyl_miss")
  expect_true(length(result$missing_probes) == 10)  # Last 10 probes
  expect_true(as.numeric(end_time - start_time) < 1)  # Should be fast
})
