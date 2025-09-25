# Test file for reference_fill.R

test_that("reference_fill() function correctly handles 'probes' type with reference vector", {
  # load sample data
  data(methyl_surro_miss)
  data(ref_vec_mean)

  # generate expected values
  expected_vals <- methyl_surro_miss
  for (probe in c("cg03", "cg06", "cg11", "cg15", "cg18")) {
    for (samp in 1:5) {
      expected_vals$methyl[probe, samp] <- ref_vec_mean[probe]
    }
  }

  # run function
  function_vals <- reference_fill(methyl_surro = methyl_surro_miss,
                                  reference = ref_vec_mean,
                                  type = "probes")

  # compare result
  expect_equal(function_vals, expected_vals)
})

test_that("reference_fill() function correctly handles 'obs' type with reference vector", {
  # load sample data
  data(methyl_surro_miss)
  data(ref_vec_mean)

  # generate expected values
  expected_vals <- methyl_surro_miss
  for (probe in c("cg02", "cg07")) {
    for (samp in 1:5) {
      if (is.na(expected_vals$methyl[probe, samp])) {
        expected_vals$methyl[probe, samp] <- ref_vec_mean[probe]
      }
    }
  }

  # run function
  function_vals <- reference_fill(methyl_surro = methyl_surro_miss,
                                  reference = ref_vec_mean,
                                  type = "obs")
  # compare results
  expect_equal(function_vals, expected_vals)
})

test_that("reference_fill() function correctly handles 'all' type with reference vector", {
  # load sample data
  data(methyl_surro_miss)
  data(ref_vec_mean)

  # generate expected values
  expected_vals <- methyl_surro_miss
  for (probe in rownames(methyl_surro_miss$methyl)) {
    for (samp in 1:5) {
      if (is.na(expected_vals$methyl[probe, samp])) {
        expected_vals$methyl[probe, samp] <- ref_vec_mean[probe]
      }
    }
  }

  # run function
  function_vals <- reference_fill(methyl_surro = methyl_surro_miss,
                                  reference = ref_vec_mean,
                                  type = "all")
  # compare results
  expect_equal(function_vals, expected_vals)
})

# NEW TESTS FOR COMPLETE COVERAGE

test_that("reference_fill() input validation works", {
  data(methyl_surro_miss)
  data(ref_vec_mean)

  # Test non-methyl_surro object
  expect_error(reference_fill(list(methyl = matrix(1:4, 2, 2)), ref_vec_mean),
               "Input must be an object of class 'methyl_surro'")

  # Test missing methyl component
  bad_surro <- list(weights = c(1, 2), intercept = 0)
  class(bad_surro) <- "methyl_surro"
  expect_error(reference_fill(bad_surro, ref_vec_mean),
               "must contain a 'methyl' matrix component")

  # Test non-matrix methyl component
  bad_surro2 <- list(methyl = c(1, 2, 3), weights = c(1, 2), intercept = 0)
  class(bad_surro2) <- "methyl_surro"
  expect_error(reference_fill(bad_surro2, ref_vec_mean),
               "must contain a 'methyl' matrix component")

  # Test empty matrix
  empty_surro <- list(
    methyl = matrix(numeric(0), nrow = 0, ncol = 0),
    weights = c(cg1 = 1.0),
    intercept = 0.5
  )
  class(empty_surro) <- "methyl_surro"
  expect_error(reference_fill(empty_surro, ref_vec_mean), "cannot be empty")

  # Test invalid reference data
  expect_error(reference_fill(methyl_surro_miss, c(1, 2, 3)),
               "Reference must be a named numeric vector")

  expect_error(reference_fill(methyl_surro_miss, matrix(1:4, 2, 2)),
               "Reference must be a named numeric vector")

  # Test reference with no names
  unnamed_ref <- c(0.5, 0.6, 0.7)
  expect_error(reference_fill(methyl_surro_miss, unnamed_ref),
               "Reference must be a named numeric vector")

  # Test invalid type
  expect_error(reference_fill(methyl_surro_miss, ref_vec_mean, type = "invalid"),
               "'arg' should be one of")
})

test_that("reference_fill() handles edge cases correctly", {
  data(methyl_surro_miss)
  data(ref_vec_mean)

  # Test with no missing values
  data(beta_matrix_comp)
  complete_surro <- list(
    methyl = beta_matrix_comp,
    weights = setNames(rep(1, nrow(beta_matrix_comp)), rownames(beta_matrix_comp)),
    intercept = 0
  )
  class(complete_surro) <- "methyl_surro"

  expect_message(result <- reference_fill(complete_surro, ref_vec_mean, verbose = TRUE),
                 "No missing values found")
  expect_identical(result$methyl, complete_surro$methyl)

  # Test with reference containing NA values
  ref_with_na <- ref_vec_mean
  ref_with_na[1:3] <- NA

  expect_message(reference_fill(methyl_surro_miss, ref_with_na, verbose = TRUE),
                 "Removed.*reference probes with NA values")
})

test_that("reference_fill() validates reference value ranges", {
  data(methyl_surro_miss)

  # Test with extreme values (should warn)
  extreme_ref <- setNames(c(-20, 25, 0.5), c("cg02", "cg03", "cg07"))
  expect_warning(reference_fill(methyl_surro_miss, extreme_ref),
                 "outside typical methylation range")

  # Test with beta-like values (should message if verbose)
  beta_ref <- setNames(c(0.2, 0.8, 0.5), c("cg02", "cg03", "cg07"))
  expect_message(reference_fill(methyl_surro_miss, beta_ref, verbose = TRUE),
                 "appear to be beta values")

  # Test with M-like values (should message if verbose)
  m_ref <- setNames(c(-2.1, 1.5, 0.2), c("cg02", "cg03", "cg07"))
  expect_message(reference_fill(methyl_surro_miss, m_ref, verbose = TRUE),
                 "appear to be M-values")
})

test_that("reference_fill() verbose messaging works", {
  data(methyl_surro_miss)
  data(ref_vec_mean)

  expect_message(
    reference_fill(methyl_surro_miss, ref_vec_mean, type = "all", verbose = TRUE),
    "Starting reference filling analysis"
  )

  expect_message(
    reference_fill(methyl_surro_miss, ref_vec_mean, type = "probes", verbose = TRUE),
    "Filling.*completely missing probes"
  )

  # Test non-verbose mode
  expect_message(
    reference_fill(methyl_surro_miss, ref_vec_mean, type = "all", verbose = FALSE),
    "Filled.*values using reference data"
  )
})

test_that("reference_fill() return_stats functionality works", {
  data(methyl_surro_miss)
  data(ref_vec_mean)

  result <- reference_fill(methyl_surro_miss, ref_vec_mean,
                           type = "all", return_stats = TRUE)

  expect_true("reference_fill_stats" %in% names(result))
  expect_s3_class(result$reference_fill_stats, "reference_fill_stats")

  stats <- result$reference_fill_stats
  required_fields <- c("type", "timestamp", "n_total_probes", "n_samples",
                       "n_reference_probes", "n_probes_filled_complete",
                       "n_probes_filled_partial", "n_values_filled_complete",
                       "n_values_filled_partial", "n_total_values_filled",
                       "n_missing_before_filling", "n_missing_after_filling",
                       "fill_rate", "probes_filled_complete", "probes_filled_partial",
                       "reference_range")

  expect_true(all(required_fields %in% names(stats)))
  expect_equal(stats$type, "all")
  expect_true(stats$fill_rate >= 0 && stats$fill_rate <= 1)
})

test_that("reference_fill_stats print method works", {
  data(methyl_surro_miss)
  data(ref_vec_mean)

  result <- reference_fill(methyl_surro_miss, ref_vec_mean,
                           type = "probes", return_stats = TRUE)

  expect_output(print(result$reference_fill_stats), "Reference Filling Statistics")
  expect_output(print(result$reference_fill_stats), "Type: probes")
  expect_output(print(result$reference_fill_stats), "Matrix Summary:")
  expect_output(print(result$reference_fill_stats), "Reference Data:")
  expect_output(print(result$reference_fill_stats), "Filling Results:")
})

test_that("reference_fill() preserves object structure", {
  data(methyl_surro_miss)
  data(ref_vec_mean)

  result <- reference_fill(methyl_surro_miss, ref_vec_mean, type = "all")
