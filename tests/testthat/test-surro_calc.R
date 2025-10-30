# Test file for surro_calc.R

test_that("surro_calc() function works with complete data and linear transformation", {
  # load sample data
  data(beta_matrix_comp)
  data(wts_df)
  data(ref_df)

  # create vectors from data frames
  wts_vec_lin <- setNames(wts_df$wt_lin, rownames(wts_df))
  ref_vec_mean <- setNames(ref_df$mean, rownames(ref_df))

  # generate methyl_surro object
  methyl_surro_comp_lin <- surro_set(methyl = beta_matrix_comp,
                                     weights = wts_vec_lin,
                                     intercept = "Intercept") |>
    reference_fill(reference = ref_vec_mean,
                   type = "probes")

  # generate expected values
  expected_vals <- as.numeric(
    t(methyl_surro_comp_lin$methyl) %*%
      methyl_surro_comp_lin$weights[rownames(methyl_surro_comp_lin$methyl)]
  ) |>
    `names<-`(paste0("samp", 1:5))
  expected_vals <- expected_vals + methyl_surro_comp_lin$intercept

  # run function
  function_vals <- surro_calc(methyl_surro = methyl_surro_comp_lin,
                              transform = "linear")

  # compare result
  expect_equal(function_vals, expected_vals)
})

test_that("surro_calc() function works with complete data and count transformation", {
  # load sample data
  data(beta_matrix_comp)
  data(wts_df)
  data(ref_df)

  # create vectors from data frames
  wts_vec_cnt <- setNames(wts_df$wt_cnt, rownames(wts_df))
  ref_vec_mean <- setNames(ref_df$mean, rownames(ref_df))

  # generate methyl_surro object
  methyl_surro_comp_cnt <- surro_set(methyl = beta_matrix_comp,
                                     weights = wts_vec_cnt,
                                     intercept = "Intercept") |>
    reference_fill(reference = ref_vec_mean,
                   type = "probes")

  # generate expected values
  expected_vals <- as.numeric(
    t(methyl_surro_comp_cnt$methyl) %*%
      methyl_surro_comp_cnt$weights[rownames(methyl_surro_comp_cnt$methyl)]
  ) |>
    `names<-`(paste0("samp", 1:5))
  expected_vals <- expected_vals + methyl_surro_comp_cnt$intercept
  expected_vals <- exp(expected_vals)

  # run function
  function_vals <- surro_calc(methyl_surro = methyl_surro_comp_cnt,
                              transform = "count")

  # compare result
  expect_equal(function_vals, expected_vals)
})

test_that("surro_calc() function works with complete data and probability transformation", {
  # load sample data
  data(beta_matrix_comp)
  data(wts_df)
  data(ref_df)

  # create vectors from data frames
  wts_vec_prb <- setNames(wts_df$wt_prb, rownames(wts_df))
  ref_vec_mean <- setNames(ref_df$mean, rownames(ref_df))

  # generate methyl_surro object
  methyl_surro_comp_prb <- surro_set(methyl = beta_matrix_comp,
                                     weights = wts_vec_prb,
                                     intercept = "Intercept") |>
    reference_fill(reference = ref_vec_mean,
                   type = "probes")

  # generate expected values
  expected_vals <- as.numeric(
    t(methyl_surro_comp_prb$methyl) %*%
      methyl_surro_comp_prb$weights[rownames(methyl_surro_comp_prb$methyl)]
  ) |>
    `names<-`(paste0("samp", 1:5))
  expected_vals <- expected_vals + methyl_surro_comp_prb$intercept
  expected_vals <- 1 / (1 + exp(-expected_vals))

  # run function
  function_vals <- surro_calc(methyl_surro = methyl_surro_comp_prb,
                              transform = "probability")

  # compare result
  expect_equal(function_vals, expected_vals)
})

test_that("surro_calc() function works with missing samples and linear transformation", {
  # load sample data
  data(beta_matrix_miss)
  data(wts_df)
  data(ref_df)

  # create vectors from data frames
  wts_vec_lin <- setNames(wts_df$wt_lin, rownames(wts_df))
  ref_vec_mean <- setNames(ref_df$mean, rownames(ref_df))

  # generate methyl_surro object
  methyl_surro_miss_lin <- surro_set(methyl = beta_matrix_miss,
                                     weights = wts_vec_lin,
                                     intercept = "Intercept") |>
    reference_fill(reference = ref_vec_mean,
                   type = "probes") |>
    impute_obs(method = "mean", min_nonmiss_prop = 0.5)

  # generate expected values
  methyl_sub <- methyl_surro_miss_lin$methyl[, -c(2, 3, 5)]
  expected_vals <- as.numeric(
    t(methyl_sub) %*%
      methyl_surro_miss_lin$weights[rownames(methyl_sub)]
  ) |>
    `names<-`(paste0("samp", c(1, 4)))
  expected_vals <- expected_vals + methyl_surro_miss_lin$intercept

  # run function
  function_vals <- surro_calc(methyl_surro = methyl_surro_miss_lin,
                              transform = "linear")

  # compare result
  expect_equal(function_vals, expected_vals)
})

test_that("surro_calc() function works with missing probes and linear transformation", {
  # load sample data
  data(beta_matrix_comp)
  data(wts_df)

  # create vectors from data frames
  wts_vec_lin <- setNames(wts_df$wt_lin, rownames(wts_df))

  # generate methyl_surro object
  methyl_surro_miss_lin <- surro_set(methyl = beta_matrix_comp,
                                     weights = wts_vec_lin,
                                     intercept = "Intercept")

  # generate expected values
  methyl_sub <- methyl_surro_miss_lin$methyl
  for (probe in rownames(methyl_sub)) {
    for (samp in 1:5) {
      if (is.na(methyl_sub[probe, samp])) {
        methyl_sub[probe, samp] <- 0
      }
    }
  }
  expected_vals <- as.numeric(
    t(methyl_sub) %*%
      methyl_surro_miss_lin$weights[rownames(methyl_sub)]
  ) |>
    `names<-`(paste0("samp", c(1:5)))
  expected_vals <- expected_vals + methyl_surro_miss_lin$intercept

  # run function
  function_vals <- surro_calc(methyl_surro = methyl_surro_miss_lin,
                              transform = "linear")

  # compare result
  expect_equal(function_vals, expected_vals)
})

test_that("surro_calc() function works with missing probes, missing samples, and linear transformation", {
  # load sample data
  data(beta_matrix_miss)
  data(wts_df)

  # create vectors from data frames
  wts_vec_lin <- setNames(wts_df$wt_lin, rownames(wts_df))

  # generate methyl_surro object
  methyl_surro_miss_lin <- surro_set(methyl = beta_matrix_miss,
                                     weights = wts_vec_lin,
                                     intercept = "Intercept") |>
    impute_obs(method = "mean", min_nonmiss_prop = 0.5)

  # generate expected values
  methyl_sub <- methyl_surro_miss_lin$methyl[, -c(2, 3, 5)]
  for (probe in rownames(methyl_sub)) {
    for (samp in 1:2) {
      if (is.na(methyl_sub[probe, samp])) {
        methyl_sub[probe, samp] <- 0
      }
    }
  }
  expected_vals <- as.numeric(
    t(methyl_sub) %*%
      methyl_surro_miss_lin$weights[rownames(methyl_sub)]
  ) |>
    `names<-`(paste0("samp", c(1, 4)))
  expected_vals <- expected_vals + methyl_surro_miss_lin$intercept

  # run function
  function_vals <- surro_calc(methyl_surro = methyl_surro_miss_lin,
                              transform = "linear")

  # compare result
  expect_equal(function_vals, expected_vals)
})

test_that("surro_calc() input validation works", {
  data(beta_matrix_comp)
  data(wts_df)

  # create vectors from data frames
  wts_vec_lin <- setNames(wts_df$wt_lin, rownames(wts_df))

  # Test non-methyl_surro object
  expect_error(surro_calc(list(methyl = matrix(1:4, 2, 2))),
               "Input must be an object of class 'methyl_surro'")

  # Test missing methyl component
  bad_surro <- list(weights = wts_vec_lin, intercept = 0)
  class(bad_surro) <- "methyl_surro"
  expect_error(surro_calc(bad_surro), "must contain a 'methyl' matrix component")

  # Test non-matrix methyl component
  bad_surro2 <- list(methyl = c(1, 2, 3), weights = wts_vec_lin, intercept = 0)
  class(bad_surro2) <- "methyl_surro"
  expect_error(surro_calc(bad_surro2), "must contain a 'methyl' matrix component")

  # Test missing weights component
  bad_surro3 <- list(methyl = beta_matrix_comp, intercept = 0)
  class(bad_surro3) <- "methyl_surro"
  expect_error(surro_calc(bad_surro3), "must contain a numeric 'weights' component")

  # Test non-numeric weights
  bad_surro4 <- list(methyl = beta_matrix_comp, weights = c("a", "b"), intercept = 0)
  class(bad_surro4) <- "methyl_surro"
  expect_error(surro_calc(bad_surro4), "must contain a numeric 'weights' component")

  # Test empty matrix
  empty_surro <- list(
    methyl = matrix(numeric(0), nrow = 0, ncol = 0),
    weights = c(cg1 = 1.0),
    intercept = 0.5
  )
  class(empty_surro) <- "methyl_surro"
  expect_error(surro_calc(empty_surro), "cannot be empty")

  # Test empty weights
  empty_weights_surro <- list(
    methyl = beta_matrix_comp,
    weights = numeric(0),
    intercept = 0
  )
  class(empty_weights_surro) <- "methyl_surro"
  expect_error(surro_calc(empty_weights_surro), "cannot be empty")

  # Test invalid transform
  valid_surro <- surro_set(beta_matrix_comp, wts_vec_lin, "Intercept")
  expect_error(surro_calc(valid_surro, transform = "invalid"),
               "'arg' should be one of")
})

test_that("surro_calc() probe-weight alignment validation works", {
  data(beta_matrix_comp)

  # Test weights without names
  unnamed_weights <- c(1, 2, 3, 4, 5)
  unnamed_surro <- list(methyl = beta_matrix_comp, weights = unnamed_weights, intercept = 0)
  class(unnamed_surro) <- "methyl_surro"
  expect_error(surro_calc(unnamed_surro), "must be a named vector")

  # Test methyl matrix without row names
  no_names_matrix <- beta_matrix_comp
  rownames(no_names_matrix) <- NULL
  no_names_surro <- list(
    methyl = no_names_matrix,
    weights = setNames(1:15, paste0("cg", 1:15)),
    intercept = 0
  )
  class(no_names_surro) <- "methyl_surro"
  expect_error(surro_calc(no_names_surro), "must have row names")

  # Test no common probes between weights and methylation
  no_common_weights <- setNames(c(1, 2, 3), c("probe1", "probe2", "probe3"))
  no_common_surro <- list(methyl = beta_matrix_comp, weights = no_common_weights, intercept = 0)
  class(no_common_surro) <- "methyl_surro"
  expect_error(surro_calc(no_common_surro), "No common probes found")
})

test_that("surro_calc() verbose messaging works", {
  data(beta_matrix_comp)
  data(wts_df)

  # create vectors from data frames
  wts_vec_lin <- setNames(wts_df$wt_lin, rownames(wts_df))

  surro <- surro_set(beta_matrix_comp, wts_vec_lin, "Intercept")

  expect_message(surro_calc(surro, verbose = TRUE), "Starting surrogate calculation")
  expect_message(surro_calc(surro, verbose = TRUE), "Input:")
  expect_message(surro_calc(surro, verbose = TRUE), "Weights:")
  expect_message(surro_calc(surro, verbose = TRUE), "Calculation completed successfully")

  # Test non-verbose mode - capture messages and check they're minimal
  messages <- capture.output({
    result <- surro_calc(surro, verbose = FALSE)
  }, type = "message")

  # Should only have essential warnings about missing probes, not verbose output
  expect_false(any(grepl("Starting surrogate calculation", messages)))
  expect_false(any(grepl("Calculation completed", messages)))
})

test_that("surro_calc() return_diagnostics functionality works", {
  data(beta_matrix_comp)
  data(wts_df)

  # create vectors from data frames
  wts_vec_lin <- setNames(wts_df$wt_lin, rownames(wts_df))

  surro <- surro_set(beta_matrix_comp, wts_vec_lin, "Intercept")
  result <- surro_calc(surro, return_diagnostics = TRUE)

  expect_type(result, "list")
  expect_true(all(c("predictions", "diagnostics") %in% names(result)))

  # Test predictions component
  expect_type(result$predictions, "double")
  expect_true(!is.null(names(result$predictions)))

  # Test diagnostics component
  diag <- result$diagnostics
  required_diag_fields <- c("transform", "n_samples_original", "n_probes_original",
                            "n_weights", "intercept_used", "n_common_probes",
                            "n_completely_missing_probes", "n_samples_with_missing",
                            "n_final_samples", "final_result_range")

  expect_true(all(required_diag_fields %in% names(diag)))
  expect_equal(diag$transform, "linear")
  expect_true(diag$intercept_used)
})

test_that("surro_calc() handles intercept correctly", {
  data(beta_matrix_comp)
  data(wts_df)

  # create vectors from data frames
  wts_vec_lin <- setNames(wts_df$wt_lin, rownames(wts_df))

  # Test with intercept
  surro_with_int <- surro_set(beta_matrix_comp, wts_vec_lin, "Intercept")
  result_with_int <- surro_calc(surro_with_int, verbose = TRUE)

  expect_message(surro_calc(surro_with_int, verbose = TRUE), "Intercept:")

  # Test without intercept
  weights_no_int <- wts_vec_lin[names(wts_vec_lin) != "Intercept"]
  surro_no_int <- surro_set(beta_matrix_comp, weights_no_int)
  result_no_int <- surro_calc(surro_no_int, return_diagnostics = TRUE)

  expect_false(result_no_int$diagnostics$intercept_used)

  # Results should be different
  expect_false(all(result_with_int == result_no_int$predictions))
})

test_that("surro_calc() handles missing data patterns correctly", {
  data(wts_df)

  # create vectors from data frames
  wts_vec_lin <- setNames(wts_df$wt_lin, rownames(wts_df))

  # Create matrix with specific missing pattern
  mixed_missing_matrix <- matrix(c(
    0.1, 0.2, 0.3, 0.4, 0.5,  # complete (samp1-5)
    NA, NA, NA, NA, NA,       # completely missing
    0.2, NA, 0.4, 0.5, 0.6,   # partially missing (missing samp2)
    0.3, 0.4, 0.5, 0.6, 0.7   # complete
  ), nrow = 4, byrow = TRUE,
  dimnames = list(c("cg02", "cg03", "cg07", "cg08"), paste0("samp", 1:5)))

  mixed_surro <- list(
    methyl = mixed_missing_matrix,
    weights = wts_vec_lin[c("cg02", "cg03", "cg07", "cg08")],
    intercept = wts_vec_lin["Intercept"]
  )
  class(mixed_surro) <- "methyl_surro"

  expect_message(
    result <- surro_calc(mixed_surro, verbose = FALSE),
    "probes set to zero.*samples omitted"
  )

  # Samples with ANY missing values are removed
  # samp2 has missing in cg07, so should be removed
  # Only samp1, samp3, samp4, samp5 should remain (4 samples)
  expect_length(result, 4)
  expect_true(all(names(result) %in% c("samp1", "samp3", "samp4", "samp5")))
})

test_that("surro_calc() handles edge cases with single row/column", {
  data(wts_df)

  # create vectors from data frames
  wts_vec_lin <- setNames(wts_df$wt_lin, rownames(wts_df))

  # Single probe
  single_probe_matrix <- matrix(c(0.3, 0.4, 0.5), nrow = 1,
                                dimnames = list("cg02", c("s1", "s2", "s3")))
  single_probe_surro <- list(
    methyl = single_probe_matrix,
    weights = wts_vec_lin["cg02"],
    intercept = wts_vec_lin["Intercept"]
  )
  class(single_probe_surro) <- "methyl_surro"

  result1 <- surro_calc(single_probe_surro)
  expect_length(result1, 3)
  expect_true(all(names(result1) == c("s1", "s2", "s3")))

  # Single sample
  single_sample_matrix <- matrix(c(0.2, 0.3, 0.4), ncol = 1,
                                 dimnames = list(c("cg02", "cg07", "cg08"), "s1"))
  single_sample_surro <- list(
    methyl = single_sample_matrix,
    weights = wts_vec_lin[c("cg02", "cg07", "cg08")],
    intercept = wts_vec_lin["Intercept"]
  )
  class(single_sample_surro) <- "methyl_surro"

  result2 <- surro_calc(single_sample_surro)
  expect_length(result2, 1)
  expect_equal(names(result2), "s1")
})

test_that("surro_calc() transformation functions work correctly", {
  # Create a simple test case
  simple_matrix <- matrix(c(0.3, 0.4, 0.5, 0.6), nrow = 2,
                          dimnames = list(c("cg02", "cg07"), c("s1", "s2")))
  simple_surro <- list(
    methyl = simple_matrix,
    weights = c(cg02 = 1.0, cg07 = 2.0),
    intercept = 0.5
  )
  class(simple_surro) <- "methyl_surro"

  # Test linear transformation - get actual results and check they make sense
  linear_result <- surro_calc(simple_surro, transform = "linear", verbose = FALSE)

  # Verify we get sensible numeric results
  expect_type(linear_result, "double")
  expect_length(linear_result, 2)
  expect_true(all(linear_result > 0))  # Should be positive given our test data

  # Test count transformation (exp of linear)
  count_result <- surro_calc(simple_surro, transform = "count", verbose = FALSE)
  expect_true(all(count_result > linear_result))  # exp should be larger

  # Test probability transformation (sigmoid of linear)
  prob_result <- surro_calc(simple_surro, transform = "probability", verbose = FALSE)
  expect_true(all(prob_result > 0 & prob_result < 1))  # Should be between 0 and 1
})

test_that("surro_calc() handles warnings and recommendations correctly", {
  data(wts_df)

  # create vectors from data frames
  wts_vec_lin <- setNames(wts_df$wt_lin, rownames(wts_df))

  # Test case with completely missing probes
  missing_probe_matrix <- matrix(c(
    0.1, 0.2, 0.3,
    NA, NA, NA,
    0.4, 0.5, 0.6
  ), nrow = 3, byrow = TRUE,
  dimnames = list(c("cg02", "cg03", "cg07"), c("s1", "s2", "s3")))

  missing_probe_surro <- list(
    methyl = missing_probe_matrix,
    weights = wts_vec_lin[c("cg02", "cg03", "cg07")],
    intercept = wts_vec_lin["Intercept"]
  )
  class(missing_probe_surro) <- "methyl_surro"

  expect_message(
    surro_calc(missing_probe_surro),
    "probes set to zero.*reference_fill"
  )

  # Test case with missing samples
  missing_sample_matrix <- matrix(c(
    0.1, NA, 0.3,
    0.2, 0.4, 0.5,
    0.3, 0.6, 0.7
  ), nrow = 3, byrow = TRUE,
  dimnames = list(c("cg02", "cg03", "cg07"), c("s1", "s2", "s3")))

  missing_sample_surro <- list(
    methyl = missing_sample_matrix,
    weights = wts_vec_lin[c("cg02", "cg03", "cg07")],
    intercept = wts_vec_lin["Intercept"]
  )
  class(missing_sample_surro) <- "methyl_surro"

  expect_message(
    surro_calc(missing_sample_surro),
    "samples omitted.*impute_obs"
  )
})

test_that("surro_calc() probe filtering and alignment works", {
  data(beta_matrix_comp)
  data(wts_df)

  # create vectors from data frames
  wts_vec_lin <- setNames(wts_df$wt_lin, rownames(wts_df))

  # Add extra probes to methylation matrix that aren't in weights
  extended_matrix <- rbind(beta_matrix_comp,
                           matrix(runif(10), nrow = 2,
                                  dimnames = list(c("extra1", "extra2"), NULL)))

  # Add extra weights that aren't in methylation matrix
  extended_weights <- c(wts_vec_lin, c(missing1 = 1.5, missing2 = 2.5))

  extended_surro <- list(
    methyl = extended_matrix,
    weights = extended_weights,
    intercept = NULL
  )
  class(extended_surro) <- "methyl_surro"

  expect_message(
    result <- surro_calc(extended_surro, verbose = TRUE),
    "probes in methylation data are not in weights"
  )

  expect_message(
    surro_calc(extended_surro, verbose = TRUE),
    "weight probes are not in methylation data"
  )

  # Should still work with common probes
  expect_type(result, "double")
  expect_true(length(result) > 0)
})

test_that("surro_calc() diagnostics contain correct information", {
  data(beta_matrix_comp)
  data(wts_df)

  # create vectors from data frames
  wts_vec_lin <- setNames(wts_df$wt_lin, rownames(wts_df))

  surro <- surro_set(beta_matrix_comp, wts_vec_lin, "Intercept")
  result <- surro_calc(surro, transform = "probability", return_diagnostics = TRUE)

  diag <- result$diagnostics

  # Test specific diagnostic values - use the surro object's dimensions
  expect_equal(diag$n_samples_original, ncol(beta_matrix_comp))
  expect_equal(diag$n_probes_original, nrow(surro$methyl))
  expect_equal(diag$n_weights, length(wts_vec_lin) - 1)  # Excluding intercept
  expect_equal(diag$transform, "probability")
  expect_true(diag$intercept_used)
  expect_equal(diag$n_final_samples, ncol(beta_matrix_comp))

  # Test result range for probability transformation
  expect_true(diag$final_result_range["min"] >= 0)
  expect_true(diag$final_result_range["max"] <= 1)
})
