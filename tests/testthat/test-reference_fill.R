# Test file for reference_fill.R

test_that("reference_fill() handles edge cases with single row/column", {
  data(ref_df)

  ref_vec_mean <- setNames(ref_df$mean, rownames(ref_df))

  # Single row matrix with missing values
  single_row_matrix <- matrix(c(0.3, NA, 0.5), nrow = 1,
                              dimnames = list("cg02", c("s1", "s2", "s3")))
  single_row_surro <- list(
    methyl = single_row_matrix,
    weights = c(cg02 = 1.0),
    intercept = 0
  )
  class(single_row_surro) <- "methyl_surro"

  result <- reference_fill(single_row_surro, ref_vec_mean[c("cg02")], type = "obs")
  expect_equal(as.numeric(result$methyl[1, 2]), as.numeric(ref_vec_mean["cg02"]))

  # Single column matrix with completely missing probe
  single_col_matrix <- matrix(c(0.2, NA, 0.8), ncol = 1,
                              dimnames = list(c("cg02", "cg03", "cg07"), "s1"))
  single_col_surro <- list(
    methyl = single_col_matrix,
    weights = setNames(c(1, 2, 3), c("cg02", "cg03", "cg07")),
    intercept = 0
  )
  class(single_col_surro) <- "methyl_surro"

  result2 <- reference_fill(single_col_surro, ref_vec_mean[c("cg02", "cg03", "cg07")], type = "all")
  expect_equal(as.numeric(result2$methyl["cg03", "s1"]), as.numeric(ref_vec_mean["cg03"]))
})

test_that("reference_fill() numeric precision is maintained", {
  data(beta_matrix_miss)
  data(wts_df)

  wts_vec_lin <- setNames(wts_df$wt_lin, rownames(wts_df))

  methyl_surro_miss <- surro_set(
    methyl = beta_matrix_miss,
    weights = wts_vec_lin,
    intercept = "Intercept"
  )

  precision_ref <- setNames(c(0.123456789, 0.987654321), c("cg03", "cg06"))

  result <- reference_fill(methyl_surro_miss, precision_ref, type = "probes")

  expect_equal(as.numeric(result$methyl["cg03", 1]), as.numeric(precision_ref["cg03"]), tolerance = 1e-10)
  expect_equal(as.numeric(result$methyl["cg06", 1]), as.numeric(precision_ref["cg06"]), tolerance = 1e-10)
})


test_that("reference_fill() validates methyl_surro class", {
  data(ref_df)
  ref_vec <- setNames(ref_df$mean, rownames(ref_df))

  # Not a methyl_surro object
  fake_object <- list(methyl = matrix(1:9, nrow = 3))
  expect_error(
    reference_fill(fake_object, ref_vec, type = "probes"),
    "Input must be an object of class 'methyl_surro'"
  )
})

test_that("reference_fill() validates methyl matrix exists", {
  data(ref_df)
  ref_vec <- setNames(ref_df$mean, rownames(ref_df))

  # Missing methyl component
  bad_surro <- list(weights = c(1, 2, 3), intercept = 0)
  class(bad_surro) <- "methyl_surro"

  expect_error(
    reference_fill(bad_surro, ref_vec, type = "probes"),
    "methyl_surro object must contain a 'methyl' matrix component"
  )

  # Methyl is not a matrix
  bad_surro2 <- list(methyl = c(1, 2, 3), weights = c(1, 2, 3), intercept = 0)
  class(bad_surro2) <- "methyl_surro"

  expect_error(
    reference_fill(bad_surro2, ref_vec, type = "probes"),
    "methyl_surro object must contain a 'methyl' matrix component"
  )
})

test_that("reference_fill() validates empty matrix", {
  data(ref_df)
  ref_vec <- setNames(ref_df$mean, rownames(ref_df))

  # Empty matrix (0 rows)
  empty_matrix <- matrix(numeric(0), nrow = 0, ncol = 3,
                         dimnames = list(NULL, c("s1", "s2", "s3")))
  empty_surro <- list(methyl = empty_matrix, weights = numeric(0), intercept = 0)
  class(empty_surro) <- "methyl_surro"

  expect_error(
    reference_fill(empty_surro, ref_vec, type = "probes"),
    "Methylation matrix cannot be empty"
  )

  # Empty matrix (0 columns)
  empty_matrix2 <- matrix(numeric(0), nrow = 3, ncol = 0,
                          dimnames = list(c("cg01", "cg02", "cg03"), NULL))
  empty_surro2 <- list(methyl = empty_matrix2, weights = c(1, 2, 3), intercept = 0)
  class(empty_surro2) <- "methyl_surro"

  expect_error(
    reference_fill(empty_surro2, ref_vec, type = "probes"),
    "Methylation matrix cannot be empty"
  )
})

test_that("reference_fill() validates reference format", {
  data(beta_matrix_miss)
  data(wts_df)

  wts_vec <- setNames(wts_df$wt_lin, rownames(wts_df))
  methyl_surro <- surro_set(beta_matrix_miss, wts_vec, "Intercept")

  # Not numeric
  expect_error(
    reference_fill(methyl_surro, c("a", "b", "c"), type = "probes"),
    "Reference must be a named numeric vector"
  )

  # Not a vector (matrix instead)
  ref_matrix <- matrix(1:9, nrow = 3)
  expect_error(
    reference_fill(methyl_surro, ref_matrix, type = "probes"),
    "Reference must be a named numeric vector"
  )

  # No names
  expect_error(
    reference_fill(methyl_surro, c(0.5, 0.6, 0.7), type = "probes"),
    "Reference must be a named numeric vector"
  )
})

test_that("reference_fill() handles reference with NA values", {
  data(beta_matrix_miss)
  data(wts_df)
  data(ref_df)

  wts_vec <- setNames(wts_df$wt_lin, rownames(wts_df))
  methyl_surro <- surro_set(beta_matrix_miss, wts_vec, "Intercept")

  # Reference with some NAs
  ref_with_na <- setNames(ref_df$mean, rownames(ref_df))
  ref_with_na[1:2] <- NA

  expect_message(
    result <- reference_fill(methyl_surro, ref_with_na, type = "probes", verbose = TRUE),
    "Removed.*reference probes with NA values"
  )

  expect_s3_class(result, "methyl_surro")
})

test_that("reference_fill() stops when all reference values are NA", {
  data(beta_matrix_miss)
  data(wts_df)

  wts_vec <- setNames(wts_df$wt_lin, rownames(wts_df))
  methyl_surro <- surro_set(beta_matrix_miss, wts_vec, "Intercept")

  # All NAs
  ref_all_na <- setNames(rep(NA_real_, 5), paste0("cg", 1:5))

  expect_error(
    reference_fill(methyl_surro, ref_all_na, type = "probes"),
    "No valid reference values available after removing NAs"
  )
})

test_that("reference_fill() warns about out-of-range reference values", {
  data(beta_matrix_miss)
  data(wts_df)

  wts_vec <- setNames(wts_df$wt_lin, rownames(wts_df))
  methyl_surro <- surro_set(beta_matrix_miss, wts_vec, "Intercept")

  # Values way outside normal range
  extreme_ref <- setNames(c(-20, 25, 0.5), c("cg01", "cg02", "cg03"))

  expect_warning(
    reference_fill(methyl_surro, extreme_ref, type = "all"),
    "Reference values outside typical methylation range"
  )
})

test_that("reference_fill() type='probes' only fills completely missing", {
  # Create matrix with completely and partially missing probes
  test_matrix <- matrix(
    c(0.1, 0.2, 0.3,  # cg01 - complete
      NA,  NA,  NA,   # cg02 - completely missing
      0.4, NA,  0.5,  # cg03 - partially missing
      NA,  NA,  NA),  # cg04 - completely missing
    nrow = 4, byrow = TRUE,
    dimnames = list(paste0("cg0", 1:4), paste0("s", 1:3))
  )

  test_surro <- list(
    methyl = test_matrix,
    weights = setNames(1:4, paste0("cg0", 1:4)),
    intercept = 0
  )
  class(test_surro) <- "methyl_surro"

  ref_vec <- setNames(c(0.6, 0.7, 0.8, 0.9), paste0("cg0", 1:4))

  result <- reference_fill(test_surro, ref_vec, type = "probes")

  # cg02 and cg04 should be filled
  expect_equal(as.numeric(result$methyl["cg02", ]), rep(0.7, 3))
  expect_equal(as.numeric(result$methyl["cg04", ]), rep(0.9, 3))

  # cg03 should still have NA (partially missing)
  expect_true(is.na(result$methyl["cg03", "s2"]))
  expect_equal(result$methyl["cg03", "s1"], 0.4)
})

test_that("reference_fill() type='obs' only fills partially missing", {
  # Create matrix with completely and partially missing probes
  test_matrix <- matrix(
    c(0.1, 0.2, 0.3,  # cg01 - complete
      NA,  NA,  NA,   # cg02 - completely missing
      0.4, NA,  0.5,  # cg03 - partially missing
      0.6, 0.7, NA),  # cg04 - partially missing
    nrow = 4, byrow = TRUE,
    dimnames = list(paste0("cg0", 1:4), paste0("s", 1:3))
  )

  test_surro <- list(
    methyl = test_matrix,
    weights = setNames(1:4, paste0("cg0", 1:4)),
    intercept = 0
  )
  class(test_surro) <- "methyl_surro"

  ref_vec <- setNames(c(0.6, 0.7, 0.8, 0.9), paste0("cg0", 1:4))

  result <- reference_fill(test_surro, ref_vec, type = "obs")

  # cg03 and cg04 partially missing should be filled
  expect_equal(as.numeric(result$methyl["cg03", "s2"]), 0.8)
  expect_equal(as.numeric(result$methyl["cg04", "s3"]), 0.9)

  # cg02 completely missing should still be all NA
  expect_true(all(is.na(result$methyl["cg02", ])))
})

test_that("reference_fill() type='all' fills both types", {
  test_matrix <- matrix(
    c(0.1, 0.2, 0.3,  # cg01 - complete
      NA,  NA,  NA,   # cg02 - completely missing
      0.4, NA,  0.5), # cg03 - partially missing
    nrow = 3, byrow = TRUE,
    dimnames = list(paste0("cg0", 1:3), paste0("s", 1:3))
  )

  test_surro <- list(
    methyl = test_matrix,
    weights = setNames(1:3, paste0("cg0", 1:3)),
    intercept = 0
  )
  class(test_surro) <- "methyl_surro"

  ref_vec <- setNames(c(0.6, 0.7, 0.8), paste0("cg0", 1:3))

  result <- reference_fill(test_surro, ref_vec, type = "all")

  # Both types should be filled
  expect_equal(as.numeric(result$methyl["cg02", ]), rep(0.7, 3))
  expect_equal(as.numeric(result$methyl["cg03", "s2"]), 0.8)
})

test_that("reference_fill() handles matrix with no missing values", {
  complete_matrix <- matrix(
    runif(12, 0, 1), nrow = 4,
    dimnames = list(paste0("cg0", 1:4), paste0("s", 1:3))
  )

  complete_surro <- list(
    methyl = complete_matrix,
    weights = setNames(1:4, paste0("cg0", 1:4)),
    intercept = 0
  )
  class(complete_surro) <- "methyl_surro"

  ref_vec <- setNames(rep(0.5, 4), paste0("cg0", 1:4))

  expect_message(
    result <- reference_fill(complete_surro, ref_vec, type = "all", verbose = TRUE),
    "No missing values found"
  )

  # Matrix should be unchanged
  expect_equal(result$methyl, complete_matrix)
})

test_that("reference_fill() with no missing returns stats correctly", {
  complete_matrix <- matrix(
    runif(9, 0, 1), nrow = 3,
    dimnames = list(paste0("cg0", 1:3), paste0("s", 1:3))
  )

  complete_surro <- list(
    methyl = complete_matrix,
    weights = setNames(1:3, paste0("cg0", 1:3)),
    intercept = 0
  )
  class(complete_surro) <- "methyl_surro"

  ref_vec <- setNames(rep(0.5, 3), paste0("cg0", 1:3))

  result <- reference_fill(complete_surro, ref_vec, type = "all",
                           return_stats = TRUE, verbose = FALSE)

  expect_true("reference_fill_stats" %in% names(result))
  expect_equal(result$reference_fill_stats$n_total_values_filled, 0)
  expect_equal(result$reference_fill_stats$n_missing_before_filling, 0)
})

test_that("reference_fill() verbose output works correctly", {
  data(beta_matrix_miss)
  data(wts_df)
  data(ref_df)

  wts_vec <- setNames(wts_df$wt_lin, rownames(wts_df))
  methyl_surro <- surro_set(beta_matrix_miss, wts_vec, "Intercept")

  ref_vec <- setNames(ref_df$mean, rownames(ref_df))

  expect_message(
    reference_fill(methyl_surro, ref_vec, type = "all", verbose = TRUE),
    "Starting reference filling analysis"
  )

  expect_message(
    reference_fill(methyl_surro, ref_vec, type = "all", verbose = TRUE),
    "Matrix dimensions"
  )

  expect_message(
    reference_fill(methyl_surro, ref_vec, type = "all", verbose = TRUE),
    "Reference filling completed successfully"
  )
})

test_that("reference_fill() shows remaining missing probes message", {
  # Create a matrix where some probes will remain missing
  test_matrix <- matrix(
    c(0.1, 0.2, 0.3,
      NA,  NA,  NA,   # Will remain missing (no reference)
      NA,  NA,  NA),  # Will be filled
    nrow = 3, byrow = TRUE,
    dimnames = list(c("cg01", "cg02", "cg03"), paste0("s", 1:3))
  )

  test_surro <- list(
    methyl = test_matrix,
    weights = setNames(1:3, c("cg01", "cg02", "cg03")),
    intercept = 0
  )
  class(test_surro) <- "methyl_surro"

  # Reference only for cg03
  ref_vec <- setNames(0.5, "cg03")

  expect_message(
    reference_fill(test_surro, ref_vec, type = "probes", verbose = FALSE),
    "probes remain missing.*methyl_miss"
  )
})

test_that("reference_fill() return_stats=TRUE adds stats component", {
  data(beta_matrix_miss)
  data(wts_df)
  data(ref_df)

  wts_vec <- setNames(wts_df$wt_lin, rownames(wts_df))
  methyl_surro <- surro_set(beta_matrix_miss, wts_vec, "Intercept")

  ref_vec <- setNames(ref_df$mean, rownames(ref_df))

  result <- reference_fill(methyl_surro, ref_vec, type = "all", return_stats = TRUE)

  expect_true("reference_fill_stats" %in% names(result))
  expect_s3_class(result$reference_fill_stats, "reference_fill_stats")

  stats <- result$reference_fill_stats
  expect_true("type" %in% names(stats))
  expect_true("timestamp" %in% names(stats))
  expect_true("n_total_values_filled" %in% names(stats))
  expect_true("fill_rate" %in% names(stats))
})

test_that("reference_fill() return_stats=FALSE does not add stats", {
  data(beta_matrix_miss)
  data(wts_df)
  data(ref_df)

  wts_vec <- setNames(wts_df$wt_lin, rownames(wts_df))
  methyl_surro <- surro_set(beta_matrix_miss, wts_vec, "Intercept")

  ref_vec <- setNames(ref_df$mean, rownames(ref_df))

  result <- reference_fill(methyl_surro, ref_vec, type = "all", return_stats = FALSE)

  expect_false("reference_fill_stats" %in% names(result))
})

test_that("print.reference_fill_stats works correctly", {
  data(beta_matrix_miss)
  data(wts_df)
  data(ref_df)

  wts_vec <- setNames(wts_df$wt_lin, rownames(wts_df))
  methyl_surro <- surro_set(beta_matrix_miss, wts_vec, "Intercept")

  ref_vec <- setNames(ref_df$mean, rownames(ref_df))

  result <- reference_fill(methyl_surro, ref_vec, type = "all", return_stats = TRUE)

  # Capture print output
  output <- capture.output(print(result$reference_fill_stats))

  expect_true(any(grepl("Reference Filling Statistics", output)))
  expect_true(any(grepl("Matrix Summary", output)))
  expect_true(any(grepl("Filling Results", output)))
})

test_that("print.reference_fill_stats handles different scenarios", {
  data(beta_matrix_miss)
  data(wts_df)
  data(ref_df)

  wts_vec <- setNames(wts_df$wt_lin, rownames(wts_df))
  methyl_surro <- surro_set(beta_matrix_miss, wts_vec, "Intercept")

  ref_vec <- setNames(ref_df$mean, rownames(ref_df))

  # Test with type="probes"
  result_probes <- reference_fill(methyl_surro, ref_vec, type = "probes", return_stats = TRUE)
  output_probes <- capture.output(print(result_probes$reference_fill_stats))
  expect_true(any(grepl("Type: probes", output_probes)))

  # Test with type="obs"
  result_obs <- reference_fill(methyl_surro, ref_vec, type = "obs", return_stats = TRUE)
  output_obs <- capture.output(print(result_obs$reference_fill_stats))
  expect_true(any(grepl("Type: obs", output_obs)))
})

test_that("reference_fill() does not modify original object", {
  data(beta_matrix_miss)
  data(wts_df)
  data(ref_df)

  wts_vec <- setNames(wts_df$wt_lin, rownames(wts_df))
  methyl_surro_original <- surro_set(beta_matrix_miss, wts_vec, "Intercept")

  # Store original values
  original_matrix <- methyl_surro_original$methyl

  ref_vec <- setNames(ref_df$mean, rownames(ref_df))

  # Run reference_fill
  result <- reference_fill(methyl_surro_original, ref_vec, type = "all")

  # Original should be unchanged
  expect_equal(methyl_surro_original$methyl, original_matrix)

  # Result should be different
  expect_false(identical(result$methyl, original_matrix))
})
