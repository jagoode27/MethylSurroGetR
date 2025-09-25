# Test file for conversion_functions.R

test_that("convert_beta_to_m() function works on complete data", {
  # load sample beta matrix
  data(beta_matrix_comp)

  # generate expected values
  expected_mvals_comp <- log2(beta_matrix_comp / (1 - beta_matrix_comp))

  # run conversion function
  function_mvals_comp <- convert_beta_to_m(beta_matrix_comp)

  # compare result
  expect_equal(function_mvals_comp, expected_mvals_comp, tolerance = 0.0001)
})

test_that("convert_beta_to_m() function works on missing data", {
  # load sample beta matrix
  data(beta_matrix_miss)

  # generate expected values
  expected_mvals_miss <- log2(beta_matrix_miss / (1 - beta_matrix_miss))

  # run conversion function
  function_mvals_miss <- convert_beta_to_m(beta_matrix_miss)

  # compare result
  expect_equal(function_mvals_miss, expected_mvals_miss, tolerance = 0.0001)
})

test_that("convert_m_to_beta() function works on complete data", {
  # load sample M-value matrix
  data(mval_matrix_comp)

  # generate expected values
  expected_betas_comp <- 2^mval_matrix_comp / (2^mval_matrix_comp + 1)

  # run conversion function
  function_betas_comp <- convert_m_to_beta(mval_matrix_comp)

  # compare result
  expect_equal(function_betas_comp, expected_betas_comp, tolerance = 0.0001)
})

test_that("convert_m_to_beta() function works on missing data", {
  # load sample M-value matrix
  data(mval_matrix_miss)

  # generate expected values
  expected_betas_miss <- 2^mval_matrix_miss / (2^mval_matrix_miss + 1)

  # run conversion function
  function_betas_miss <- convert_m_to_beta(mval_matrix_miss)

  # compare result
  expect_equal(function_betas_miss, expected_betas_miss, tolerance = 0.0001)
})

# NEW TESTS FOR COMPLETE COVERAGE

test_that("convert_beta_to_m() handles edge cases correctly", {
  # Test extreme beta values that would cause infinite M-values
  extreme_beta <- matrix(c(0, 1, 0.5, 0.001, 0.999), nrow = 5, ncol = 1,
                         dimnames = list(paste0("cg", 1:5), "samp1"))

  result <- convert_beta_to_m(extreme_beta)

  # Check that no infinite values are produced
  expect_true(all(is.finite(result[!is.na(result)])))

  # Check that extreme values are clamped appropriately
  expect_true(all(result[!is.na(result)] > -50))  # reasonable M-value range
  expect_true(all(result[!is.na(result)] < 50))
})

test_that("convert_beta_to_m() with in_place = TRUE modifies original matrix", {
  data(beta_matrix_comp)
  original_copy <- beta_matrix_comp

  # Test in_place modification
  expect_message(
    result <- convert_beta_to_m(beta_matrix_comp, in_place = TRUE),
    "permanently modify"
  )

  # Original matrix should be modified
  expect_false(identical(beta_matrix_comp, original_copy))

  # Result should be the same as modified original
  expect_identical(result, beta_matrix_comp)
})

test_that("convert_m_to_beta() with in_place = TRUE modifies original matrix", {
  data(mval_matrix_comp)
  original_copy <- mval_matrix_comp

  # Test in_place modification
  expect_message(
    result <- convert_m_to_beta(mval_matrix_comp, in_place = TRUE),
    "permanently modify"
  )

  # Original matrix should be modified
  expect_false(identical(mval_matrix_comp, original_copy))

  # Result should be the same as modified original
  expect_identical(result, mval_matrix_comp)
})

test_that("convert_beta_to_m() input validation works", {
  # Test non-matrix input
  expect_error(convert_beta_to_m(c(0.5, 0.6, 0.7)), "Input must be a matrix")

  # Test non-numeric matrix
  char_matrix <- matrix(c("a", "b", "c", "d"), nrow = 2)
  expect_error(convert_beta_to_m(char_matrix), "must be numeric")

  # Test empty matrix
  empty_matrix <- matrix(numeric(0), nrow = 0, ncol = 0)
  expect_error(convert_beta_to_m(empty_matrix), "cannot be empty")

  # Test matrix with zero rows
  zero_row_matrix <- matrix(numeric(0), nrow = 0, ncol = 2)
  expect_error(convert_beta_to_m(zero_row_matrix), "cannot be empty")

  # Test matrix with zero columns
  zero_col_matrix <- matrix(numeric(0), nrow = 2, ncol = 0)
  expect_error(convert_beta_to_m(zero_col_matrix), "cannot be empty")
})

test_that("convert_m_to_beta() input validation works", {
  # Test non-matrix input
  expect_error(convert_m_to_beta(c(0.5, 0.6, 0.7)), "Input must be a matrix")

  # Test non-numeric matrix
  char_matrix <- matrix(c("a", "b", "c", "d"), nrow = 2)
  expect_error(convert_m_to_beta(char_matrix), "must be numeric")

  # Test empty matrix
  empty_matrix <- matrix(numeric(0), nrow = 0, ncol = 0)
  expect_error(convert_m_to_beta(empty_matrix), "cannot be empty")
})

test_that("convert_beta_to_m() warnings work correctly", {
  # Test missing row names warning
  no_names_matrix <- matrix(c(0.3, 0.4, 0.5, 0.6), nrow = 2)
  expect_warning(convert_beta_to_m(no_names_matrix), "lacks row names")

  # Test values outside [0,1] range warning
  out_of_range_matrix <- matrix(c(-0.1, 0.5, 1.2, 0.8), nrow = 2,
                                dimnames = list(c("cg1", "cg2"), c("s1", "s2")))
  expect_warning(convert_beta_to_m(out_of_range_matrix), "outside \\[0,1\\] range")
})

test_that("convert_m_to_beta() warnings work correctly", {
  # Test missing row names warning
  no_names_matrix <- matrix(c(1.5, -1.2, 0.8, -0.3), nrow = 2)
  expect_warning(convert_m_to_beta(no_names_matrix), "lacks row names")

  # Test out of range beta values after conversion
  extreme_m_matrix <- matrix(c(-20, 20, 0, 5), nrow = 2,
                             dimnames = list(c("cg1", "cg2"), c("s1", "s2")))
  expect_warning(convert_m_to_beta(extreme_m_matrix), "outside the expected \\[0,1\\] range")
})

test_that("convert_beta_to_m() preserves NA values correctly", {
  # Create matrix with NAs
  beta_with_na <- matrix(c(0.3, NA, 0.5, 0.6, NA, 0.8), nrow = 3,
                         dimnames = list(paste0("cg", 1:3), c("s1", "s2")))

  result <- convert_beta_to_m(beta_with_na)

  # Check that NAs are preserved in same positions
  expect_identical(is.na(beta_with_na), is.na(result))

  # Check that non-NA values are converted correctly
  non_na_positions <- !is.na(beta_with_na)
  expected_non_na <- log2(beta_with_na[non_na_positions] / (1 - beta_with_na[non_na_positions]))
  expect_equal(result[non_na_positions], expected_non_na, tolerance = 1e-10)
})

test_that("convert_m_to_beta() preserves NA values correctly", {
  # Create matrix with NAs
  m_with_na <- matrix(c(-1.5, NA, 0.8, 1.2, NA, -0.3), nrow = 3,
                      dimnames = list(paste0("cg", 1:3), c("s1", "s2")))

  result <- convert_m_to_beta(m_with_na)

  # Check that NAs are preserved in same positions
  expect_identical(is.na(m_with_na), is.na(result))

  # Check that non-NA values are converted correctly
  non_na_positions <- !is.na(m_with_na)
  expected_non_na <- 2^m_with_na[non_na_positions] / (2^m_with_na[non_na_positions] + 1)
  expect_equal(result[non_na_positions], expected_non_na, tolerance = 1e-10)
})

test_that("conversion functions are reciprocal", {
  data(beta_matrix_comp)

  # Convert beta to M and back to beta
  m_values <- convert_beta_to_m(beta_matrix_comp)
  beta_recovered <- convert_m_to_beta(m_values)

  expect_equal(beta_matrix_comp, beta_recovered, tolerance = 1e-10)

  # Convert M to beta and back to M
  data(mval_matrix_comp)
  beta_values <- convert_m_to_beta(mval_matrix_comp)
  m_recovered <- convert_beta_to_m(beta_values)

  expect_equal(mval_matrix_comp, m_recovered, tolerance = 1e-10)
})

test_that("convert_beta_to_m() numerical stability with in_place", {
  # Test that in_place conversion handles extreme values properly
  extreme_beta <- matrix(c(0, 1, 1e-10, 1-1e-10), nrow = 2, ncol = 2,
                         dimnames = list(c("cg1", "cg2"), c("s1", "s2")))

  result <- convert_beta_to_m(extreme_beta, in_place = TRUE)

  # Should not produce infinite values due to clamping
  expect_true(all(is.finite(result)))
  expect_true(all(abs(result) < 50))  # reasonable M-value range
})

test_that("functions handle single row/column matrices", {
  # Single row matrix
  single_row <- matrix(c(0.3, 0.4, 0.5), nrow = 1,
                       dimnames = list("cg1", c("s1", "s2", "s3")))

  m_result <- convert_beta_to_m(single_row)
  expect_equal(dim(m_result), c(1, 3))

  beta_result <- convert_m_to_beta(m_result)
  expect_equal(single_row, beta_result, tolerance = 1e-10)

  # Single column matrix
  single_col <- matrix(c(0.2, 0.6, 0.8), ncol = 1,
                       dimnames = list(c("cg1", "cg2", "cg3"), "s1"))

  m_result2 <- convert_beta_to_m(single_col)
  expect_equal(dim(m_result2), c(3, 1))

  beta_result2 <- convert_m_to_beta(m_result2)
  expect_equal(single_col, beta_result2, tolerance = 1e-10)
})
