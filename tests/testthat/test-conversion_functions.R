test_that("convert_beta_to_m() function works on complete data", {

  # load sample beta matrix
  data(beta_matrix_comp)

  # generate expected values
  expected_mvals_comp <- log2(beta_matrix_comp /
                           (1 - beta_matrix_comp))

  # run conversion function
  function_mvals_comp <- convert_beta_to_m(beta_matrix_comp)

  # compare result
  expect_equal(function_mvals_comp, expected_mvals_comp, tolerance = 0.0001)
})

test_that("convert_beta_to_m() function works on missing data", {

  # load sample beta matrix
  data(beta_matrix_miss)

  # generate expected values
  expected_mvals_miss <- log2(beta_matrix_miss /
                           (1 - beta_matrix_miss))

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
