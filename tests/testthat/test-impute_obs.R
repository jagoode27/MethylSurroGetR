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
