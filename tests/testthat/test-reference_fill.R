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

test_that("reference_fill() function correctly handles 'cases' type with reference vector", {

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
                                type = "cases")
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

test_that("reference_fill() function correctly handles 'probes' type with reference matrix", {

  # load sample data
  data(methyl_surro_miss)
  data(ref_mat)

  # generate expected values
  expected_vals <- methyl_surro_miss
  for (probe in c("cg03", "cg06", "cg11", "cg15", "cg18")) {
    for (samp in 1:5) {
      expected_vals$methyl[probe, samp] <- ref_mat[probe, "mean"]
    }
  }

  # run function
  function_vals <- reference_fill(methyl_surro = methyl_surro_miss,
                                  reference = ref_mat,
                                  col_name = "mean",
                                  type = "probes")

  # compare result
  expect_equal(function_vals, expected_vals)
})

test_that("reference_fill() function correctly handles 'cases' type with reference matrix", {

  # load sample data
  data(methyl_surro_miss)
  data(ref_mat)

  # generate expected values
  expected_vals <- methyl_surro_miss
  for (probe in c("cg02", "cg07")) {
    for (samp in 1:5) {
      if (is.na(expected_vals$methyl[probe, samp])) {
        expected_vals$methyl[probe, samp] <- ref_mat[probe, "mean"]
      }
    }
  }

  # run function
  function_vals <- reference_fill(methyl_surro = methyl_surro_miss,
                                  reference = ref_mat,
                                  col_name = "mean",
                                  type = "cases")
  # compare results
  expect_equal(function_vals, expected_vals)
})

test_that("reference_fill() function correctly handles 'all' type with reference matrix", {

  # load sample data
  data(methyl_surro_miss)
  data(ref_mat)

  # generate expected values
  expected_vals <- methyl_surro_miss
  for (probe in rownames(methyl_surro_miss$methyl)) {
    for (samp in 1:5) {
      if (is.na(expected_vals$methyl[probe, samp])) {
        expected_vals$methyl[probe, samp] <- ref_mat[probe, "mean"]
      }
    }
  }

  # run function
  function_vals <- reference_fill(methyl_surro = methyl_surro_miss,
                                  reference = ref_mat,
                                  col_name = "mean",
                                  type = "all")
  # compare results
  expect_equal(function_vals, expected_vals)
})

test_that("reference_fill() function correctly handles 'probes' type with reference data frame", {

  # load sample data
  data(methyl_surro_miss)
  data(ref_df)

  # generate expected values
  expected_vals <- methyl_surro_miss
  for (probe in c("cg03", "cg06", "cg11", "cg15", "cg18")) {
    for (samp in 1:5) {
      expected_vals$methyl[probe, samp] <- ref_df[probe, "mean"]
    }
  }

  # run function
  function_vals <- reference_fill(methyl_surro = methyl_surro_miss,
                                  reference = ref_df,
                                  col_name = "mean",
                                  type = "probes")

  # compare result
  expect_equal(function_vals, expected_vals)
})

test_that("reference_fill() function correctly handles 'cases' type with reference data frame", {

  # load sample data
  data(methyl_surro_miss)
  data(ref_df)

  # generate expected values
  expected_vals <- methyl_surro_miss
  for (probe in c("cg02", "cg07")) {
    for (samp in 1:5) {
      if (is.na(expected_vals$methyl[probe, samp])) {
        expected_vals$methyl[probe, samp] <- ref_df[probe, "mean"]
      }
    }
  }

  # run function
  function_vals <- reference_fill(methyl_surro = methyl_surro_miss,
                                  reference = ref_df,
                                  col_name = "mean",
                                  type = "cases")
  # compare results
  expect_equal(function_vals, expected_vals)
})

test_that("reference_fill() function correctly handles 'all' type with reference data frame", {

  # load sample data
  data(methyl_surro_miss)
  data(ref_df)

  # generate expected values
  expected_vals <- methyl_surro_miss
  for (probe in rownames(methyl_surro_miss$methyl)) {
    for (samp in 1:5) {
      if (is.na(expected_vals$methyl[probe, samp])) {
        expected_vals$methyl[probe, samp] <- ref_df[probe, "mean"]
      }
    }
  }

  # run function
  function_vals <- reference_fill(methyl_surro = methyl_surro_miss,
                                  reference = ref_df,
                                  col_name = "mean",
                                  type = "all")
  # compare results
  expect_equal(function_vals, expected_vals)
})
