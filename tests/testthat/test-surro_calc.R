test_that("surro_calc() function works with complete data and linear transformation", {

  # load sample data
  data(beta_matrix_comp)
  data(wts_vec_lin)
  data(ref_vec_mean)

  # generate methyl_surro object
  methyl_surro_comp_lin <- surro_set(methyl = beta_matrix_comp,
                                     weights = wts_vec_lin,
                                     intercept = "Intercept") |>
    reference_fill(reference = ref_vec_mean,
                   type = "probes")

  # generate expected values
  expected_vals <- as.numeric(
      t(methyl_surro_comp_lin$methy) %*%
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
  data(wts_vec_cnt)
  data(ref_vec_mean)

  # generate methyl_surro object
  methyl_surro_comp_cnt <- surro_set(methyl = beta_matrix_comp,
                                     weights = wts_vec_cnt,
                                     intercept = "Intercept") |>
    reference_fill(reference = ref_vec_mean,
                   type = "probes")

  # generate expected values
  expected_vals <- as.numeric(
    t(methyl_surro_comp_cnt$methy) %*%
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
  data(wts_vec_prb)
  data(ref_vec_mean)

  # generate methyl_surro object
  methyl_surro_comp_prb <- surro_set(methyl = beta_matrix_comp,
                                     weights = wts_vec_prb,
                                     intercept = "Intercept") |>
    reference_fill(reference = ref_vec_mean,
                   type = "probes")

  # generate expected values
  expected_vals <- as.numeric(
    t(methyl_surro_comp_prb$methy) %*%
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
  data(beta_matrix_comp)
  data(wts_vec_lin)
  data(ref_vec_mean)

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
  data(wts_vec_lin)
  data(ref_vec_mean)

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
  data(wts_vec_lin)
  data(ref_vec_mean)

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
