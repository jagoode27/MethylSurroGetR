test_that("methyl_miss() function works", {

  # load sample data
  data(beta_matrix_miss)
  data(wts_vec_lin)

  # generate expected values
  expected_vals <- list(
    missing_cases = c(0.6, 0.2) |>
      `names<-`(c("cg02", "cg07")),
    missing_probes = c("cg03", "cg06", "cg11", "cg15", "cg18")
  )

  # run function
  function_vals <- surro_set(
    methyl = beta_matrix_miss,
    weights = wts_vec_cnt,
    intercept = "Intercept"
  ) |>
    methyl_miss()

  # compare result
  expect_identical(function_vals, expected_vals)
})
