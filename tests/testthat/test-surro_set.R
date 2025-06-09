# Create Expected Values ----

# expected methylation matrix
expected_methyl <- rbind(
  beta_matrix_comp[c("cg02", "cg07", "cg08", "cg13", "cg17"), ],
  matrix(rep(NA, 25), nrow = 5) |>
    `colnames<-`(paste0("samp", 1:5)) |>
    `rownames<-`(paste0("cg", c("03", "06", "11", "15", "18")))
)

# Run Tests ----
test_that("surro_set() function works with weights vector", {

  # load sample data
  data(beta_matrix_comp)
  data(wts_vec_cnt)

  # generate expected values
  expected_vals <- list(
    methyl = expected_methyl,
    weights = wts_vec_cnt[-11],
    intercept = wts_vec_cnt[11]
  )
  class(expected_vals) <- "methyl_surro"

  # run function
  function_vals <- surro_set(
    methyl = beta_matrix_comp,
    weights = wts_vec_cnt,
    intercept = "Intercept"
  )

  # compare result
  expect_identical(function_vals, expected_vals)
})

test_that("surro_set() function works with weights matrix", {

  # load sample data
  data(beta_matrix_comp)
  data(wts_mat)

  # generate expected values
  expected_vals <- list(
    methyl = expected_methyl,
    weights = wts_mat[-11, "wt_cnt"],
    intercept = wts_mat[11, "wt_cnt"] |>
      `names<-`("Intercept")
  )
  class(expected_vals) <- "methyl_surro"

  # run function
  function_vals <- surro_set(
    methyl = beta_matrix_comp,
    weights = wts_mat[, "wt_cnt"],
    intercept = "Intercept"
  )

  # compare result
  expect_identical(function_vals, expected_vals)
})

test_that("surro_set() function works with weights data frame", {

  # load sample data
  data(beta_matrix_comp)
  data(wts_df)

  # generate expected values
  wts_df_vector <- wts_df[, "wt_cnt", drop = FALSE]
  wts_df_vector <- wts_df_vector[["wt_cnt"]] |>
    `names<-`(rownames(wts_df))
  expected_vals <- list(
    methyl = expected_methyl,
    weights = wts_df_vector[-11],
    intercept = wts_df_vector[11] |>
      `names<-`("Intercept")
  )
  class(expected_vals) <- "methyl_surro"

  # run function
  function_vals <- surro_set(
    methyl = beta_matrix_comp,
    weights = wts_df["wt_cnt"],
    intercept = "Intercept"
  )

  # compare result
  expect_identical(function_vals, expected_vals)
})
