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
