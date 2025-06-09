# Create Expected Values ----

# complete beta matrix
expected_beta_comp <- c(
  0.1028, 0.3203, 0.4829, 0.7205, 0.3694,
  0.2875, 0.7883, 0.4089, 0.8830, 0.9404,
  0.4348, 0.1876, 0.8903, 0.1422, 0.9842,
  0.9849, 0.7822, 0.9144, 0.5492, 0.1542,
  0.8998, 0.2460, 0.0420, 0.3279, 0.9545,
  0.8895, 0.6928, 0.6405, 0.9942, 0.6557,
  0.8930, 0.0935, 0.6087, 0.9540, 0.0910,
  0.8864, 0.4667, 0.4106, 0.5854, 0.1419,
  0.1750, 0.5115, 0.1470, 0.4045, 0.6900,
  0.9630, 0.9022, 0.6907, 0.7954, 0.0246,
  0.1306, 0.5999, 0.9352, 0.6478, 0.6192,
  0.6531, 0.3328, 0.3012, 0.3198, 0.8913,
  0.1428, 0.4145, 0.4137, 0.3688, 0.1524,
  0.3435, 0.4886, 0.0607, 0.3077, 0.6729,
  0.6567, 0.9544, 0.9477, 0.2197, 0.7370
) |>
  matrix(nrow = 15, byrow = TRUE) |>
  `colnames<-`(paste0("samp", 1:5)) |>
  `rownames<-`(paste0("cg", c("01", "02", "04",
                              "05", "07", "08",
                              "09", "10", "12",
                              "13", "14", "16",
                              "17", "19", "20")))

# missing beta matrix
miss <- data.frame(
  row = c(1, 2, 2, 2, 5,
          7, 7, 8, 8, 9,
          9, 11, 11, 14, 15),
  col = c(4, 2, 3, 5, 3,
          1, 5, 2, 3, 2,
          4, 2, 5, 5, 5)
)
expected_beta_miss <- expected_beta_comp
for (i in seq_len(nrow(miss))) {
  expected_beta_miss[miss$row[i], miss$col[i]] <- NA
}
rm("miss")

# complete M-value matrix
expected_mval_comp <- log2(beta_matrix_comp /
                           (1 - beta_matrix_comp))

# missing M-value matrix
expected_mval_miss <- log2(beta_matrix_miss /
                             (1 - beta_matrix_miss))

# weights
wts_names <- c("cg02", "cg03", "cg06", "cg07",
               "cg08", "cg11", "cg13", "cg15",
               "cg17", "cg18", "Intercept")
wts_vec_lin_expected <- c(
  -0.009083377, -0.001155999, 0.005978497,
  -0.007562015, 0.001218960, -0.005869372,
  -0.007449367, 0.005066157, 0.007900907,
  -0.002510744, 1.211000000
) |>
  `names<-`(wts_names)
wts_vec_prb_expected <- c(
  0.16511519, -0.40515934, -0.11603036,
  -0.22561636, 0.31464004, -0.05148366,
  0.31006435, 0.31238951, 0.29434232,
  -0.06016831, 0.01900000
) |>
  `names<-`(wts_names)
wts_vec_cnt_expected <- c(
  0.050895032, 0.025844226, 0.042036480,
  -0.099875045, -0.004936685, -0.055976223,
  -0.024036692, 0.022554201, -0.029640418,
  -0.077772915, 0.937000000
) |>
  `names<-`(wts_names)
wts_mat_expected <- c(
  wts_vec_lin_expected,
  wts_vec_prb_expected,
  wts_vec_cnt_expected
) |>
  matrix(ncol = 3, byrow = FALSE) |>
  `row.names<-`(wts_names) |>
  `colnames<-`(c("wt_lin", "wt_prb", "wt_cnt"))
wts_df_expected <- data.frame(wts_mat_expected)
rm("wts_names")

# reference data frame
reference_expected = data.frame(
  mean = c(0.3992, 0.6616, 0.4948, 0.5278,
           0.6770, 0.5526, 0.4940, 0.7745,
           0.5281, 0.4982, 0.4566, 0.3856,
           0.6752, 0.5866, 0.4004, 0.4996,
           0.2984, 0.3923, 0.3747, 0.7031),
  median = c(0.3694, 0.7883, 0.5281, 0.4348,
             0.7822, 0.5726, 0.3279, 0.6928,
             0.6087, 0.4667, 0.5440, 0.4045,
             0.7954, 0.6192, 0.3181, 0.3328,
             0.3688, 0.2659, 0.3435, 0.7370)
) |>
  `rownames<-`(paste0("cg", c("01", "02", "03", "04",
                              "05", "06", "07", "08",
                              "09", "10", "11", "12",
                              "13", "14", "15", "16",
                              "17", "18", "19", "20")))

# Run Tests ----

test_that("beta_matrix_comp is correct", {

  # load sample data
  data(beta_matrix_comp)

  # compare result
  expect_equal(beta_matrix_comp, expected_beta_comp, tolerance = 0.0001)
})

test_that("beta_matrix_miss is correct", {

  # load sample data
  data(beta_matrix_miss)

  # compare result
  expect_equal(beta_matrix_miss, expected_beta_miss, tolerance = 0.0001)
})

test_that("mval_matrix_comp is correct", {

  # load sample data
  data(mval_matrix_comp)

  # compare result
  expect_equal(mval_matrix_comp, expected_mval_comp, tolerance = 0.0001)
})

test_that("mval_matrix_comp is correct", {

  # load sample data
  data(mval_matrix_miss)

  # compare result
  expect_equal(mval_matrix_miss, expected_mval_miss, tolerance = 0.0001)
})

test_that("refernce data frame is correct", {

  # load sample data
  data(reference)

  # compare result
  expect_equal(reference, reference_expected, tolerance = 0.0001)
})

test_that("wts_vec_lin is correct", {

  # load sample data
  data(wts_vec_lin)

  # compare result
  expect_equal(wts_vec_lin, wts_vec_lin_expected, tolerance = 0.0001)
})

test_that("wts_vec_prb is correct", {

  # load sample data
  data(wts_vec_prb)

  # compare result
  expect_equal(wts_vec_prb, wts_vec_prb_expected, tolerance = 0.0001)
})

test_that("wts_vec_cnt is correct", {

  # load sample data
  data(wts_vec_cnt)

  # compare result
  expect_equal(wts_vec_cnt, wts_vec_cnt_expected, tolerance = 0.0001)
})

test_that("wts_mat is correct", {

  # load sample data
  data(wts_mat)

  # compare result
  expect_equal(wts_mat, wts_mat_expected, tolerance = 0.0001)
})

test_that("wts_df is correct", {

  # load sample data
  data(wts_df)

  # compare result
  expect_equal(wts_df, wts_df_expected, tolerance = 0.0001)
})




