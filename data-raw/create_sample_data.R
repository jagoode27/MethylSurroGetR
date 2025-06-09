# Step 1: Set Up ----

# Clear environment
rm(list = ls())
gc()

# Load packages
library(tidyverse)

# Set the seed for reproducibility
set.seed(123)

# Set simulation parameters
sim <- list(n_cases = 5,
            n_wts = 10,
            n_ref = 20,
            n_beta = 15,
            prb_weights = c("cg02", "cg03", "cg06", "cg07", "cg08",
                            "cg11", "cg13", "cg15", "cg17", "cg18"),
            prb_beta = c("cg01", "cg02", "cg04", "cg05", "cg07",
                         "cg08", "cg09", "cg10", "cg12", "cg13",
                         "cg14", "cg16", "cg17", "cg19", "cg20"),
            miss_prop = 0.2
)

# Step 2: Simulate a methylation beta matrix ----

# Generate the methylation beta matrix
beta_matrix <- matrix(runif(sim$n_cases * sim$n_wts, 0, 1),
                      nrow = sim$n_cases, ncol = sim$n_wts) |>
  `colnames<-`(paste0("probe", 1:sim$n_wts)) |>
  `rownames<-`(paste0("case", 1:sim$n_cases))

# Step 3: Simulate linear regression outcome and store weights -----

# Coefficients for linear regression
coeffs_linear <- runif(sim$n_wts, -0.01, 0.01)

# Linear regression outcome
outcome_linear <- beta_matrix %*% coeffs_linear

# Scale the outcome to be between -2 and 2
outcome_linear <- scales::rescale(outcome_linear, to = c(-2, 2))

# Step 4: Simulate logistic regression outcome and store weights -----

# Coefficients for logistic regression
coeffs_logistic <- runif(sim$n_wts, -0.5, 0.5)

# Linear predictor for logistic regression
logit_scores <- beta_matrix %*% coeffs_logistic

# Convert to probabilities using the logistic function
probabilities_logistic <- 1 / (1 + exp(-logit_scores))

# Scale probabilities between 0.1 and 0.9
probabilities_logistic <- scales::rescale(probabilities_logistic, to = c(0.1, 0.9))

# Step 5: Simulate poisson regression outcome and store weights -----

# Coefficients for Poisson regression
coeffs_poisson <- runif(sim$n_wts, -0.1, 0.1)

# Linear predictor for Poisson regression
log_lambda <- beta_matrix %*% coeffs_poisson

# Convert to counts using the exponential function and rounding to the nearest integer
lambda <- exp(log_lambda)

# Scale lambda between 0 and 5
lambda <- scales::rescale(lambda, to = c(0, 5))

# Generate Poisson-distributed counts
outcome_poisson <- rpois(sim$n_cases, lambda)

# Step 6: Compile results ----

# Combine all outcomes into a data frame
regression_weights <- data.frame(
  wt_lin = coeffs_linear,
  wt_prb = coeffs_logistic,
  wt_cnt = coeffs_poisson
) |>
  `row.names<-`(c("cg02", "cg03", "cg06", "cg07", "cg08",
                  "cg11", "cg13", "cg15", "cg17", "cg18"))

# Clean up
rm(list = setdiff(ls(), c("regression_weights", "beta_matrix", "sim")))

# Step 7: Build final objects ----

# Storing regression coefficients into a list
wts_df <- rbind(
  regression_weights |>
    `rownames<-`(sim$prb_weights),
  data.frame(wt_lin = 1.211,
             wt_prb = 0.019,
             wt_cnt = 0.937) |>
  `rownames<-`("Intercept")
)

# Reformat beta matrix
beta_matrix <- beta_matrix |>
  `colnames<-`(sim$prb_weights) |>
  `rownames<-`(paste0("samp", 1:sim$n_cases)) |>
  t()

# Define full set of probes to subset into final objects
full_matrix <- runif(((sim$n_ref - sim$n_wts) * sim$n_cases), 0, 1) |>
  matrix(ncol = sim$n_cases) |>
  `rownames<-`(setdiff(paste0("cg", sprintf(paste0("%0", 2, "d"), 1:sim$n_ref)),
                       rownames(beta_matrix))) |>
  `colnames<-`(paste0("samp", 1:sim$n_cases)) |>
  rbind(beta_matrix) |>
  data.frame() |>
  rownames_to_column(var = "probe") |>
  arrange(probe)

# Create the complete beta matrix for sample
beta_matrix_comp <- full_matrix |>
  filter(probe %in% sim$prb_beta) |>
  column_to_rownames(var = "probe") |>
  as.matrix()

# Create the beta matrix for sample w/ missing values
miss_index <- sample(1:length(beta_matrix_comp),
                     round(length(beta_matrix_comp) * sim$miss_prop, 0))
beta_matrix_miss <- beta_matrix_comp
beta_matrix_miss[miss_index] <-NA
rm(miss_index)

# Create the m-value versions matrix for sample
mval_matrix_comp <- log2(beta_matrix_comp / (1 - beta_matrix_comp))
mval_matrix_miss <- log2(beta_matrix_miss / (1 - beta_matrix_miss))

# Generate data frame version of reference values
ref_df <- full_matrix |>
  column_to_rownames(var = "probe") |>
  apply(1, function(x) {
    cbind(mean = mean(x), median = median(x)) |>
      data.frame()
  }) |>
  bind_rows() |>
  `rownames<-`(paste0("cg", sprintf(paste0("%0", 2, "d"), 1:sim$n_ref)))

# Generate matrix version of reference values
ref_mat <- as.matrix(ref_df)

# Generate vector versions of reference values
ref_vec_mean <- ref_df |>
  pull(mean) |>
  `names<-`(rownames(ref_df))
ref_vec_median <- ref_df |>
  pull(median) |>
  `names<-`(rownames(ref_df))

# Generate matrix version of weights
wts_mat <- as.matrix(wts_df)

# Generate vector versions of weights
wts_vec_lin <- wts_df |>
  pull(wt_lin) |>
  `names<-`(rownames(wts_df))
wts_vec_prb <- wts_df |>
  pull(wt_prb) |>
  `names<-`(rownames(wts_df))
wts_vec_cnt <- wts_df |>
  pull(wt_cnt) |>
  `names<-`(rownames(wts_df))

# Create a complete methyl_surro object
methyl_surro_comp <- surro_set(methyl = beta_matrix_comp,
                               weights = wts_vec_lin,
                               intercept = "Intercept")

# Create a methyl_surro object with missing values
methyl_surro_miss <- surro_set(methyl = beta_matrix_miss,
                               weights = wts_vec_lin,
                               intercept = "Intercept")

# Clean up
rm(list = setdiff(ls(), c("beta_matrix_comp", "beta_matrix_miss",
                          "mval_matrix_comp", "mval_matrix_miss",
                          "wts_df", "wts_mat", "wts_vec_lin",
                          "wts_vec_prb", "wts_vec_cnt",
                          "ref_df", "ref_mat",
                          "ref_vec_mean", "ref_vec_median",
                          "methyl_surro_miss", "methyl_surro_comp")))

# Export to package
usethis::use_data(beta_matrix_comp,
                  beta_matrix_miss,
                  mval_matrix_comp,
                  mval_matrix_miss,
                  wts_df,
                  wts_mat,
                  wts_vec_lin,
                  wts_vec_prb,
                  wts_vec_cnt,
                  ref_df,
                  ref_mat,
                  ref_vec_mean,
                  ref_vec_median,
                  methyl_surro_comp,
                  methyl_surro_miss,
                  overwrite = TRUE)
