# Tests for compute_sample_size function

test_that("compute_sample_size returns correct structure", {
  result <- compute_sample_size(
    phi_prior = 0.5,
    alpha_prior = 1,
    beta_prior = 1,
    alpha_intro = 1,
    beta_intro = 1,
    p_intro = 0.04,
    growth_rate = 0.01,
    rho = 0.9,
    pi = 0,
    dconf = 0.95,
    n_steps = 5
  )

  # Check class
  expect_s3_class(result, "zipfreet")

  # Check required components exist
  expect_true("n" %in% names(result))
  expect_true("p_freedom_prior" %in% names(result))
  expect_true("p_freedom_post" %in% names(result))
  expect_true("sensitivity" %in% names(result))
  expect_true("f_posterior" %in% names(result))

  # Check lengths match n_steps
  expect_length(result$n, 5)
  expect_length(result$p_freedom_prior, 5)
  expect_length(result$p_freedom_post, 5)
})

test_that("compute_sample_size single step reaches desired confidence", {
  dconf <- 0.95

  result <- compute_sample_size(
    phi_prior = 0.5,
    alpha_prior = 1,
    beta_prior = 8,
    alpha_intro = 1,
    beta_intro = 1,
    p_intro = 0.0,
    growth_rate = 0.0,
    rho = 0.9,
    pi = 0,
    dconf = dconf,
    n_steps = 1,
    method = "restore"
  )

  expect_s3_class(result, "zipfreet")
  # Posterior should meet or exceed desired confidence
  expect_gte(result$p_freedom_post[1], dconf - 0.01)  # Small tolerance for numerical precision
})

test_that("compute_sample_size works with 'restore' method", {
  result <- compute_sample_size(
    phi_prior = 0.5,
    alpha_prior = 1,
    beta_prior = 8,
    alpha_intro = 1,
    beta_intro = 1,
    p_intro = 0.03,
    growth_rate = 0.1,
    rho = 0.9,
    pi = 0,
    dconf = 0.95,
    n_steps = 3,
    method = "restore"
  )

  expect_s3_class(result, "zipfreet")
  expect_length(result$n, 3)
  # Each step should achieve the desired confidence
  expect_true(all(result$p_freedom_post >= 0.94))  # Allow small tolerance
})

test_that("compute_sample_size works with 'maintain' method", {
  result <- compute_sample_size(
    phi_prior = 0.5,
    alpha_prior = 1,
    beta_prior = 8,
    alpha_intro = 1,
    beta_intro = 1,
    p_intro = 0.03,
    growth_rate = 0.1,
    rho = 0.9,
    pi = 0,
    dconf = 0.95,
    n_steps = 3,
    method = "maintain"
  )

  expect_s3_class(result, "zipfreet")
  expect_length(result$n, 3)
})

test_that("compute_sample_size returns non-negative sample sizes", {
  result <- compute_sample_size(
    phi_prior = 0.5,
    alpha_prior = 1,
    beta_prior = 8,
    alpha_intro = 1,
    beta_intro = 1,
    p_intro = 0.03,
    growth_rate = 0.1,
    rho = 0.9,
    pi = 0,
    dconf = 0.95,
    n_steps = 5
  )

  expect_true(all(result$n >= 0))
})

test_that("compute_sample_size with higher confidence requires more samples", {
  result_low <- compute_sample_size(
    phi_prior = 0.5,
    alpha_prior = 1,
    beta_prior = 8,
    alpha_intro = 1,
    beta_intro = 1,
    p_intro = 0.0,
    growth_rate = 0.0,
    rho = 0.9,
    pi = 0,
    dconf = 0.90,
    n_steps = 1
  )

  result_high <- compute_sample_size(
    phi_prior = 0.5,
    alpha_prior = 1,
    beta_prior = 8,
    alpha_intro = 1,
    beta_intro = 1,
    p_intro = 0.0,
    growth_rate = 0.0,
    rho = 0.9,
    pi = 0,
    dconf = 0.99,
    n_steps = 1
  )

  # Higher confidence should require more samples
  expect_gt(result_high$n[1], result_low$n[1])
})

test_that("compute_sample_size with lower sensitivity requires more samples", {
  result_high_sens <- compute_sample_size(
    phi_prior = 0.5,
    alpha_prior = 1,
    beta_prior = 8,
    alpha_intro = 1,
    beta_intro = 1,
    p_intro = 0.0,
    growth_rate = 0.0,
    rho = 0.95,
    pi = 0,
    dconf = 0.95,
    n_steps = 1
  )

  result_low_sens <- compute_sample_size(
    phi_prior = 0.5,
    alpha_prior = 1,
    beta_prior = 8,
    alpha_intro = 1,
    beta_intro = 1,
    p_intro = 0.0,
    growth_rate = 0.0,
    rho = 0.7,
    pi = 0,
    dconf = 0.95,
    n_steps = 1
  )

  # Lower test sensitivity should require more samples
  expect_gt(result_low_sens$n[1], result_high_sens$n[1])
})

test_that("compute_sample_size handles vector dconf", {
  result <- compute_sample_size(
    phi_prior = 0.5,
    alpha_prior = 1,
    beta_prior = 8,
    alpha_intro = 1,
    beta_intro = 1,
    p_intro = 0.03,
    growth_rate = 0.1,
    rho = 0.9,
    pi = 0,
    dconf = c(0.90, 0.95, 0.99),  # Different confidence levels at each step
    n_steps = 3
  )

  expect_s3_class(result, "zipfreet")
  expect_length(result$n, 3)
})

test_that("compute_sample_size with high prior requires fewer samples", {
  result_low_prior <- compute_sample_size(
    phi_prior = 0.3,  # Low prior freedom probability
    alpha_prior = 1,
    beta_prior = 8,
    alpha_intro = 1,
    beta_intro = 1,
    p_intro = 0.0,
    growth_rate = 0.0,
    rho = 0.9,
    pi = 0,
    dconf = 0.95,
    n_steps = 1
  )

  result_high_prior <- compute_sample_size(
    phi_prior = 0.8,  # High prior freedom probability
    alpha_prior = 1,
    beta_prior = 8,
    alpha_intro = 1,
    beta_intro = 1,
    p_intro = 0.0,
    growth_rate = 0.0,
    rho = 0.9,
    pi = 0,
    dconf = 0.95,
    n_steps = 1
  )

  # Higher prior freedom should require fewer samples
  expect_lt(result_high_prior$n[1], result_low_prior$n[1])
})

test_that("compute_sample_size respects n_max constraint", {
  result <- compute_sample_size(
    phi_prior = 0.1,  # Very low prior
    alpha_prior = 1,
    beta_prior = 8,
    alpha_intro = 1,
    beta_intro = 1,
    p_intro = 0.0,
    growth_rate = 0.0,
    rho = 0.5,  # Low sensitivity
    pi = 0,
    dconf = 0.99,  # Very high confidence
    n_steps = 1,
    n_max = 500
  )

  # Sample size should not exceed n_max
  expect_lte(result$n[1], 500)
})

test_that("compute_sample_size probabilities are valid", {
  result <- compute_sample_size(
    phi_prior = 0.5,
    alpha_prior = 1,
    beta_prior = 8,
    alpha_intro = 1,
    beta_intro = 1,
    p_intro = 0.03,
    growth_rate = 0.1,
    rho = 0.9,
    pi = 0,
    dconf = 0.95,
    n_steps = 3
  )

  # All probabilities should be between 0 and 1
  expect_true(all(result$p_freedom_prior >= 0 & result$p_freedom_prior <= 1))
  expect_true(all(result$p_freedom_post >= 0 & result$p_freedom_post <= 1))
  expect_true(all(result$sensitivity >= 0 & result$sensitivity <= 1))
})

# Tests for sampling_schedule parameter

test_that("sampling_schedule forces n=0 at specified timesteps", {
  result <- compute_sample_size(
    phi_prior = 0.5,
    alpha_prior = 1,
    beta_prior = 8,
    alpha_intro = 1,
    beta_intro = 1,
    p_intro = 0.03,
    growth_rate = 0.1,
    rho = 0.9,
    pi = 0,
    dconf = 0.95,
    n_steps = 5,
    sampling_schedule = c(TRUE, FALSE, TRUE, FALSE, TRUE)
  )

  expect_s3_class(result, "zipfreet")
  expect_length(result$n, 5)

  # Timesteps 2 and 4 should have n = 0
  expect_equal(result$n[2], 0)
  expect_equal(result$n[4], 0)

  # Timesteps 1, 3, and 5 should have n > 0 (computed sample sizes)
  expect_gt(result$n[1], 0)
  expect_gt(result$n[3], 0)
  expect_gt(result$n[5], 0)
})

test_that("sampling_schedule with all TRUE behaves like default", {
  result_default <- compute_sample_size(
    phi_prior = 0.5,
    alpha_prior = 1,
    beta_prior = 8,
    alpha_intro = 1,
    beta_intro = 1,
    p_intro = 0.03,
    growth_rate = 0.1,
    rho = 0.9,
    pi = 0,
    dconf = 0.95,
    n_steps = 3
  )

  result_explicit <- compute_sample_size(
    phi_prior = 0.5,
    alpha_prior = 1,
    beta_prior = 8,
    alpha_intro = 1,
    beta_intro = 1,
    p_intro = 0.03,
    growth_rate = 0.1,
    rho = 0.9,
    pi = 0,
    dconf = 0.95,
    n_steps = 3,
    sampling_schedule = c(TRUE, TRUE, TRUE)
  )

  # Results should be identical
  expect_equal(result_default$n, result_explicit$n)
  expect_equal(result_default$p_freedom_post, result_explicit$p_freedom_post)
})

test_that("sampling_schedule accepts numeric 0/1 values", {
  result <- compute_sample_size(
    phi_prior = 0.5,
    alpha_prior = 1,
    beta_prior = 8,
    alpha_intro = 1,
    beta_intro = 1,
    p_intro = 0.03,
    growth_rate = 0.1,
    rho = 0.9,
    pi = 0,
    dconf = 0.95,
    n_steps = 3,
    sampling_schedule = c(1, 0, 1)  # Numeric instead of logical
  )

  expect_s3_class(result, "zipfreet")
  expect_equal(result$n[2], 0)
  expect_gt(result$n[1], 0)
  expect_gt(result$n[3], 0)
})

test_that("sampling_schedule errors on length mismatch", {
  expect_error(
    compute_sample_size(
      phi_prior = 0.5,
      alpha_prior = 1,
      beta_prior = 8,
      alpha_intro = 1,
      beta_intro = 1,
      p_intro = 0.03,
      growth_rate = 0.1,
      rho = 0.9,
      pi = 0,
      dconf = 0.95,
      n_steps = 5,
      sampling_schedule = c(TRUE, FALSE, TRUE)  # Length 3, but n_steps = 5
    ),
    "sampling_schedule must have length"
  )
})

test_that("sampling_schedule errors on invalid values", {
  expect_error(
    compute_sample_size(
      phi_prior = 0.5,
      alpha_prior = 1,
      beta_prior = 8,
      alpha_intro = 1,
      beta_intro = 1,
      p_intro = 0.03,
      growth_rate = 0.1,
      rho = 0.9,
      pi = 0,
      dconf = 0.95,
      n_steps = 3,
      sampling_schedule = c(1, 2, 1)  # Invalid value 2
    ),
    "sampling_schedule must be a logical vector"
  )
})

test_that("sampling_schedule with skip at start still computes correctly", {
  result <- compute_sample_size(
    phi_prior = 0.5,
    alpha_prior = 1,
    beta_prior = 8,
    alpha_intro = 1,
    beta_intro = 1,
    p_intro = 0.03,
    growth_rate = 0.1,
    rho = 0.9,
    pi = 0,
    dconf = 0.95,
    n_steps = 3,
    sampling_schedule = c(FALSE, TRUE, TRUE)
  )

  expect_equal(result$n[1], 0)
  expect_gt(result$n[2], 0)
  expect_gt(result$n[3], 0)
})

test_that("sampling_schedule with all FALSE returns all zeros", {
  result <- compute_sample_size(
    phi_prior = 0.5,
    alpha_prior = 1,
    beta_prior = 8,
    alpha_intro = 1,
    beta_intro = 1,
    p_intro = 0.03,
    growth_rate = 0.1,
    rho = 0.9,
    pi = 0,
    dconf = 0.95,
    n_steps = 3,
    sampling_schedule = c(FALSE, FALSE, FALSE)
  )

  expect_equal(result$n, c(0, 0, 0))
})

test_that("sampling_schedule works with maintain method", {
  result <- compute_sample_size(
    phi_prior = 0.5,
    alpha_prior = 1,
    beta_prior = 8,
    alpha_intro = 1,
    beta_intro = 1,
    p_intro = 0.03,
    growth_rate = 0.1,
    rho = 0.9,
    pi = 0,
    dconf = 0.95,
    n_steps = 5,
    method = "maintain",
    sampling_schedule = c(TRUE, FALSE, TRUE, FALSE, TRUE)
  )

  expect_s3_class(result, "zipfreet")
  expect_equal(result$n[2], 0)
  expect_equal(result$n[4], 0)
})
