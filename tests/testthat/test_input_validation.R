# Tests for input validation

# compute_probability_of_freedom validation tests
test_that("compute_probability_of_freedom errors on vector phi_prior", {
  expect_error(
    compute_probability_of_freedom(
      n = c(50, 50),
      phi_prior = c(0.5, 0.6),  # Should be length 1
      alpha_prior = 1,
      beta_prior = 8,
      alpha_intro = 1,
      beta_intro = 1,
      p_intro = 0.0,
      growth_rate = 0.0,
      rho = 0.9,
      pi = 0.0
    ),
    "phi_prior must be length 1"
  )
})

test_that("compute_probability_of_freedom errors on vector alpha_prior", {
  expect_error(
    compute_probability_of_freedom(
      n = c(50, 50),
      phi_prior = 0.5,
      alpha_prior = c(1, 2),  # Should be length 1
      beta_prior = 8,
      alpha_intro = 1,
      beta_intro = 1,
      p_intro = 0.0,
      growth_rate = 0.0,
      rho = 0.9,
      pi = 0.0
    ),
    "alpha_prior must be length 1"
  )
})

test_that("compute_probability_of_freedom errors on vector beta_prior", {
  expect_error(
    compute_probability_of_freedom(
      n = c(50, 50),
      phi_prior = 0.5,
      alpha_prior = 1,
      beta_prior = c(8, 9),  # Should be length 1
      alpha_intro = 1,
      beta_intro = 1,
      p_intro = 0.0,
      growth_rate = 0.0,
      rho = 0.9,
      pi = 0.0
    ),
    "beta_prior must be length 1"
  )
})

test_that("compute_probability_of_freedom errors on mismatched vector lengths", {
  expect_error(
    compute_probability_of_freedom(
      n = c(50, 50, 50),
      phi_prior = 0.5,
      alpha_prior = 1,
      beta_prior = 8,
      alpha_intro = c(1, 2),  # Length 2, but n has length 3
      beta_intro = 1,
      p_intro = 0.0,
      growth_rate = 0.0,
      rho = 0.9,
      pi = 0.0
    ),
    "Vector length mismatch"
  )
})

test_that("compute_probability_of_freedom errors on mismatched rho length", {
  expect_error(
    compute_probability_of_freedom(
      n = c(50, 50, 50),
      phi_prior = 0.5,
      alpha_prior = 1,
      beta_prior = 8,
      alpha_intro = 1,
      beta_intro = 1,
      p_intro = 0.0,
      growth_rate = 0.0,
      rho = c(0.9, 0.8),  # Length 2, but n has length 3
      pi = 0.0
    ),
    "Vector length mismatch"
  )
})

# compute_sample_size validation tests
test_that("compute_sample_size errors on vector phi_prior", {
  expect_error(
    compute_sample_size(
      phi_prior = c(0.5, 0.6),  # Should be length 1
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
    ),
    "phi_prior must be length 1"
  )
})

test_that("compute_sample_size errors on vector alpha_prior", {
  expect_error(
    compute_sample_size(
      phi_prior = 0.5,
      alpha_prior = c(1, 2),  # Should be length 1
      beta_prior = 8,
      alpha_intro = 1,
      beta_intro = 1,
      p_intro = 0.0,
      growth_rate = 0.0,
      rho = 0.9,
      pi = 0,
      dconf = 0.95,
      n_steps = 1
    ),
    "alpha_prior must be length 1"
  )
})

test_that("compute_sample_size errors on vector beta_prior", {
  expect_error(
    compute_sample_size(
      phi_prior = 0.5,
      alpha_prior = 1,
      beta_prior = c(8, 9),  # Should be length 1
      alpha_intro = 1,
      beta_intro = 1,
      p_intro = 0.0,
      growth_rate = 0.0,
      rho = 0.9,
      pi = 0,
      dconf = 0.95,
      n_steps = 1
    ),
    "beta_prior must be length 1"
  )
})

test_that("compute_sample_size errors on invalid method", {
  expect_error(
    compute_sample_size(
      phi_prior = 0.5,
      alpha_prior = 1,
      beta_prior = 8,
      alpha_intro = 1,
      beta_intro = 1,
      p_intro = 0.0,
      growth_rate = 0.0,
      rho = 0.9,
      pi = 0,
      dconf = 0.95,
      n_steps = 1,
      method = "invalid_method"  # Should be "restore" or "maintain"
    ),
    "Invalid method"
  )
})

test_that("compute_sample_size errors on mismatched vector lengths", {
  expect_error(
    compute_sample_size(
      phi_prior = 0.5,
      alpha_prior = 1,
      beta_prior = 8,
      alpha_intro = c(1, 2),  # Length 2
      beta_intro = c(1, 2, 3),  # Length 3 - mismatch
      p_intro = 0.0,
      growth_rate = 0.0,
      rho = 0.9,
      pi = 0,
      dconf = 0.95
    ),
    "Vector length mismatch"
  )
})

# is_valid_length tests (internal function, accessed via package internals)
test_that("is_valid_length correctly validates vector lengths", {
  # This tests the internal function behavior via the exported functions
  # Length 1 should always be valid
  expect_no_error(
    compute_probability_of_freedom(
      n = c(50, 50, 50),
      phi_prior = 0.5,
      alpha_prior = 1,
      beta_prior = 8,
      alpha_intro = 1,  # Length 1 - valid
      beta_intro = 1,
      p_intro = 0.0,
      growth_rate = 0.0,
      rho = 0.9,
      pi = 0.0
    )
  )

  # Length matching n should be valid
  expect_no_error(
    compute_probability_of_freedom(
      n = c(50, 50, 50),
      phi_prior = 0.5,
      alpha_prior = 1,
      beta_prior = 8,
      alpha_intro = c(1, 1, 1),  # Length 3 - matches n
      beta_intro = c(1, 1, 1),
      p_intro = c(0.0, 0.0, 0.0),
      growth_rate = c(0.0, 0.0, 0.0),
      rho = c(0.9, 0.9, 0.9),
      pi = c(0.0, 0.0, 0.0)
    )
  )
})

# Edge case tests
test_that("compute_probability_of_freedom handles extreme phi_prior values", {
  # Very low prior
  result_low <- compute_probability_of_freedom(
    n = 100,
    phi_prior = 0.01,
    alpha_prior = 1,
    beta_prior = 8,
    alpha_intro = 1,
    beta_intro = 1,
    p_intro = 0.0,
    growth_rate = 0.0,
    rho = 0.9,
    pi = 0.0
  )
  expect_s3_class(result_low, "zipfreet")

  # Very high prior
  result_high <- compute_probability_of_freedom(
    n = 100,
    phi_prior = 0.99,
    alpha_prior = 1,
    beta_prior = 8,
    alpha_intro = 1,
    beta_intro = 1,
    p_intro = 0.0,
    growth_rate = 0.0,
    rho = 0.9,
    pi = 0.0
  )
  expect_s3_class(result_high, "zipfreet")
})

test_that("compute_sample_size handles extreme desired confidence", {
  # Very low confidence
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
    dconf = 0.60,  # Low confidence
    n_steps = 1
  )
  expect_s3_class(result_low, "zipfreet")

  # Very high confidence
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
    dconf = 0.999,  # Very high confidence
    n_steps = 1
  )
  expect_s3_class(result_high, "zipfreet")
})
