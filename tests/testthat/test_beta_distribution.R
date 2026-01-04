# Tests for beta_distribution_statistics function

test_that("beta_distribution_statistics computes from alpha and beta", {
  result <- beta_distribution_statistics(alpha = 2, beta = 5)

  expect_named(result, c("alpha", "beta", "mu", "sigma"))
  expect_equal(result["alpha"], c(alpha = 2))
  expect_equal(result["beta"], c(beta = 5))

  # Expected mean for Beta(2, 5): alpha / (alpha + beta) = 2/7
  expected_mu <- 2 / (2 + 5)
  expect_equal(unname(result["mu"]), expected_mu, tolerance = 1e-6)

  # Expected variance: alpha*beta / ((alpha+beta)^2 * (alpha+beta+1))
  expected_var <- (2 * 5) / ((2 + 5)^2 * (2 + 5 + 1))
  expected_sigma <- sqrt(expected_var)
  expect_equal(unname(result["sigma"]), expected_sigma, tolerance = 1e-6)
})

test_that("beta_distribution_statistics computes from alpha and mu", {
  alpha <- 3
  mu <- 0.4

  result <- beta_distribution_statistics(alpha = alpha, mu = mu)

  expect_named(result, c("alpha", "beta", "mu", "sigma"))
  expect_equal(unname(result["alpha"]), alpha)

  # Expected beta: alpha * (1-mu) / mu
  expected_beta <- alpha * (1 - mu) / mu
  expect_equal(unname(result["beta"]), expected_beta, tolerance = 1e-6)
})

test_that("beta_distribution_statistics computes from mu and sigma", {
  mu <- 0.3
  sigma <- 0.1

  result <- beta_distribution_statistics(mu = mu, sigma = sigma)

  expect_named(result, c("alpha", "beta", "mu", "sigma"))

  # Verify the computed mu and sigma match inputs
  computed_alpha <- result["alpha"]
  computed_beta <- result["beta"]

  # Recalculate mu from alpha and beta
  recalc_mu <- computed_alpha / (computed_alpha + computed_beta)
  expect_equal(unname(recalc_mu), mu, tolerance = 1e-4)

  # Recalculate sigma from alpha and beta
  recalc_var <- (computed_alpha * computed_beta) /
    ((computed_alpha + computed_beta)^2 * (computed_alpha + computed_beta + 1))
  recalc_sigma <- sqrt(recalc_var)
  expect_equal(unname(recalc_sigma), sigma, tolerance = 1e-4)
})

test_that("beta_distribution_statistics round-trip alpha/beta -> mu/sigma -> alpha/beta", {
  # Start with known alpha and beta
  original_alpha <- 5
  original_beta <- 10

  # Get mu and sigma
  step1 <- beta_distribution_statistics(alpha = original_alpha, beta = original_beta)

  # Use mu and sigma to recover alpha and beta
  step2 <- beta_distribution_statistics(mu = step1["mu"], sigma = step1["sigma"])

  expect_equal(unname(step2["alpha"]), original_alpha, tolerance = 1e-4)
  expect_equal(unname(step2["beta"]), original_beta, tolerance = 1e-4)
})

test_that("beta_distribution_statistics errors on invalid input combinations", {
  # Should error when providing wrong combinations
  expect_error(beta_distribution_statistics(alpha = 1))  # Only alpha
  expect_error(beta_distribution_statistics(beta = 1))   # Only beta
  expect_error(beta_distribution_statistics(mu = 0.5))   # Only mu
  expect_error(beta_distribution_statistics(sigma = 0.1))  # Only sigma
  expect_error(beta_distribution_statistics(alpha = 1, beta = 2, mu = 0.3))  # Three params
  expect_error(beta_distribution_statistics(alpha = 1, beta = 2, mu = 0.3, sigma = 0.1))  # All four
})

test_that("beta_distribution_statistics handles Beta(1,1) - uniform distribution", {
  result <- beta_distribution_statistics(alpha = 1, beta = 1)

  expect_equal(unname(result["mu"]), 0.5)
  # Variance of Uniform(0,1) = 1/12, so sigma = sqrt(1/12)
  expect_equal(unname(result["sigma"]), sqrt(1/12), tolerance = 1e-6)
})

test_that("beta_distribution_statistics handles symmetric distributions", {
  # Beta(a, a) should have mu = 0.5
  for (a in c(1, 2, 5, 10)) {
    result <- beta_distribution_statistics(alpha = a, beta = a)
    expect_equal(unname(result["mu"]), 0.5, tolerance = 1e-6)
  }
})

test_that("beta_distribution_statistics sigma decreases with larger parameters", {
  # Higher alpha and beta (with same ratio) means lower variance
  result_low <- beta_distribution_statistics(alpha = 1, beta = 1)
  result_mid <- beta_distribution_statistics(alpha = 5, beta = 5)
  result_high <- beta_distribution_statistics(alpha = 50, beta = 50)

  expect_gt(result_low["sigma"], result_mid["sigma"])
  expect_gt(result_mid["sigma"], result_high["sigma"])
})

test_that("beta_distribution_statistics handles extreme parameters", {
  # Very small alpha and beta
  result_small <- beta_distribution_statistics(alpha = 0.1, beta = 0.1)
  expect_true(!is.na(result_small["mu"]))
  expect_true(!is.na(result_small["sigma"]))

  # Large alpha and beta
  result_large <- beta_distribution_statistics(alpha = 100, beta = 200)
  expect_equal(unname(result_large["mu"]), 100 / 300, tolerance = 1e-6)
})

test_that("beta_distribution_statistics mu is always between 0 and 1", {
  test_cases <- list(
    c(1, 1),
    c(1, 10),
    c(10, 1),
    c(0.5, 2),
    c(2, 0.5),
    c(100, 100)
  )

  for (params in test_cases) {
    result <- beta_distribution_statistics(alpha = params[1], beta = params[2])
    expect_true(result["mu"] > 0 && result["mu"] < 1)
    expect_true(result["sigma"] > 0)
  }
})
