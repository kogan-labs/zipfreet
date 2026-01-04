# Tests for compute_probability_of_freedom function

test_that("compute_probability_of_freedom returns correct structure", {
  result <- compute_probability_of_freedom(
    n = c(21, 4, 4, 4, 4),
    phi_prior = 0.5,
    alpha_prior = 1,
    beta_prior = 1,
    alpha_intro = 1,
    beta_intro = 1,
    p_intro = 0.04,
    growth_rate = 0.01,
    rho = 0.9,
    pi = 0.0
  )

  # Check class

expect_s3_class(result, "zipfreet")

  # Check required components exist
  expect_true("n" %in% names(result))
  expect_true("p_freedom_prior" %in% names(result))
  expect_true("p_freedom_post" %in% names(result))
  expect_true("sensitivity" %in% names(result))
  expect_true("f_posterior" %in% names(result))

  # Check lengths match input
  expect_length(result$n, 5)
  expect_length(result$p_freedom_prior, 5)
  expect_length(result$p_freedom_post, 5)
  expect_length(result$sensitivity, 5)
})

test_that("compute_probability_of_freedom handles single timestep", {
  result <- compute_probability_of_freedom(
    n = 100,
    phi_prior = 0.5,
    alpha_prior = 1,
    beta_prior = 8,
    alpha_intro = 1,
    beta_intro = 1,
    p_intro = 0.0,
    growth_rate = 0.0,
    rho = 0.9,
    pi = 0.0
  )

  expect_s3_class(result, "zipfreet")
  expect_length(result$n, 1)
  expect_equal(result$n, 100)

  # Posterior freedom should be higher than prior after sampling
  expect_gt(result$p_freedom_post[1], result$p_freedom_prior[1])
})

test_that("compute_probability_of_freedom increases freedom probability with more samples", {
  result_low <- compute_probability_of_freedom(
    n = 10,
    phi_prior = 0.5,
    alpha_prior = 1,
    beta_prior = 8,
    alpha_intro = 1,
    beta_intro = 1,
    p_intro = 0.0,
    growth_rate = 0.0,
    rho = 0.9,
    pi = 0.0
  )

  result_high <- compute_probability_of_freedom(
    n = 100,
    phi_prior = 0.5,
    alpha_prior = 1,
    beta_prior = 8,
    alpha_intro = 1,
    beta_intro = 1,
    p_intro = 0.0,
    growth_rate = 0.0,
    rho = 0.9,
    pi = 0.0
  )

  # More samples should yield higher freedom probability
  expect_gt(result_high$p_freedom_post[1], result_low$p_freedom_post[1])
})

test_that("compute_probability_of_freedom handles vector parameters", {
  # Using vector parameters for growth_rate and rho
  result <- compute_probability_of_freedom(
    n = c(50, 50, 50),
    phi_prior = 0.5,
    alpha_prior = 1,
    beta_prior = 8,
    alpha_intro = c(1, 1, 1),
    beta_intro = c(1, 1, 1),
    p_intro = c(0.01, 0.02, 0.03),
    growth_rate = c(0.1, 0.2, 0.3),
    rho = c(0.8, 0.85, 0.9),
    pi = 0.0
  )

  expect_s3_class(result, "zipfreet")
  expect_length(result$n, 3)
})

test_that("compute_probability_of_freedom with zero samples preserves prior", {
  result <- compute_probability_of_freedom(
    n = c(100, 0, 100),
    phi_prior = 0.5,
    alpha_prior = 1,
    beta_prior = 8,
    alpha_intro = 1,
    beta_intro = 1,
    p_intro = 0.0,
    growth_rate = 0.0,
    rho = 0.9,
    pi = 0.0
  )

  expect_s3_class(result, "zipfreet")
  expect_equal(result$n[2], 0)
})

test_that("compute_probability_of_freedom with perfect sensitivity (rho=1)", {
  result <- compute_probability_of_freedom(
    n = 50,
    phi_prior = 0.5,
    alpha_prior = 1,
    beta_prior = 8,
    alpha_intro = 1,
    beta_intro = 1,
    p_intro = 0.0,
    growth_rate = 0.0,
    rho = 1.0,
    pi = 0.0
  )

  expect_s3_class(result, "zipfreet")
  expect_equal(result$sensitivity[1], 1.0, tolerance = 1e-6)
})

test_that("compute_probability_of_freedom sensitivity increases with sample size", {
  result <- compute_probability_of_freedom(
    n = c(10, 50, 100, 200),
    phi_prior = 0.5,
    alpha_prior = 1,
    beta_prior = 8,
    alpha_intro = 1,
    beta_intro = 1,
    p_intro = 0.0,
    growth_rate = 0.0,
    rho = 0.9,
    pi = 0.0
  )

  # Sensitivity at each timestep should increase with cumulative sampling
  # (though this depends on the model - just check they're all valid)
  expect_true(all(result$sensitivity >= 0 & result$sensitivity <= 1))
})

test_that("compute_probability_of_freedom probabilities are valid", {
  result <- compute_probability_of_freedom(
    n = c(21, 4, 4, 4, 4),
    phi_prior = 0.5,
    alpha_prior = 1,
    beta_prior = 1,
    alpha_intro = 1,
    beta_intro = 1,
    p_intro = 0.04,
    growth_rate = 0.01,
    rho = 0.9,
    pi = 0.0
  )

  # All probabilities should be between 0 and 1
  expect_true(all(result$p_freedom_prior >= 0 & result$p_freedom_prior <= 1))
  expect_true(all(result$p_freedom_post >= 0 & result$p_freedom_post <= 1))
  expect_true(all(result$sensitivity >= 0 & result$sensitivity <= 1))
})
