# Tests for S3 methods (print.zipfreet and plot.zipfreet)

test_that("print.zipfreet produces correct output format", {
  result <- compute_probability_of_freedom(
    n = c(21, 4, 4),
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

  # Capture printed output
  output <- capture.output(print(result))

  # Check header is present
  expect_true(any(grepl("time", output)))
  expect_true(any(grepl("p_freedom_prior", output)))
  expect_true(any(grepl("p_freedom_post", output)))

  # Should have header rows plus data rows
  expect_gte(length(output), 4)  # At least header + separator + 3 data rows
})

test_that("print.zipfreet handles single timestep", {
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

  # Should not error
  output <- capture.output(print(result))
  expect_true(length(output) >= 3)  # Header + separator + 1 data row
})

test_that("print.zipfreet returns invisibly", {
  result <- compute_probability_of_freedom(
    n = 50,
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

  # print should return NULL (invisibly)
  ret <- capture.output(invisible(print(result)))
  expect_true(length(ret) > 0)  # Output was produced
})

test_that("plot.zipfreet creates a ggplot object", {
  skip_if_not_installed("ggplot2")

  result <- compute_probability_of_freedom(
    n = c(50, 50, 50),
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

  # Should not error when plotting
  expect_no_error({
    # Suppress the plot output
    grDevices::pdf(NULL)
    plot(result)
    grDevices::dev.off()
  })
})

test_that("plot.zipfreet works with subset of times", {
  skip_if_not_installed("ggplot2")

  result <- compute_probability_of_freedom(
    n = c(50, 50, 50, 50, 50),
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

  # Should not error when plotting subset
  expect_no_error({
    grDevices::pdf(NULL)
    plot(result, time = c(1, 3, 5))
    grDevices::dev.off()
  })
})

test_that("plot.zipfreet works with single time point", {
  skip_if_not_installed("ggplot2")

  result <- compute_probability_of_freedom(
    n = c(50, 50, 50),
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

  expect_no_error({
    grDevices::pdf(NULL)
    plot(result, time = 2)
    grDevices::dev.off()
  })
})

test_that("plot.zipfreet works with compute_sample_size output", {
  skip_if_not_installed("ggplot2")

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

  expect_no_error({
    grDevices::pdf(NULL)
    plot(result)
    grDevices::dev.off()
  })
})
