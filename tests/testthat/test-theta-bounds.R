# Parameter values are validated against params_bounds (open intervals) by every
# generic, so out-of-domain input errors clearly instead of yielding NaN.

test_that("out-of-domain parameters are rejected with an informative message", {
  d <- gaussian_distrib()

  expect_error(distrib_pdf(d, 0, list(mu = 0, sigma = -1)), "outside its domain")
  expect_error(distrib_pdf(d, 0, list(mu = 0, sigma = -1)), "sigma")
  expect_error(distrib_pdf(d, 0, list(mu = 0, sigma = 0)), "\\(0, Inf\\)")
  expect_error(distrib_pdf(d, 0, list(mu = NA, sigma = 1)), "mu")
  expect_error(distrib_pdf(d, 0, list(mu = Inf, sigma = 1)), "mu")

  # bounded parameters: open interval on both sides
  b <- bernoulli_distrib()
  expect_error(distrib_pdf(b, 1, list(mu = 0)), "outside its domain")
  expect_error(distrib_pdf(b, 1, list(mu = 1)), "outside its domain")
  expect_error(distrib_pdf(b, 1, list(mu = 1.2)), "outside its domain")
  expect_silent(distrib_pdf(b, 1, list(mu = 0.5)))
})

test_that("the check applies to every generic, not just the density", {
  d <- gaussian_distrib()
  bad <- list(mu = 0, sigma = -2)

  expect_error(distrib_cdf(d, 0, bad), "outside its domain")
  expect_error(distrib_quantile(d, 0.5, bad), "outside its domain")
  expect_error(distrib_rng(d, 1, bad), "outside its domain")
  expect_error(distrib_gradient(d, 0, bad), "outside its domain")
  expect_error(distrib_hessian(d, 0, bad), "outside its domain")
  expect_error(distrib_expected_hessian(d, 0, bad), "outside its domain")
  expect_error(distrib_deriv3(d, 0, bad), "outside its domain")
  expect_error(distrib_deriv4(d, 0, bad), "outside its domain")
  expect_error(distrib_grad_y(d, 0, bad), "outside its domain")
})

test_that("check_theta_bounds is usable directly and reports every offender", {
  d <- student_t_distrib()
  expect_invisible(check_theta_bounds(d, list(mu = 0, sigma = 1, nu = 5)))

  err <- tryCatch(
    check_theta_bounds(d, list(mu = 0, sigma = -1, nu = -2)),
    error = function(e) conditionMessage(e)
  )
  expect_match(err, "sigma")
  expect_match(err, "nu")
  expect_match(err, "student t")
})

test_that("vectorized parameters are validated element-wise", {
  d <- gaussian_distrib()
  expect_silent(distrib_pdf(d, c(0, 1), list(mu = c(0, 1), sigma = c(1, 2))))
  expect_error(
    distrib_pdf(d, c(0, 1), list(mu = c(0, 1), sigma = c(1, -2))),
    "outside its domain"
  )
})
