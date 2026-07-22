# check_distrib(): the user-facing numerical validator.

test_that("built-in distributions pass every check", {
  set.seed(1)
  cases <- list(
    list(d = gaussian_distrib(), th = list(mu = 1.5, sigma = 2)),
    list(d = poisson_distrib(),  th = list(mu = 4)),
    list(d = laplace_distrib(),  th = list(mu = 1, b = 2))
  )
  for (case in cases) {
    res <- check_distrib(case$d, theta = case$th, n = 40, nsim = 5e4,
                         orders = 1:2, verbose = FALSE)
    expect_true(all(res$status == "OK"),
                label = paste(case$d@distrib_name, "-",
                              paste(res$check[res$status != "OK"], collapse = "; ")))
  }
})

test_that("a deterministic score product does not fake an information mismatch", {
  # Regression: the Laplace score for mu is sign(y - mu)/b, so its square is
  # exactly 1/b^2 on every draw. The Monte Carlo standard error is then pure
  # floating-point dust, and standardising the (essentially exact) difference by
  # it used to report a z around 1500 -- a failure raised precisely because the
  # two sides agreed to eleven decimals.
  CheckLap <- S7::new_class("CheckLap", parent = continuous_distrib, package = NULL)
  S7::method(distrib_pdf, CheckLap) <- function(distrib, y, theta, log = FALSE) {
    ld <- -log(2 * theta[[2]]) - abs(y - theta[[1]]) / theta[[2]]
    if (log) ld else exp(ld)
  }
  d <- CheckLap(
    distrib_name = "check laplace", dimension = "univariate", bounds = c(-Inf, Inf),
    params = c("mu", "b"), params_interpretation = c(mu = "loc", b = "scale"),
    n_params = 2, params_bounds = list(mu = c(-Inf, Inf), b = c(0, Inf)),
    link_params = list(
      mu = linkfunctions7::identity_link(),
      b = linkfunctions7::log_link()
    ),
    params_smooth = c(mu = FALSE, b = TRUE)
  )

  for (s in 1:4) {
    set.seed(s)
    res <- check_distrib(d, theta = list(mu = 1, b = 2), n = 40, nsim = 5e4,
                         orders = 1:2, verbose = FALSE)
    row <- res[res$check == "expected information vs Monte Carlo", ]
    expect_equal(row$status, "OK", label = sprintf("seed %d (z = %.4g)", s, row$statistic))
  }
})

test_that("check_distrib returns a tidy data frame and prints a report", {
  set.seed(2)
  res <- check_distrib(gaussian_distrib(), theta = list(mu = 0, sigma = 1),
                       n = 30, nsim = 2e4, orders = 1:2, verbose = FALSE)
  expect_s3_class(res, "data.frame")
  expect_named(res, c("check", "status", "statistic", "detail"))
  expect_true(all(res$status %in% c("OK", "FAIL")))
  expect_true(any(grepl("density", res$check)))
  expect_true(any(grepl("link-scale", res$check)))

  out <- utils::capture.output(
    check_distrib(gaussian_distrib(), theta = list(mu = 0, sigma = 1),
                  n = 20, nsim = 1e4, orders = 1, verbose = TRUE)
  )
  expect_true(any(grepl("Distribution:", out)))
  expect_true(any(grepl("checks passed", out)))
})

test_that("check_distrib detects a wrong analytical gradient", {
  # a distribution whose density is right but whose score is deliberately wrong
  BadGrad <- S7::new_class("BadGrad", parent = continuous_distrib, package = NULL)
  S7::method(distrib_pdf, BadGrad) <- function(distrib, y, theta, log = FALSE) {
    stats::dnorm(y, theta[[1]], theta[[2]], log = log)
  }
  S7::method(distrib_gradient, BadGrad) <- function(distrib, y, theta,
                                                    scale = c("parameter", "link"), ...) {
    # correct value would be (y - mu)/sigma^2
    list(mu = 2 * (y - theta[[1]]) / theta[[2]]^2,
         sigma = ((y - theta[[1]])^2 - theta[[2]]^2) / theta[[2]]^3)
  }
  bad <- BadGrad(
    distrib_name = "bad gradient", dimension = "univariate", bounds = c(-Inf, Inf),
    params = c("mu", "sigma"), params_interpretation = c(mu = "m", sigma = "s"),
    n_params = 2, params_bounds = list(mu = c(-Inf, Inf), sigma = c(0, Inf)),
    link_params = list(mu = linkfunctions7::identity_link(),
                       sigma = linkfunctions7::log_link())
  )

  set.seed(3)
  res <- check_distrib(bad, theta = list(mu = 0, sigma = 1), n = 30, nsim = 2e4,
                       orders = 1, verbose = FALSE)
  grad_row <- res[grepl("^gradient", res$check), ]
  expect_equal(nrow(grad_row), 1L)
  expect_equal(grad_row$status, "FAIL")

  # the density itself is fine, so those checks still pass
  expect_equal(res$status[res$check == "density integrates to 1"], "OK")
})

test_that("check_distrib draws its own theta when none is supplied", {
  set.seed(4)
  res <- check_distrib(gaussian_distrib(), n = 20, nsim = 1e4,
                       orders = 1, verbose = FALSE)
  expect_s3_class(res, "data.frame")
  expect_gt(nrow(res), 0)
})
