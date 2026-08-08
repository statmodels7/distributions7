# check_distrib(): the user-facing numerical validator.

test_that("built-in distributions pass every check", {
  set.seed(1)
  cases <- list(
    list(d = gaussian1_distrib(), th = list(mu = 1.5, sigma = 2)),
    list(d = poisson_distrib(),  th = list(mu = 4)),
    list(d = laplace_distrib(),  th = list(mu = 1, sigma = 2))
  )
  for (case in cases) {
    res <- check_distrib(case$d, theta = case$th, n = 40, nsim = 5e4,
                         orders = 1:2, verbose = FALSE)
    expect_true(all(res$status == "OK"),
                label = paste(case$d@distrib_name, "-",
                              paste(res$check[res$status != "OK"], collapse = "; ")))
  }
})

test_that("check_distrib is reproducible and survives a draw on a kink", {
  # It used to call set.seed(NULL), discarding whatever seed the caller had set,
  # so two runs never agreed and a failing check could not be reproduced.
  set.seed(42)
  a <- check_distrib(laplace_distrib(), list(mu = 1, sigma = 2), n = 40, nsim = 2e4,
                     orders = 1:2, verbose = FALSE)
  set.seed(42)
  b <- check_distrib(laplace_distrib(), list(mu = 1, sigma = 2), n = 40, nsim = 2e4,
                     orders = 1:2, verbose = FALSE)
  expect_equal(a$statistic, b$statistic)

  # The Laplace has no derivative at y = mu. An observation landing within a
  # finite-difference step of it makes the FD *reference* wrong, not the
  # analytical value, and the check used to report a failure for correct code --
  # rarely and unpredictably, since the draws are random. This forces the worst
  # case: one observation exactly on the kink and one a nanometre away.
  Kinked <- S7::new_class("KinkedLap", parent = continuous_distrib, package = NULL)
  S7::method(distrib_pdf, Kinked) <- function(distrib, y, theta, log = FALSE) {
    ld <- -log(2 * theta[[2]]) - abs(y - theta[[1]]) / theta[[2]]
    if (log) ld else exp(ld)
  }
  S7::method(distrib_gradient, Kinked) <- function(distrib, y, theta,
                                                   scale = c("parameter", "link"), ...) {
    r <- y - theta[[1]]
    list(mu = sign(r) / theta[[2]], b = (abs(r) / theta[[2]] - 1) / theta[[2]])
  }
  S7::method(distrib_quantile, Kinked) <- function(distrib, p, theta,
                                                   lower.tail = TRUE, log.p = FALSE) {
    mu <- theta[[1]]; bb <- theta[[2]]
    mu - bb * sign(p - 0.5) * log(1 - 2 * abs(p - 0.5))
  }
  S7::method(distrib_rng, Kinked) <- function(distrib, n, theta) {
    y <- distrib_quantile(distrib, stats::runif(n), theta)
    y[1] <- theta[[1]]           # exactly on the kink
    y[2] <- theta[[1]] + 1e-9    # inside any finite-difference step
    y
  }
  kink <- Kinked(
    distrib_name = "kinked", dimension = "univariate", bounds = c(-Inf, Inf),
    params = c("mu", "b"), params_interpretation = c(mu = "loc", b = "scale"),
    n_params = 2, params_bounds = list(mu = c(-Inf, Inf), b = c(0, Inf)),
    link_params = list(
      mu = linkfunctions7::identity_link(),
      b  = linkfunctions7::log_link()
    ),
    params_smooth = c(mu = FALSE, b = TRUE)
  )

  set.seed(3)
  res <- check_distrib(kink, list(mu = 1, b = 2), n = 40, nsim = 2e4,
                       orders = 1:2, verbose = FALSE)
  expect_true(all(res$status == "OK"),
    label = paste("kinked:", paste(res$check[res$status != "OK"], collapse = "; ")))

  # and the guard must not blunt the check: the same kinked distribution with a
  # gradient that is 5% wrong is still caught
  KinkedBad <- S7::new_class("KinkedLapBad", parent = Kinked, package = NULL)
  S7::method(distrib_gradient, KinkedBad) <- function(distrib, y, theta,
                                                      scale = c("parameter", "link"), ...) {
    r <- y - theta[[1]]
    list(mu = 1.05 * sign(r) / theta[[2]], b = (abs(r) / theta[[2]] - 1) / theta[[2]])
  }
  kink_bad <- KinkedBad(
    distrib_name = "kinked bad", dimension = "univariate", bounds = c(-Inf, Inf),
    params = c("mu", "b"), params_interpretation = c(mu = "loc", b = "scale"),
    n_params = 2, params_bounds = list(mu = c(-Inf, Inf), b = c(0, Inf)),
    link_params = list(
      mu = linkfunctions7::identity_link(),
      b  = linkfunctions7::log_link()
    ),
    params_smooth = c(mu = FALSE, b = TRUE)
  )
  set.seed(5)
  bad <- check_distrib(kink_bad, list(mu = 1, b = 2), n = 40, nsim = 2e4,
                       orders = 1:2, verbose = FALSE)
  expect_true(any(grepl("gradient", bad$check[bad$status != "OK"])))
})

test_that("check_distrib catches a defect injected into each component", {
  # A validator that never fails is worth nothing. Break one thing at a time in
  # an otherwise-correct gaussian and confirm the corresponding check goes red.
  build <- function(nm, methods) {
    cls <- S7::new_class(nm, parent = continuous_distrib, package = NULL)
    S7::method(distrib_pdf, cls) <- function(distrib, y, theta, log = FALSE) {
      stats::dnorm(y, theta[[1]], theta[[2]], log = log)
    }
    for (g in names(methods)) {
      gen <- get(g)                      # S7's method<- mutates the generic
      S7::method(gen, cls) <- methods[[g]]
    }
    cls(
      distrib_name = nm, dimension = "univariate", bounds = c(-Inf, Inf),
      params = c("mu", "sigma"), params_interpretation = c(mu = "m", sigma = "s"),
      n_params = 2, params_bounds = list(mu = c(-Inf, Inf), sigma = c(0, Inf)),
      link_params = list(
        mu = linkfunctions7::identity_link(),
        sigma = linkfunctions7::log_link()
      )
    )
  }

  defects <- list(
    gradient = list(
      m = list(distrib_gradient = function(distrib, y, theta, scale = c("parameter", "link"), ...) {
        r <- y - theta[[1]]
        list(mu = 1.05 * r / theta[[2]]^2, sigma = (r^2 / theta[[2]]^2 - 1) / theta[[2]])
      }),
      hits = "gradient"),
    # A cdf shifted by a constant satisfies every other check: it stays in [0,1],
    # it is non-decreasing, and the quantile round-trip cannot see it because the
    # numerical quantile is derived from that same cdf. Only the comparison with
    # the density pins it down.
    cdf = list(
      m = list(distrib_cdf = function(distrib, q, theta, lower.tail = TRUE, log.p = FALSE) {
        stats::pnorm(q, theta[[1]] + 0.05, theta[[2]], lower.tail = lower.tail, log.p = log.p)
      }),
      hits = "cdf agrees with the density"),
    rng = list(
      m = list(distrib_rng = function(distrib, n, theta) {
        stats::rnorm(n, theta[[1]] + 0.25 * theta[[2]], theta[[2]])
      }),
      hits = "rng"),
    response = list(
      m = list(distrib_grad_y = function(distrib, y, theta) {
        -1.4 * (y - theta[[1]]) / theta[[2]]^2
      }),
      hits = "response")
  )

  for (nm in names(defects)) {
    set.seed(51)
    res <- check_distrib(build(paste0("Broken", nm), defects[[nm]]$m),
                         list(mu = 1.5, sigma = 2), n = 60, nsim = 4e4, verbose = FALSE)
    failed <- res$check[res$status != "OK"]
    expect_true(any(grepl(defects[[nm]]$hits, failed)),
      label = sprintf("defect '%s' should be caught (failed: %s)", nm,
                      paste(failed, collapse = "; ")))
  }

  # and the same construction with nothing broken must pass everything
  set.seed(51)
  ok <- check_distrib(build("BrokenNone", list()), list(mu = 1.5, sigma = 2),
                      n = 60, nsim = 4e4, verbose = FALSE)
  expect_true(all(ok$status == "OK"),
    label = paste("reference:", paste(ok$check[ok$status != "OK"], collapse = "; ")))
})

test_that("a deterministic score product does not fake an information mismatch", {
  # Regression: the Laplace score for mu is sign(y - mu)/b, so its square is
  # exactly 1/b^2 on every draw. The Monte Carlo standard error is then pure
  # floating-point dust, and standardizing the (essentially exact) difference by
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
  res <- check_distrib(gaussian1_distrib(), theta = list(mu = 0, sigma = 1),
                       n = 30, nsim = 2e4, orders = 1:2, verbose = FALSE)
  expect_s3_class(res, "data.frame")
  expect_named(res, c("check", "status", "statistic", "detail"))
  expect_true(all(res$status %in% c("OK", "FAIL")))
  expect_true(any(grepl("density", res$check)))
  expect_true(any(grepl("link-scale", res$check)))

  out <- utils::capture.output(
    check_distrib(gaussian1_distrib(), theta = list(mu = 0, sigma = 1),
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
  res <- check_distrib(gaussian1_distrib(), n = 20, nsim = 1e4,
                       orders = 1, verbose = FALSE)
  expect_s3_class(res, "data.frame")
  expect_gt(nrow(res), 0)
})
