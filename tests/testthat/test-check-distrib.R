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

test_that("a check whose statistic has no value is reported as not run", {
  # The generalized Pareto's expected information does not exist for
  # xi < -1/2, so that row has no statistic to judge. It used to be a FAIL on
  # a correct family, in five runs of five; it is listed as not run instead,
  # and every row that is in the table still says OK or FAIL.
  set.seed(11)
  res <- check_distrib(gpd_distrib(), list(sigma = 1, xi = -0.7), n = 40,
                       nsim = 2e4, orders = 1:2, verbose = FALSE)
  expect_true(all(res$status %in% c("OK", "FAIL")))
  expect_false("expected information vs Monte Carlo" %in% res$check)
  sk <- attr(res, "skipped")
  expect_named(sk, c("check", "reason"))
  expect_equal(sk$check, "expected information vs Monte Carlo")
  expect_match(sk$reason, "the statistic is NA")
  # measured at three seeds, every row that is emitted passes
  expect_true(all(res$status == "OK"),
    label = paste("gpd:", paste(res$check[res$status != "OK"], collapse = "; ")))

  out <- utils::capture.output(
    check_distrib(gpd_distrib(), list(sigma = 1, xi = -0.7), n = 20,
                  nsim = 1e4, orders = 1, verbose = TRUE))
  expect_true(any(grepl("expected information vs Monte Carlo +not run", out)))
  expect_true(any(grepl("1 check not run", out)))

  # a run with nothing skipped carries no attribute
  set.seed(2)
  ok <- check_distrib(gaussian1_distrib(), list(mu = 0, sigma = 1), n = 30,
                      nsim = 2e4, orders = 1:2, verbose = FALSE)
  expect_null(attr(ok, "skipped"))
})

test_that("an infinite statistic is a defect and stays a failure", {
  # Only NaN and NA mean there is nothing to judge. A component that overflows
  # where its reference does not is exactly what the check exists to catch, so
  # reading an infinite statistic as not run would hide it.
  BadInf <- S7::new_class("BadInfScore", parent = continuous_distrib, package = NULL)
  S7::method(distrib_pdf, BadInf) <- function(distrib, y, theta, log = FALSE) {
    stats::dnorm(y, theta[[1]], theta[[2]], log = log)
  }
  S7::method(distrib_gradient, BadInf) <- function(distrib, y, theta,
                                                   scale = c("parameter", "link"), ...) {
    r <- y - theta[[1]]
    g <- list(mu = r / theta[[2]]^2, sigma = (r^2 - theta[[2]]^2) / theta[[2]]^3)
    g$mu[1] <- -Inf
    g
  }
  bad <- BadInf(
    distrib_name = "bad inf", dimension = "univariate", bounds = c(-Inf, Inf),
    params = c("mu", "sigma"), params_interpretation = c(mu = "m", sigma = "s"),
    n_params = 2, params_bounds = list(mu = c(-Inf, Inf), sigma = c(0, Inf)),
    link_params = list(mu = linkfunctions7::identity_link(),
                       sigma = linkfunctions7::log_link()))
  set.seed(3)
  res <- check_distrib(bad, list(mu = 0, sigma = 1), n = 30, nsim = 2e4,
                       orders = 1, verbose = FALSE)
  row <- res[res$check == "gradient vs finite differences", ]
  expect_equal(row$status, "FAIL")
  expect_equal(row$statistic, Inf)
})

test_that("a Monte Carlo draw exactly on a bound is left out and counted", {
  # The Weibull's score in its shape carries log(y), which is -Inf at y = 0, so
  # a single draw landing exactly on the lower bound made the whole
  # expected-information row NaN. Three such draws are put among the large
  # sample; the derivative checks draw fewer than a thousand and see none.
  wb <- weibull1_distrib()
  W0 <- S7::new_class("WeibullOnBound", parent = S7::S7_class(wb), package = NULL)
  S7::method(distrib_rng, W0) <- function(distrib, n, theta, ...) {
    y <- stats::rweibull(n, shape = theta[[2]], scale = theta[[1]])
    if (n >= 1000) y[1:3] <- 0
    y
  }
  w0 <- do.call(W0, S7::props(wb))
  set.seed(21)
  res <- check_distrib(w0, list(mu = 2, sigma = 1.5), n = 40, nsim = 2e4,
                       orders = 1:2, verbose = FALSE)
  row <- res[res$check == "expected information vs Monte Carlo", ]
  expect_equal(nrow(row), 1L)
  expect_equal(row$status, "OK")
  expect_match(row$detail, "3 of 20000 draws fell exactly on a bound")
  expect_null(attr(res, "skipped"))

  # A non-finite score at an INTERIOR point is not a bound, and still leaves
  # the row without a value: a gaussian whose score is NaN at one interior
  # point, drawn three times.
  gs <- gaussian1_distrib()
  G0 <- S7::new_class("GaussianNaNScore", parent = S7::S7_class(gs), package = NULL)
  S7::method(distrib_rng, G0) <- function(distrib, n, theta, ...) {
    y <- stats::rnorm(n, theta[[1]], theta[[2]])
    if (n >= 1000) y[1:3] <- theta[[1]] + 0.5
    y
  }
  S7::method(distrib_gradient, G0) <- function(distrib, y, theta,
                                               scale = c("parameter", "link"), ...) {
    r <- y - theta[[1]]
    g <- list(mu = r / theta[[2]]^2, sigma = (r^2 - theta[[2]]^2) / theta[[2]]^3)
    g$mu[r == 0.5] <- NaN
    g
  }
  g0 <- do.call(G0, S7::props(gs))
  set.seed(21)
  res2 <- check_distrib(g0, list(mu = 0, sigma = 1), n = 40, nsim = 2e4,
                        orders = 1:2, verbose = FALSE)
  expect_false("expected information vs Monte Carlo" %in% res2$check)
  expect_equal(attr(res2, "skipped")$check, "expected information vs Monte Carlo")
})

test_that("the response reference chooses its step near a bound", {
  # A gamma with shape 1/2 carries (a - 1) log y near zero. A central
  # difference whose step is cut to 0.49 of the distance to the bound has a
  # relative error of (h/d)^2/3 once the cut binds, 9.4e-02 on the first
  # derivative, which halving the step cannot expose; the quotient that also
  # tries a step scaled on the distance reads the closed form to 2e-10.
  d <- gamma1_distrib()
  th <- list(mu = 2, phi = 2)
  y <- c(1e-5, 1e-6, 1e-7)
  eps <- .Machine$double.eps
  lp <- function(v) distrib_pdf(d, v, th, log = TRUE)
  q1 <- function(h) (lp(y + h) - lp(y - h)) / (2 * h)
  g <- distrib_grad_y(d, y, th)
  cut <- abs(q1(fd_steps_y(y, d@bounds, eps^(1 / 3))) / g - 1)
  chosen <- abs(fd_stable_quotient(q1, y, d@bounds, eps^(1 / 3)) / g - 1)
  expect_true(all(cut > 1e-2))
  expect_true(all(chosen < 1e-8))

  # where no bound is finite the two steps coincide, and the cut quotient
  # comes back as it is without the three further differences
  y2 <- c(-3, 0.2, 5)
  gs <- gaussian1_distrib()
  lq <- function(v) distrib_pdf(gs, v, list(mu = 0, sigma = 1), log = TRUE)
  calls <- 0L
  q2 <- function(h) {
    calls <<- calls + 1L
    (lq(y2 + h) - lq(y2 - h)) / (2 * h)
  }
  expect_identical(fd_stable_quotient(q2, y2, c(-Inf, Inf), eps^(1 / 3)),
                   q2(fd_steps_y(y2, c(-Inf, Inf), eps^(1 / 3))))
  expect_equal(calls, 2L)
})

test_that("the cdf row keeps its step inside the support", {
  # A gamma with shape 1/5 puts the lowest grid point 3.3e-05 above zero. The
  # old step of 1e-5 max(1, |x|) is not kept inside the support, and its
  # difference is out by 2.4e-02 against a threshold of 1e-4.
  d <- gamma1_distrib()
  th <- list(mu = 1, phi = 5)
  grid <- distrib_quantile(d, seq(0.1, 0.9, length.out = 15), th)
  f <- distrib_pdf(d, grid, th)
  old_h <- pmax(abs(grid), 1) * 1e-5
  old <- (distrib_cdf(d, grid + old_h, th) - distrib_cdf(d, grid - old_h, th)) / (2 * old_h)
  expect_gt(max(abs(old - f) / pmax(1, abs(f))), 1e-4)

  set.seed(1)
  res <- check_distrib(d, th, n = 100, nsim = 2e4, orders = 1:2, verbose = FALSE)
  expect_equal(res$status[res$check == "cdf agrees with the density"], "OK")
  expect_equal(res$status[res$check == "response derivatives vs finite differences"], "OK")

  # and a 5% error in the response score near the bound is still caught
  Bad <- S7::new_class("Gamma1BadResponseScore", parent = S7::S7_class(d),
                       package = NULL)
  S7::method(distrib_grad_y, Bad) <- function(distrib, y, theta, ...) {
    a <- 1 / theta[[2]]
    1.05 * ((a - 1) / y - a / theta[[1]])
  }
  bad <- do.call(Bad, S7::props(d))
  set.seed(1)
  res_bad <- check_distrib(bad, list(mu = 3, phi = 2), n = 100, nsim = 2e4,
                           orders = 1:2, verbose = FALSE)
  expect_equal(res_bad$status[res_bad$check == "response derivatives vs finite differences"],
               "FAIL")
})

test_that("near a non-zero bound the reference divides by the steps taken", {
  # beta1 at mu 0.5, phi 0.5 is singular at 1, where the spacing of doubles is
  # absolute. A step that is not a whole number of units of that spacing is
  # not the step the evaluation points lie at, and one below a unit rounds to
  # zero; the reference on the steps taken, with no step below two units,
  # reads the closed form where the nominal quotients sit on a plateau.
  d <- beta1_distrib()
  th <- list(mu = 0.5, phi = 0.5)
  y <- 1 - c(1e-10, 1e-12, 1e-13)
  b <- d@bounds
  eps <- .Machine$double.eps
  lp <- function(v) distrib_pdf(d, v, th, log = TRUE)
  l0 <- lp(y)
  g <- distrib_grad_y(d, y, th)
  hh <- distrib_hess_y(d, y, th)
  rel <- function(a, e) abs(a - e) / pmax(1, abs(e))
  eg <- rel(fd_stable_quotient(function(h) fd_first_taken(lp, y, h), y, b, eps^(1 / 3)), g)
  eh <- rel(fd_stable_quotient(function(h) fd_second_taken(lp, y, h, l0), y, b, eps^(1 / 4)), hh)
  expect_true(all(eg < 1e-4))
  expect_true(all(eh < 1e-4))
  ng <- rel(fd_stable_quotient(function(h) (lp(y + h) - lp(y - h)) / (2 * h), y, b, eps^(1 / 3)), g)
  nh <- rel(fd_stable_quotient(function(h) (lp(y + h) - 2 * l0 + lp(y - h)) / h^2, y, b, eps^(1 / 4)), hh)
  expect_gt(max(ng, nh), 1e-2)

  # where the steps are exact the two quotients are the uniform ones, bit for bit
  x <- c(0.3, 2, 7)
  f <- function(v) sin(v) + v^3
  h <- 2^-10
  expect_identical(fd_first_taken(f, x, h), (f(x + h) - f(x - h)) / (2 * h))
  expect_identical(fd_second_taken(f, x, h), (f(x + h) - 2 * f(x) + f(x - h)) / h^2)
})

test_that("a draw in the last places of a non-zero bound is left out of the response row", {
  # Three draws are put within 128 |b| eps of 1, at 3, 16 and 200 units of the
  # spacing there, where no central difference in double precision compares a
  # derivative: at 16 units the reference is out by more than the tolerance.
  d <- beta1_distrib()
  B1 <- S7::new_class("Beta1NearOne", parent = S7::S7_class(d), package = NULL)
  S7::method(distrib_rng, B1) <- function(distrib, n, theta, ...) {
    y <- stats::rbeta(n, theta[[1]] * theta[[2]], (1 - theta[[1]]) * theta[[2]])
    if (n < 1000) y[1:3] <- 1 - c(3, 16, 200) * 2^-53
    y
  }
  b1 <- do.call(B1, S7::props(d))
  th <- list(mu = 0.5, phi = 0.5)
  set.seed(4)
  res <- check_distrib(b1, th, n = 60, nsim = 2e4, orders = 1:2, verbose = FALSE)
  row <- res[res$check == "response derivatives vs finite differences", ]
  expect_equal(row$status, "OK")
  expect_match(row$detail, "3 of 60 draws lay within 128 |b| eps of a non-zero bound",
               fixed = TRUE)

  yy <- 1 - 16 * 2^-53
  lp <- function(v) distrib_pdf(d, v, th, log = TRUE)
  ref <- fd_stable_quotient(function(h) fd_first_taken(lp, yy, h), yy, d@bounds,
                            .Machine$double.eps^(1 / 3))
  an <- distrib_grad_y(d, yy, th)
  expect_gt(abs(ref - an) / abs(an), 1e-3)
})
