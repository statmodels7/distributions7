# The second parametrisations whose derivatives are elementary at every order:
# the inverse gaussian in its shape, and the beta and beta-binomial in their
# canonical shapes. Each is checked against its twin, which shares no code.

test_that("invgauss2 is invgauss1 with lambda = 1/phi", {
  set.seed(2)
  y <- statmod::rinvgauss(300, mean = 2, dispersion = 1 / 3)
  d2 <- invgauss2_distrib(); t2 <- list(mu = 2, lambda = 3)
  d1 <- invgauss1_distrib(); t1 <- list(mu = 2, phi = 1 / 3)

  expect_identical(distrib_pdf(d2, y, t2), distrib_pdf(d1, y, t1))
  expect_identical(distrib_cdf(d2, y, t2), distrib_cdf(d1, y, t1))
  expect_equal(variance(d2, t2), 2^3 / 3)
  expect_equal(skewness(d2, t2), skewness(d1, t1))

  f2 <- fit_distrib(d2, y)
  f1 <- fit_distrib(d1, y)
  expect_true(f2@converged)
  expect_equal(as.numeric(logLik(f2)), as.numeric(logLik(f1)), tolerance = 1e-8)
  expect_equal(unname(coef(f2)[2]), unname(1 / coef(f1)[2]),
               tolerance = sqrt(eval(formals(optimizers7::crit_grad)$tol)))
})


test_that("beta2 is beta1 in the shapes", {
  set.seed(4)
  y <- rbeta(300, 2, 5)
  d2 <- beta2_distrib(); t2 <- list(alpha = 2, beta = 5)
  d1 <- beta1_distrib(); t1 <- list(mu = 2 / 7, phi = 7)

  expect_identical(distrib_pdf(d2, y, t2), distrib_pdf(d1, y, t1))
  expect_identical(distrib_cdf(d2, y, t2), distrib_cdf(d1, y, t1))
  expect_equal(mean(d2, t2), 2 / 7)
  expect_equal(variance(d2, t2), variance(d1, t1))
  expect_equal(kurtosis(d2, t2), kurtosis(d1, t1))

  f2 <- fit_distrib(d2, y)
  f1 <- fit_distrib(d1, y)
  expect_true(f2@converged)
  expect_equal(as.numeric(logLik(f2)), as.numeric(logLik(f1)), tolerance = 1e-8)
  expect_equal(sum(coef(f2)), unname(coef(f1)[2]),
               tolerance = sqrt(eval(formals(optimizers7::crit_grad)$tol)))
})


test_that("the beta in its shapes has data-free derivatives past the first", {
  # The data enter the log-density only through log y and log(1-y), both linear
  # in the parameters, so every derivative beyond the first is a constant. The
  # expected Hessian is then the observed one exactly, not to a tolerance.
  d <- beta2_distrib()
  th <- list(alpha = 2, beta = 5)
  set.seed(4)
  y <- rbeta(50, 2, 5)

  h <- distrib_hessian(d, y, th)
  for (nm in names(h)) {
    expect_length(unique(h[[nm]]), 1L)
  }
  eh <- distrib_expected_hessian(d, y[1], th)
  expect_identical(unname(unlist(eh)),
                   unname(vapply(h, function(v) v[1], numeric(1))))

  # and the printed second derivatives are the polygamma differences
  a <- 2; b <- 5
  expect_equal(h[["alpha_alpha"]][1], trigamma(a + b) - trigamma(a))
  expect_equal(h[["alpha_beta"]][1], trigamma(a + b))
  expect_equal(h[["beta_beta"]][1], trigamma(a + b) - trigamma(b))
})


test_that("betabinom2 is betabinom1 in the shapes", {
  d2 <- betabinom2_distrib(size = 10); t2 <- list(alpha = 2, beta = 3)
  d1 <- betabinom1_distrib(size = 10); t1 <- list(mu = 0.4, sigma = 0.2)

  expect_equal(distrib_pdf(d2, 0:10, t2), distrib_pdf(d1, 0:10, t1),
               tolerance = 1e-14)
  expect_equal(sum(distrib_pdf(d2, 0:10, t2)), 1, tolerance = 1e-14)
  expect_equal(mean(d2, t2), 10 * 0.4)
  expect_equal(variance(d2, t2), variance(d1, t1))
  expect_equal(skewness(d2, t2), skewness(d1, t1))

  # off the support there is no mass, and no complaint: the support is tested
  # before lchoose(), which warns on a non-integer count
  expect_silent(f <- distrib_pdf(d2, c(-1, 11, 2.5), t2))
  expect_identical(f, c(0, 0, 0))
})


test_that("the beta-binomial expectation is an exact sum over the support", {
  # A finite support means the expectation is a weighted sum of eleven terms,
  # so it is exact rather than approximated: the closed-form expected Hessian
  # and the sum of the observed one against the mass agree to rounding.
  d <- betabinom2_distrib(size = 10)
  th <- list(alpha = 2, beta = 3)
  supp <- 0:10
  w <- distrib_pdf(d, supp, th)
  hs <- distrib_hessian(d, supp, th)
  eh <- distrib_expected_hessian(d, supp[1], th)
  for (nm in names(eh)) {
    expect_equal(eh[[nm]][1], sum(w * hs[[nm]]), tolerance = 1e-12, label = nm)
  }
})


test_that("each of the three passes the whole validator", {
  cases <- list(
    list(invgauss2_distrib(), list(mu = 2, lambda = 3), 13L),
    list(beta2_distrib(), list(alpha = 2, beta = 5), 13L),
    list(betabinom2_distrib(size = 10), list(alpha = 2, beta = 3), 12L)
  )
  for (cs in cases) {
    set.seed(1)
    res <- check_distrib(cs[[1]], cs[[2]], verbose = FALSE)
    expect_identical(nrow(res), cs[[3]], label = cs[[1]]@distrib_name)
    expect_true(all(res$status == "OK"),
                info = paste(cs[[1]]@distrib_name,
                             paste(res$check[res$status != "OK"], collapse = ", ")))
  }
})


test_that("a broken gradient is still caught", {
  good <- S7::method(distrib_gradient, InvGauss2Distrib)
  on.exit(S7::method(distrib_gradient, InvGauss2Distrib) <- good, add = TRUE)
  suppressMessages(
    S7::method(distrib_gradient, InvGauss2Distrib) <- function(distrib, y, theta,
                                                               scale = c("parameter", "link"), ...) {
      g <- good(distrib, y, theta, scale = scale, ...)
      g[[1]] <- g[[1]] * 1.05
      g
    }
  )
  set.seed(1)
  res <- check_distrib(invgauss2_distrib(), list(mu = 2, lambda = 3),
                       verbose = FALSE)
  expect_equal(res$status[res$check == "gradient vs finite differences"], "FAIL")
})


test_that("the third and fourth derivatives are the analytical ones", {
  skip_if_not_installed("numDeriv")
  cases <- list(
    list(invgauss2_distrib(), list(mu = 2, lambda = 3)),
    list(beta2_distrib(), list(alpha = 2, beta = 5)),
    list(betabinom2_distrib(size = 10), list(alpha = 2, beta = 3))
  )
  for (cs in cases) {
    d <- cs[[1]]
    th <- cs[[2]]
    set.seed(3)
    y <- distrib_rng(d, 25, th)
    p <- d@params
    d3 <- distrib_deriv3(d, y, th)
    d4 <- distrib_deriv4(d, y, th)

    for (hn in hess_names(p)) {
      for (k in seq_along(p)) {
        key <- paste(sort(c(hess_pairs(p)[[hn]], p[k])), collapse = "_")
        if (is.null(d3[[key]])) next
        f <- function(v) {
          q <- unlist(th)
          q[k] <- v
          sum(distrib_hessian(d, y, as.list(stats::setNames(q, p)))[[hn]])
        }
        expect_equal(sum(d3[[key]]), numDeriv::grad(f, th[[k]]),
                     tolerance = 1e-6, label = paste(d@distrib_name, key))
      }
    }
    for (n3 in names(d3)) {
      for (k in seq_along(p)) {
        key <- paste(sort(c(strsplit(n3, "_")[[1]], p[k])), collapse = "_")
        if (is.null(d4[[key]])) next
        f <- function(v) {
          q <- unlist(th)
          q[k] <- v
          sum(distrib_deriv3(d, y, as.list(stats::setNames(q, p)))[[n3]])
        }
        expect_equal(sum(d4[[key]]), numDeriv::grad(f, th[[k]]),
                     tolerance = 1e-6, label = paste(d@distrib_name, key))
      }
    }
  }
})


test_that("the constructor of betabinom2 validates its size", {
  expect_error(betabinom2_distrib(size = 0), "positive integer")
  expect_error(betabinom2_distrib(size = 2.5), "positive integer")
})
