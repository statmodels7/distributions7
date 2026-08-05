# reparametrize() and the families built on it.
#
# The strongest check available here is that the general mechanism reproduces
# the families written out by hand, which share no code with it: gaussian2,
# gaussian3, gamma1, beta2 and invgauss2 each carry their own class, their own
# kernels and their own moments. Agreement to machine precision then needs no
# tolerance to be argued for.

il <- function() linkfunctions7::identity_link()
ll <- function() linkfunctions7::log_link()

reparam_gaussian2 <- function() {
  reparametrize(
    gaussian1_distrib(),
    map = function(psi) list(mu = psi$mu, sigma = sqrt(psi$sigma2)),
    params = c("mu", "sigma2"),
    bounds = list(mu = c(-Inf, Inf), sigma2 = c(0, Inf)),
    links = list(mu = il(), sigma2 = ll())
  )
}

reparam_beta2 <- function() {
  reparametrize(
    beta1_distrib(),
    map = function(psi) list(mu = psi$alpha / (psi$alpha + psi$beta),
                             phi = psi$alpha + psi$beta),
    params = c("alpha", "beta"),
    bounds = list(alpha = c(0, Inf), beta = c(0, Inf)),
    links = list(alpha = ll(), beta = ll())
  )
}

# Every derivative of both objects at one point, keyed alike.
all_derivs <- function(d, y, theta, expected = FALSE) {
  if (expected) {
    list(distrib_expected_hessian(d, y, theta),
         distrib_deriv3(d, y, theta, expected = TRUE),
         distrib_deriv4(d, y, theta, expected = TRUE))
  } else {
    list(distrib_gradient(d, y, theta), distrib_hessian(d, y, theta),
         distrib_deriv3(d, y, theta), distrib_deriv4(d, y, theta))
  }
}

worst_gap <- function(a, b) {
  max(vapply(seq_along(a), function(k) {
    max(vapply(names(b[[k]]), function(nm) {
      max(abs(a[[k]][[nm]] - b[[k]][[nm]]))
    }, numeric(1)))
  }, numeric(1)))
}


test_that("the mechanism reproduces the families written by hand", {
  set.seed(7)
  y <- rnorm(40, 1.2, sqrt(3.5))
  th <- list(mu = 1.2, sigma2 = 3.5)
  r <- reparam_gaussian2()
  h <- gaussian2_distrib()

  expect_identical(distrib_pdf(r, y, th), distrib_pdf(h, y, th))
  expect_lt(worst_gap(all_derivs(r, y, th), all_derivs(h, y, th)), 1e-12)

  # a map that is dense, both parent parameters moving with both new ones
  rg <- reparametrize(
    gamma2_distrib(),
    map = function(psi) list(mu = psi$mu, sigma2 = psi$phi * psi$mu^2),
    params = c("mu", "phi"),
    bounds = list(mu = c(0, Inf), phi = c(0, Inf)),
    links = list(mu = ll(), phi = ll())
  )
  hg <- gamma1_distrib()
  thg <- list(mu = 3, phi = 0.4)
  set.seed(2)
  yg <- rgamma(40, shape = 1 / 0.4, rate = 1 / (0.4 * 3))
  expect_identical(distrib_pdf(rg, yg, thg), distrib_pdf(hg, yg, thg))
  expect_lt(worst_gap(all_derivs(rg, yg, thg), all_derivs(hg, yg, thg)), 1e-10)

  # and one that is non-linear in both
  rb <- reparam_beta2()
  hb <- beta2_distrib()
  thb <- list(alpha = 2, beta = 5)
  set.seed(3)
  yb <- rbeta(40, 2, 5)
  expect_identical(distrib_pdf(rb, yb, thb), distrib_pdf(hb, yb, thb))
  expect_lt(worst_gap(all_derivs(rb, yb, thb), all_derivs(hb, yb, thb)), 1e-12)
})


test_that("the expected derivatives are carried, not dropped", {
  # A regression. The number of blocks in a partition IS the order of the
  # parent derivative that term carries, so what drops under expectation is the
  # ONE-block term, the score. Skipping the order-block term instead removes
  # the parent's own top derivative, and the expected Hessian of a
  # reparametrized gaussian comes back 0 against 1/3.5 -- a failure no observed
  # derivative can see, since those were exact throughout.
  set.seed(7)
  y <- rnorm(40, 1.2, sqrt(3.5))
  th <- list(mu = 1.2, sigma2 = 3.5)
  r <- reparam_gaussian2()
  h <- gaussian2_distrib()

  eh <- distrib_expected_hessian(r, y, th)
  expect_equal(eh[["mu_mu"]][1], -1 / 3.5, tolerance = 1e-12)
  expect_gt(abs(eh[["mu_mu"]][1]), 0.1)

  expect_lt(worst_gap(all_derivs(r, y, th, expected = TRUE),
                      all_derivs(h, y, th, expected = TRUE)), 1e-12)
})


test_that("a reparametrized family passes the validator and fits", {
  rb <- reparam_beta2()
  thb <- list(alpha = 2, beta = 5)
  set.seed(1)
  res <- check_distrib(rb, thb, verbose = FALSE)
  expect_true(all(res$status == "OK"),
              info = paste(res$check[res$status != "OK"], collapse = ", "))

  set.seed(5)
  y <- rbeta(500, 2, 5)
  f <- fit_distrib(rb, y)
  fh <- fit_distrib(beta2_distrib(), y)
  expect_true(f@converged)
  expect_equal(as.numeric(logLik(f)), as.numeric(logLik(fh)), tolerance = 1e-8)
  expect_equal(unname(coef(f)), unname(coef(fh)),
               tolerance = sqrt(eval(formals(optimizers7::crit_grad)$tol)))
})


test_that("the moments delegate rather than falling to a quadrature", {
  r <- reparam_gaussian2()
  th <- list(mu = 1.2, sigma2 = 3.5)
  expect_equal(mean(r, th), 1.2)
  expect_equal(variance(r, th), 3.5)
  expect_identical(skewness(r, th), 0)
  expect_identical(kurtosis(r, th), 0)
})


test_that("the constructor probes inside each parameter's own interval", {
  # A regression. Taking the midpoint of a half-line gives infinity; falling
  # back to zero put the probe outside the domain, and a map like
  # sqrt((nu-2)/nu) warned "NaNs produced" at construction, naming neither the
  # parameter nor the reason.
  expect_silent(student_t2_distrib())

  # and a map that genuinely cannot be evaluated says where
  expect_error(
    reparametrize(
      gaussian1_distrib(),
      map = function(psi) list(mu = psi$mu, sigma = log(psi$sigma - 100)),
      params = c("mu", "sigma"),
      bounds = list(mu = c(-Inf, Inf), sigma = c(0, Inf)),
      links = list(mu = il(), sigma = ll())
    ),
    "sigma = 1"
  )
})


test_that("reparametrize refuses what it cannot use", {
  expect_error(reparametrize(3, map = identity, params = "a",
                             bounds = list(a = c(0, 1)), links = list(a = ll())),
               "continuous_distrib")
  expect_error(
    reparametrize(gaussian1_distrib(), map = function(psi) list(mu = psi$a),
                  params = "a", bounds = list(a = c(0, 1)),
                  links = list(a = ll())),
    "did not return 'sigma'"
  )
  expect_error(
    reparametrize(gaussian1_distrib(),
                  map = function(psi) list(mu = psi$a, sigma = psi$a),
                  params = "a", bounds = list(a = c(0, 1)),
                  links = list()),
    "one entry per new parameter"
  )
})


test_that("the four families built on it are what they say", {
  cases <- list(
    list(lognormal2_distrib(), list(mean = 3, var = 2)),
    list(weibull3_distrib(), list(mean = 4, sigma = 1.7)),
    list(student_t2_distrib(), list(mu = 0, sigma = 2, nu = 8)),
    list(gengamma2_distrib(), list(mean = 5, d = 3, p = 1.5))
  )
  for (cs in cases) {
    d <- cs[[1]]
    th <- cs[[2]]
    set.seed(1)
    res <- check_distrib(d, th, verbose = FALSE)
    expect_identical(nrow(res), 13L, label = d@distrib_name)
    expect_true(all(res$status == "OK"),
                info = paste(d@distrib_name,
                             paste(res$check[res$status != "OK"], collapse = ", ")))
    # the first parameter is the mean, which is the point of each of them
    expect_equal(mean(d, th), unname(th[[1]]), tolerance = 1e-8,
                 label = d@distrib_name)
  }

  # and the two that name a variance name it
  expect_equal(variance(lognormal2_distrib(), list(mean = 3, var = 2)), 2)
  expect_equal(variance(student_t2_distrib(), list(mu = 0, sigma = 2, nu = 8)), 4)
})


test_that("the four families agree with their twins where one exists", {
  # weibull3 against weibull1 at the matching scale, and student_t2 against
  # student_t1 at the matching scale: two implementations of one law.
  m <- 4
  s <- 1.7
  scale <- m / gamma(1 + 1 / s)
  set.seed(8)
  y <- rweibull(50, shape = s, scale = scale)
  expect_equal(distrib_pdf(weibull3_distrib(), y, list(mean = m, sigma = s)),
               distrib_pdf(weibull1_distrib(), y, list(mu = scale, sigma = s)),
               tolerance = 1e-12)

  nu <- 8
  sd <- 2
  sc <- sd * sqrt((nu - 2) / nu)
  set.seed(9)
  yt <- rt(50, df = nu) * sc
  expect_equal(
    distrib_pdf(student_t2_distrib(), yt, list(mu = 0, sigma = sd, nu = nu)),
    distrib_pdf(student_t1_distrib(), yt, list(mu = 0, sigma = sc, nu = nu)),
    tolerance = 1e-12
  )
})


test_that("a broken map would be caught", {
  # The paired injection: agreement above is evidence only if a wrong map fails.
  bad <- reparametrize(
    gaussian1_distrib(),
    map = function(psi) list(mu = psi$mu, sigma = sqrt(psi$sigma2) * 1.05),
    params = c("mu", "sigma2"),
    bounds = list(mu = c(-Inf, Inf), sigma2 = c(0, Inf)),
    links = list(mu = il(), sigma2 = ll())
  )
  set.seed(7)
  y <- rnorm(40, 1.2, sqrt(3.5))
  th <- list(mu = 1.2, sigma2 = 3.5)
  expect_false(isTRUE(all.equal(distrib_pdf(bad, y, th),
                                distrib_pdf(gaussian2_distrib(), y, th))))
})
