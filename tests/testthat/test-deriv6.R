# The sixth derivative of the log-density and the fifth derivative of the
# expected information.
#
# Both are one tensor stencil (numericals7::fd_tensor()) on the highest
# analytic order of the family, taken on the scale requested: in the
# parameters, or in the linear predictors on the exact link-scale order below.
# The references below are closed forms written in the test and 50-digit
# values from sympy/mpmath (misure_quinto_ordine/ref.py).

rel_all <- function(a, b) {
  max(abs(unlist(a[names(b)]) - unlist(b))) / max(abs(unlist(b)))
}

# d^a/dmu^a d^b/dsigma^b of the gaussian log-density, from
# -log(sigma) - r^2 / (2 sigma^2) with r = y - mu
gauss_d <- function(y, mu, s, a, b) {
  r <- y - mu
  dpow <- function(k, b) {     # d^b/ds^b s^(-k)
    prod(-(k + seq_len(b) - 1)) * s^(-k - b)
  }
  logpart <- if (a == 0 && b > 0) (-1)^b * factorial(b - 1) / s^b else 0
  quad <- switch(a + 1,
    -r^2 / 2 * dpow(2, b),
    r * dpow(2, b),
    -dpow(2, b) + 0 * r,
    0 * r, 0 * r, 0 * r, 0 * r)
  logpart + quad
}


test_that("the sixth order of the gaussian matches its closed form", {
  y <- c(-1.3, 0.4, 2.2)
  mu <- 0.3; s <- 1.2
  d <- gaussian1_distrib()
  got <- distrib_deriv6(d, y, list(mu = mu, sigma = s))
  expect_identical(names(got), deriv_names(c("mu", "sigma"), 6))
  ref <- lapply(deriv_indices(c("mu", "sigma"), 6), function(t) {
    a <- sum(t == 1L)
    gauss_d(y, mu, s, a, 6L - a)
  })
  names(ref) <- names(got)
  expect_lt(rel_all(got, ref), 1e-9)
})


test_that("the components are the order-6 enumeration", {
  y <- c(-1.3, 0.4, 2.2)
  s <- distrib_deriv6(student_t1_distrib(), y,
                      list(mu = 0.2, sigma = 1.1, nu = 6))
  expect_identical(names(s), deriv_names(c("mu", "sigma", "nu"), 6))
  expect_length(s, 28L)
  expect_true(all(lengths(s) == length(y)))
})


test_that("the sixth order matches 50-digit references", {
  # sympy derivatives of the closed-form log-density at 50 digits
  st <- distrib_deriv6(student_t1_distrib(), c(-1.3, 0.4),
                       list(mu = 0.2, sigma = 1.1, nu = 6))
  expect_equal(st[["mu_mu_mu_mu_mu_mu"]],
               c(0.972358780679270819713179415249,
                 -1.94941826785844193281670982907),
               tolerance = 1e-8)
  expect_equal(st[["sigma_sigma_sigma_nu_nu_nu"]],
               c(0.0415511054026035253326374232105,
                 -0.00119627621457076133001443673576),
               tolerance = 1e-6)
  z <- distrib_deriv6(zero_inflated(poisson_distrib()), c(0, 1, 3),
                      list(mu = 2, zi = 0.3))
  expect_equal(z[["mu_mu_mu_zi_zi_zi"]],
               c(-25.4327275834241987378824065617, 0, 0),
               tolerance = 1e-7)
})


test_that("a higher-accuracy rule agrees", {
  cases <- list(
    list(d = gamma1_distrib(), y = c(0.4, 1.6), th = list(mu = 2, phi = 0.7)),
    list(d = negbin2_distrib(), y = c(0, 5), th = list(mu = 3, theta = 2.2)),
    list(d = beta1_distrib(), y = c(0.15, 0.7), th = list(mu = 0.4, phi = 6)))
  for (cs in cases) {
    a4 <- distrib_deriv6(cs$d, cs$y, cs$th)
    a6 <- numerical_deriv6(cs$d, cs$y, cs$th, accuracy = 6L)
    expect_lt(rel_all(a4, a6), 1e-7, label = cs$d@distrib_name)
  }
})


test_that("the base order is the fifth where it is analytic", {
  expect_identical(distributions7:::deriv6_base(gaussian1_distrib()), 5L)
  expect_identical(
    distributions7:::deriv6_base(zero_inflated(poisson_distrib())), 5L)
  # a reparametrized family has no fifth order of its own
  r <- reparametrize(gaussian1_distrib(),
                     map = function(psi) list(mu = psi$a, sigma = exp(psi$b)),
                     params = c("a", "b"),
                     bounds = list(a = c(-Inf, Inf), b = c(-Inf, Inf)),
                     links = list(a = linkfunctions7::identity_link(),
                                  b = linkfunctions7::identity_link()))
  expect_identical(distributions7:::deriv6_base(r), 4L)
})


test_that("the link scale is a stencil in the predictors", {
  # The log link of the poisson gives l(eta) = y eta - exp(eta) - log y!, so
  # every derivative of order two or more in eta is -mu, also at a small mean,
  # where a composition of the parameter scale loses two digits in five
  d <- poisson_distrib()
  for (mu in c(3, 0.05)) {
    got <- distrib_deriv6(d, c(0, 2, 5), list(mu = mu), scale = "link")
    expect_lt(max(abs(got[[1]] + mu)) / mu, 1e-9)
  }
  d <- gamma1_distrib()
  y <- c(0.4, 1.6)
  th <- list(mu = 2, phi = 0.7)
  expect_identical(distrib_deriv6(d, y, th, scale = "link"),
                   numerical_deriv6(d, y, th, scale = "link"))
})


test_that("a multivariate family has a sixth order", {
  d <- mvgaussian1_distrib(2)
  th <- list(mu1 = 0.1, mu2 = 0.3, sigma_log_L1 = -0.2, sigma_log_L2 = 0.1,
             sigma_L2.1 = 0.4)
  y <- matrix(c(0.1, 0.5, -0.3, 1.2), ncol = 2)
  r <- distrib_deriv6(d, y, th)
  expect_length(r, choose(5 + 5, 6))
  expect_true(all(lengths(r) == 2L))
  # the log-density is quadratic in the mean
  expect_equal(r$mu1_mu1_mu1_mu2_mu2_mu2, c(0, 0), tolerance = 1e-6)
})


test_that("the fifth derivative of the information has its closed forms", {
  # poisson: E[l_mumu] = -1/mu, so its fifth derivative is 120 / mu^6; on the
  # log link the information is mu = exp(eta), all of whose derivatives are mu
  d <- poisson_distrib()
  y <- c(0, 2, 5)
  p <- distrib_d5expected_hessian(d, y, list(mu = 3))
  expect_identical(names(p), d5expected_names("mu"))
  expect_equal(p[[1]], rep(120 / 3^6, 3), tolerance = 1e-9)
  for (mu in c(3, 0.05)) {
    l <- distrib_d5expected_hessian(d, y, list(mu = mu), scale = "link")
    expect_lt(max(abs(l[[1]] + mu)) / mu, 1e-9)
  }
  # gaussian: E[l_mumu] = -1/sigma^2 is -exp(-2 eta) on the log link, and
  # E[l_sigmasigma] = -2/sigma^2 is the constant -2 there
  g <- gaussian1_distrib()
  s <- 1.3
  pg <- distrib_d5expected_hessian(g, 0, list(mu = 0, sigma = s))
  expect_equal(pg[["sigma_sigma_sigma_sigma_sigma_sigma_sigma"]], 1440 / s^7,
               tolerance = 1e-9)
  expect_equal(pg[["mu_mu_sigma_sigma_sigma_sigma_sigma"]], 720 / s^7,
               tolerance = 1e-9)
  lg <- distrib_d5expected_hessian(g, 0, list(mu = 0, sigma = s),
                                   scale = "link")
  expect_equal(lg[["mu_mu_sigma_sigma_sigma_sigma_sigma"]], 32 / s^2,
               tolerance = 1e-9)
  expect_lt(abs(lg[["sigma_sigma_sigma_sigma_sigma_sigma_sigma"]]), 1e-9)
})


test_that("the fifth derivative of the information matches 50-digit references", {
  st <- distrib_d5expected_hessian(student_t1_distrib(), 0,
                                   list(mu = 0.2, sigma = 1.1, nu = 6))
  expect_equal(st[["nu_nu_nu_nu_nu_nu_nu"]], 0.00077745264914130479736710156973,
               tolerance = 1e-6)
  expect_equal(st[["sigma_nu_sigma_sigma_nu_nu_nu"]],
               -0.00238086254774621436237925367636, tolerance = 1e-6)
  z <- distrib_d5expected_hessian(zero_inflated(poisson_distrib()), 0,
                                  list(mu = 2, zi = 0.3))
  expect_equal(z[["mu_zi_mu_mu_zi_zi_zi"]], 23.5270520667887770801107631946,
               tolerance = 1e-6)
})


test_that("the stencil differences the base order that the family carries", {
  expect_identical(distributions7:::d5expected_base(gamma1_distrib()), 4L)
  expect_identical(
    distributions7:::d5expected_base(truncated(gaussian1_distrib(), lower = 0)),
    4L)
  expect_true(is.na(distributions7:::d5expected_base(mvgaussian1_distrib(2))))
  expect_error(
    distrib_d5expected_hessian(mvgaussian1_distrib(2),
                               matrix(c(0.1, 0.5), ncol = 2),
                               list(mu1 = 0, mu2 = 0, sigma_log_L1 = 0,
                                    sigma_log_L2 = 0, sigma_L2.1 = 0)),
    "no analytic second derivative")
})


test_that("a method object gives the derivative of its own matrix", {
  d <- poisson_distrib()
  y <- c(0, 2, 5)
  mu <- 3
  # expect_opg(): -(y/mu - 1)^2 at the data, differentiated five times
  o <- distrib_d5expected_hessian(d, y, list(mu = mu), method = expect_opg())
  expect_equal(o[[1]], 720 * y^2 / mu^7 - 240 * y / mu^6, tolerance = 1e-7)
  # on the log link it is -(y - exp(eta))^2, read with link-scale derivatives
  ol <- distrib_d5expected_hessian(d, y, list(mu = mu), scale = "link",
                                   method = expect_opg())
  expect_equal(ol[[1]], 2 * y * mu - 32 * mu^2, tolerance = 1e-8)
  # the score identity integrates to the closed form
  s <- distrib_d5expected_hessian(d, y, list(mu = mu),
                                  method = expect_integrate(identity = "score"))
  expect_equal(s[[1]], rep(120 / mu^6, 3), tolerance = 1e-6)
  # the Hessian identity would need the seventh order
  expect_error(distrib_d5expected_hessian(d, y, list(mu = mu),
                                          method = expect_integrate()),
               'identity = "score"')
})


test_that("the arguments are checked as at the order below", {
  d <- gaussian1_distrib()
  expect_error(distrib_deriv6(d, 0, list(mu = 0, sigma = -1)))
  expect_error(distrib_d5expected_hessian(d, 0, list(mu = 0, sigma = 1),
                                          method = "opg"), "ExpectMethod")
  expect_error(distrib_deriv6(d, 0, list(mu = 0, sigma = 1), scale = "eta"))
})


# A family that registers its own sixth order, or its own fifth derivative of
# the information, returns the parameter scale, and the generic carries it to
# the link scale by the Faa di Bruno composition with the inverse links
# up to the sixth. The families below exist in this file only: a subclass of a
# shipped family with closed forms written out, and a counter that records
# whether, and on which scale, the method was called.
d6_family_ <- function(base, deriv6 = NULL, d5e = NULL, log = new.env()) {
  cls <- S7::new_class(paste0("D6Test_", base@distrib_name), package = NULL,
                       parent = S7::S7_class(base))
  log$calls6 <- 0L
  log$calls5e <- 0L
  if (!is.null(deriv6)) {
    S7::method(distrib_deriv6, cls) <- function(distrib, y, theta,
        scale = c("parameter", "link"), ...) {
      log$calls6 <- log$calls6 + 1L
      log$scale6 <- scale
      deriv6(y, theta)
    }
  }
  if (!is.null(d5e)) {
    S7::method(distrib_d5expected_hessian, cls) <- function(distrib, y, theta,
        scale = c("parameter", "link"), method = NULL, ...) {
      log$calls5e <- log$calls5e + 1L
      log$scale5e <- scale
      d5e(y, theta)
    }
  }
  do.call(cls, S7::props(base))
}

d6_family <- function(...) suppressMessages(d6_family_(...))

d6_poisson <- function(log) {
  d6_family(
    poisson_distrib(),
    deriv6 = function(y, theta) {
      stats::setNames(list(-120 * y / theta[[1]]^6), deriv_names("mu", 6L))
    },
    d5e = function(y, theta) {
      stats::setNames(list(rep(120 / theta[[1]]^6, length(y))),
                      d5expected_names("mu"))
    },
    log = log)
}

d6_gaussian <- function(log) {
  d6_family(
    gaussian1_distrib(),
    deriv6 = function(y, theta) {
      nm <- deriv_names(c("mu", "sigma"), 6L)
      tok <- strsplit(nm, "_")
      stats::setNames(lapply(tok, function(t) {
        gauss_d(y, theta[[1]], theta[[2]], sum(t == "mu"), sum(t == "sigma"))
      }), nm)
    },
    d5e = function(y, theta) {
      nm <- d5expected_names(c("mu", "sigma"))
      s <- theta[[2]]
      stats::setNames(lapply(strsplit(nm, "_"), function(t) {
        if (!all(t[3:7] == "sigma")) return(rep(0, length(y)))
        pair <- paste(t[1:2], collapse = "_")
        rep(switch(pair, mu_mu = 720, sigma_sigma = 1440, 0) / s^7, length(y))
      }), nm)
    },
    log = log)
}


test_that("an own sixth order is composed with the inverse links", {
  lg <- new.env()
  d <- d6_poisson(lg)
  expect_true(distributions7:::has_exact_deriv6(d))
  expect_false(distributions7:::has_exact_deriv6(poisson_distrib()))
  y <- c(0, 2, 5)
  for (mu in c(3, 0.05)) {
    th <- list(mu = mu)
    lg$calls6 <- 0L
    got <- distrib_deriv6(d, y, th, scale = "link")
    # the family's method ran once, and was asked for the parameter scale
    expect_identical(lg$calls6, 1L)
    expect_identical(lg$scale6, "parameter")
    # on the log link every derivative of order two or more is -mu
    expect_lt(max(abs(got[[1]] + mu)) / mu, 1e-11)
    # and the parameter scale is the method's own result
    expect_equal(distrib_deriv6(d, y, th)[[1]], -120 * y / mu^6)
  }
})


test_that("without an own sixth order the stencil in the predictors is taken", {
  d <- poisson_distrib()
  y <- c(0, 2, 5)
  for (mu in c(3, 0.05)) {
    th <- list(mu = mu)
    got <- distrib_deriv6(d, y, th, scale = "link")
    expect_identical(got, numerical_deriv6(d, y, th, scale = "link"))
    expect_lt(max(abs(got[[1]] + mu)) / mu, 1e-9)
  }
})


test_that("the composition at order six agrees with the stencil in the predictors", {
  lg <- new.env()
  d <- d6_gaussian(lg)
  b <- gaussian1_distrib()
  y <- c(-0.7, 0.4, 1.9)
  th <- list(mu = 0.3, sigma = 1.4)
  comp <- distrib_deriv6(d, y, th, scale = "link")
  expect_identical(lg$calls6, 1L)
  sten <- distrib_deriv6(b, y, th, scale = "link")
  expect_identical(names(comp), names(sten))
  expect_lt(rel_all(comp, sten), 1e-7)
  # the parameter scale is the closed form, untouched by the composition
  par <- distrib_deriv6(d, y, th)
  expect_equal(par[["mu_mu_sigma_sigma_sigma_sigma"]],
               gauss_d(y, 0.3, 1.4, 2, 4))
})


test_that("an own fifth derivative of the information is composed with the inverse links", {
  lg <- new.env()
  d <- d6_poisson(lg)
  y <- c(0, 2, 5)
  for (mu in c(3, 0.05)) {
    th <- list(mu = mu)
    lg$calls5e <- 0L
    got <- distrib_d5expected_hessian(d, y, th, scale = "link")
    expect_identical(lg$calls5e, 1L)
    expect_identical(lg$scale5e, "parameter")
    expect_identical(names(got), d5expected_names("mu"))
    # the information on the log link is -mu at every order of eta
    expect_lt(max(abs(got[[1]] + mu)) / mu, 1e-11)
    expect_equal(distrib_d5expected_hessian(d, y, th)[[1]],
                 rep(120 / mu^6, 3))
    # the stencil of the family without the method, for comparison
    expect_lt(max(abs(distrib_d5expected_hessian(poisson_distrib(), y, th,
                                                 scale = "link")[[1]] + mu)) /
                mu, 1e-9)
  }
  expect_true(distributions7:::has_exact_d5expected(d))
  expect_false(distributions7:::has_exact_d5expected(poisson_distrib()))
})


test_that("the composition of the information at order five has the exact values", {
  lg <- new.env()
  g <- d6_gaussian(lg)
  s <- 1.3
  th <- list(mu = 0, sigma = s)
  comp <- distrib_d5expected_hessian(g, 0, th, scale = "link")
  expect_identical(lg$calls5e, 1L)
  # the values of the stencil-based test above, here to rounding
  expect_equal(comp[["mu_mu_sigma_sigma_sigma_sigma_sigma"]], 32 / s^2,
               tolerance = 1e-12)
  expect_lt(abs(comp[["sigma_sigma_sigma_sigma_sigma_sigma_sigma"]]), 1e-11)
  sten <- distrib_d5expected_hessian(gaussian1_distrib(), 0, th,
                                     scale = "link")
  expect_identical(names(comp), names(sten))
  expect_lt(rel_all(comp, sten), 1e-7)
})
