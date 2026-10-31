# The derivatives of the expected information under an explicit method.
#
# expect_opg(): the exact derivative of the outer-product matrix at the data.
# expect_integrate() / expect_mc(): the derivative of the expectation, which
# moves the density as well as the integrand.

dens_only_dx <- function() {
  Dens <- S7::new_class("DensDexpected", parent = continuous_distrib)
  S7::method(distrib_pdf, Dens) <- function(distrib, y, theta, log = FALSE, ...) {
    stats::dnorm(y, theta[[1]], theta[[2]], log = log)
  }
  Dens(distrib_name = "dens", dimension = "univariate",
       bounds = c(-Inf, Inf), params = c("mu", "sigma"),
       params_interpretation = c(mu = "location", sigma = "scale"),
       n_params = 2L,
       params_bounds = list(mu = c(-Inf, Inf), sigma = c(0, Inf)),
       link_params = list(mu = linkfunctions7::identity_link(),
                          sigma = linkfunctions7::log_link()),
       params_smooth = c(mu = TRUE, sigma = TRUE))
}

central <- function(f, th, p, h = 1e-5) {
  tu <- th; td <- th
  tu[[p]] <- th[[p]] + h
  td[[p]] <- th[[p]] - h
  a <- f(tu); b <- f(td)
  Map(function(u, v) (u - v) / (2 * h), a, b)
}

test_that("expect_opg() gives the exact derivative of the outer-product matrix", {
  d <- gamma1_distrib()
  th <- list(mu = 2.5, phi = 0.7)
  y <- c(1.2, 3.4, 0.8, 2.2)
  E <- function(t) distrib_expected_hessian(d, y, t, method = expect_opg())
  d1 <- distrib_dexpected_hessian(d, y, th, method = expect_opg())
  expect_identical(names(d1), dexpected_names(d@params))
  for (p in d@params) {
    num <- central(E, th, p)
    for (nm in names(num)) {
      expect_equal(d1[[paste0(nm, "_", p)]], num[[nm]], tolerance = 1e-6,
                   label = paste0(nm, "_", p))
    }
  }
  # the second derivative is the derivative of the first
  d2 <- distrib_d2expected_hessian(d, y, th, method = expect_opg())
  expect_identical(names(d2), d2expected_names(d@params))
  P <- d@params
  D1 <- function(t) distrib_dexpected_hessian(d, y, t, method = expect_opg())
  for (e in seq_along(P)) {
    num <- central(D1, th, P[e])
    for (nm in hess_names(P)) {
      ij <- strsplit(nm, "_")[[1L]]
      ia <- match(ij[1L], P); ib <- match(ij[2L], P)
      for (c in seq_along(P)) {
        key2 <- d2expected_key(P, ia, ib, c, e)
        expect_equal(d2[[key2]], num[[dexpected_key(P, ia, ib, c)]],
                     tolerance = 1e-5, label = key2)
      }
    }
  }
})

test_that("the opg derivative is not the derivative of the expectation", {
  # the density moves in the expectation, and the data do not move in the
  # outer-product matrix, so the two derivatives differ
  d <- gamma1_distrib()
  th <- list(mu = 2.5, phi = 0.7)
  y <- 1.5
  exact <- distrib_dexpected_hessian(d, y, th)
  o <- distrib_dexpected_hessian(d, y, th, method = expect_opg())
  expect_false(isTRUE(all.equal(unlist(o), unlist(exact), tolerance = 1e-3)))
})

test_that("integration reproduces the closed form of the derivatives", {
  d <- gamma1_distrib()
  th <- list(mu = 2.5, phi = 0.7)
  y <- c(1.2, 3.4)
  for (id in c("hessian", "score")) {
    m <- expect_integrate(identity = id)
    expect_equal(distrib_dexpected_hessian(d, y, th, method = m),
                 distrib_dexpected_hessian(d, y, th), tolerance = 1e-6)
    expect_equal(distrib_d2expected_hessian(d, y, th, method = m),
                 distrib_d2expected_hessian(d, y, th), tolerance = 1e-6)
    expect_equal(distrib_dexpected_hessian(d, y, th, method = m, scale = "link"),
                 distrib_dexpected_hessian(d, y, th, scale = "link"),
                 tolerance = 1e-6)
  }
  m <- expect_integrate(identity = "score")
  expect_equal(distrib_d3expected_hessian(d, y, th, method = m),
               distrib_d3expected_hessian(d, y, th), tolerance = 1e-5)
  expect_equal(distrib_d4expected_hessian(d, y, th, method = m),
               distrib_d4expected_hessian(d, y, th), tolerance = 1e-4)
  # the Hessian identity reads observed derivatives of order 6 at order 4,
  # which distrib_deriv6() supplies, and of order 7 at order 5
  expect_equal(distrib_d4expected_hessian(d, y, th, method = expect_integrate()),
               distrib_d4expected_hessian(d, y, th), tolerance = 1e-6)
  expect_error(distrib_d5expected_hessian(d, y, th, method = expect_integrate()),
               'identity = "score"')
  expect_no_error(distrib_d3expected_hessian(d, y, th, method = expect_integrate()))
})

test_that("Monte Carlo estimates the derivative of the expectation", {
  d <- gamma1_distrib()
  th <- list(mu = 2.5, phi = 0.7)
  set.seed(1)
  ex <- distrib_dexpected_hessian(d, 1, th)
  for (id in c("hessian", "score")) {
    mc <- distrib_dexpected_hessian(d, 1, th,
                                    method = expect_mc(nsim = 100000, identity = id))
    expect_equal(unlist(mc), unlist(ex[names(mc)]), tolerance = 0.05)
  }
})

test_that("the default is the closed form, and a family without one has no default", {
  d <- dens_only_dx()
  th <- list(mu = 0, sigma = 2)
  expect_error(distrib_dexpected_hessian(d, 0, th), "approximates its expected information")
  g <- gaussian1_distrib()
  expect_identical(distrib_dexpected_hessian(g, 0, th),
                   distrib_dexpected_hessian(g, 0, th, method = NULL))
  # a method gives a density-only family its derivative
  o <- distrib_dexpected_hessian(d, c(-1, 0.5), th, method = expect_opg())
  expect_identical(names(o), dexpected_names(d@params))
  i <- distrib_dexpected_hessian(d, 0, th, method = expect_integrate())
  expect_equal(i, distrib_dexpected_hessian(g, 0, th), tolerance = 1e-6)
  expect_error(distrib_dexpected_hessian(g, 0, th, method = "opg"),
               "must be NULL or an ExpectMethod object")
})
