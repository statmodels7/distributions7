# Two derivatives in the response, one in each parameter.

test_that("the closed forms agree with Richardson on the response Hessian", {
  skip_if_not_installed("numDeriv")
  y <- c(-2.1, -0.4, 0.3, 1.7, 3.2)
  cases <- list(
    list(d = gaussian1_distrib(), th = list(mu = 0.4, sigma = 1.3)),
    list(d = student_t1_distrib(),
         th = list(mu = -0.2, sigma = 0.9, nu = 6))
  )
  for (cs in cases) {
    d <- cs$d
    th <- cs$th
    got <- distrib_cross2_y(d, y, th)
    # the reference differentiates the ANALYTIC response Hessian in each
    # parameter by Richardson extrapolation, which shares no arithmetic with
    # the closed forms
    for (p in d@params) {
      ref <- vapply(seq_along(y), function(i) {
        f <- function(v) {
          t2 <- th
          t2[[p]] <- v
          distrib_hess_y(d, y[i], t2)
        }
        numDeriv::grad(f, th[[p]])
      }, numeric(1))
      expect_equal(got[[p]], ref, tolerance = 1e-6,
                   info = paste(d@distrib_name, p))
    }
  }
})

test_that("the fallback agrees with the closed form where both exist", {
  # the two routes share nothing: one is the written-out derivative, the other
  # a difference of the response Hessian
  y <- c(-1.3, 0.2, 2.4)
  d <- gaussian1_distrib()
  th <- list(mu = 0.1, sigma = 1.7)
  expect_equal(distrib_cross2_y(d, y, th), numerical_cross2_y(d, y, th),
               tolerance = 1e-7)

  d2 <- student_t1_distrib()
  th2 <- list(mu = 0.3, sigma = 1.1, nu = 8)
  expect_equal(distrib_cross2_y(d2, y, th2), numerical_cross2_y(d2, y, th2),
               tolerance = 1e-6)
})

test_that("a family with no closed form still answers", {
  y <- c(0.4, 1.2, 2.9)
  d <- gamma1_distrib()
  th <- list(mu = 1.5, phi = 0.8)
  got <- distrib_cross2_y(d, y, th)
  expect_named(got, d@params)
  expect_true(all(vapply(got, function(v) all(is.finite(v)), logical(1))))
  expect_true(all(lengths(got) == length(y)))
})

test_that("the link scale is the first-order diagonal chain rule", {
  # the response derivatives are untouched by a reparametrization of theta, so
  # only h'(eta) enters, exactly as for the gradient and for cross_y
  y <- c(-0.7, 0.5, 1.9)
  d <- gaussian1_distrib()
  th <- list(mu = 0.2, sigma = 1.4)
  par <- distrib_cross2_y(d, y, th, scale = "parameter")
  lnk <- distrib_cross2_y(d, y, th, scale = "link")
  # the log link on sigma: h'(eta) = sigma
  expect_equal(lnk$sigma, par$sigma * th$sigma, tolerance = 1e-12)
  expect_equal(lnk$mu, par$mu, tolerance = 1e-12)
})

test_that("the gaussian's curvature does not move with the location", {
  # -1/sigma^2 has no mu in it, and a component that came back non-zero would
  # say the location changed the shape
  y <- rnorm(20)
  got <- distrib_cross2_y(gaussian1_distrib(), y,
                          list(mu = 0.9, sigma = 2.2))
  expect_true(all(got$mu == 0))
  expect_equal(got$sigma, rep(2 / 2.2^3, 20), tolerance = 1e-14)
})

test_that("fixed() delegates it, subset to the free parameters", {
  # a penalty is built on fixed(gaussian1_distrib(), mu = 0), so a wrapper
  # that did not delegate would send every ridge through the fallback while
  # the closed form sat one class away
  d <- gaussian1_distrib()
  f <- fixed(d, mu = 0)
  y <- c(-1.2, 0.3, 2.1)
  got <- distrib_cross2_y(f, y, list(sigma = 1.6))
  expect_named(got, "sigma")
  expect_equal(got$sigma, distrib_cross2_y(d, y, list(mu = 0, sigma = 1.6))$sigma,
               tolerance = 1e-14)
  # and it is the closed form, not a difference of it
  expect_equal(got$sigma, rep(2 / 1.6^3, 3), tolerance = 1e-14)
})
