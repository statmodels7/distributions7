# Three derivatives in the response, one in each parameter.

test_that("the location identity agrees with Richardson on deriv3_y", {
  skip_if_not_installed("numDeriv")
  y <- c(-2.1, -0.4, 0.3, 1.7, 3.2)
  cases <- list(
    list(d = student_t1_distrib(),
         th = list(mu = -0.2, sigma = 0.9, nu = 6)),
    list(d = logistic_distrib(), th = list(mu = 0.3, sigma = 1.4)),
    list(d = gaussian1_distrib(), th = list(mu = 0.4, sigma = 1.3))
  )
  for (cs in cases) {
    d <- cs$d
    th <- cs$th
    got <- distrib_cross3_y(d, y, th)
    expect_named(got, d@params)
    # the reference differentiates the ANALYTIC third response derivative by
    # Richardson extrapolation, which shares no arithmetic with deriv4
    for (p in d@params) {
      ref <- vapply(seq_along(y), function(i) {
        f <- function(v) {
          t2 <- th
          t2[[p]] <- v
          distrib_deriv3_y(d, y[i], t2)
        }
        numDeriv::grad(f, th[[p]])
      }, numeric(1))
      expect_equal(got[[p]], ref, tolerance = 1e-6,
                   info = paste(d@distrib_name, p))
    }
  }
})

test_that("a location family takes the identity and not the difference", {
  # a tolerance cannot tell an accurate difference from the identity, so the
  # route is asserted by IDENTITY against the fourth derivative it reads
  d <- student_t1_distrib()
  y <- c(-1.3, 0.2, 2.4)
  th <- list(mu = 0.1, sigma = 1.7, nu = 5)
  got <- distrib_cross3_y(d, y, th)
  d4 <- distrib_deriv4(d, y, th)
  expect_identical(got$mu, -d4[["mu_mu_mu_mu"]])
  expect_identical(got$sigma, -d4[["mu_mu_mu_sigma"]])
  expect_identical(got$nu, -d4[["mu_mu_mu_nu"]])
  # and the fallback, which is the difference, agrees with it
  expect_equal(got, numerical_cross3_y(d, y, th), tolerance = 1e-6)
})

test_that("a family with no identity takes the difference and answers", {
  y <- c(0.4, 1.2, 2.9)
  d <- gamma1_distrib()
  th <- list(mu = 1.5, phi = 0.8)
  got <- distrib_cross3_y(d, y, th)
  expect_named(got, d@params)
  expect_true(all(vapply(got, function(v) all(is.finite(v)), logical(1))))
  expect_true(all(lengths(got) == length(y)))
  expect_identical(got, numerical_cross3_y(d, y, th))
})

test_that("the link scale is the first-order diagonal chain rule", {
  y <- c(-0.7, 0.5, 1.9)
  d <- student_t1_distrib()
  th <- list(mu = 0.2, sigma = 1.4, nu = 7)
  par <- distrib_cross3_y(d, y, th, scale = "parameter")
  lnk <- distrib_cross3_y(d, y, th, scale = "link")
  # log links on sigma and nu: h'(eta) is the parameter itself
  expect_equal(lnk$sigma, par$sigma * th$sigma, tolerance = 1e-12)
  expect_equal(lnk$nu, par$nu * th$nu, tolerance = 1e-12)
  expect_equal(lnk$mu, par$mu, tolerance = 1e-12)
})

test_that("fixed() delegates it, subset to the free parameters", {
  # the heavy-tailed prior of penalties7 is fixed(student_t1_distrib(), mu = 0),
  # so a wrapper that did not delegate would difference a closed form
  d <- student_t1_distrib()
  f <- fixed(d, mu = 0)
  y <- c(-1.2, 0.3, 2.1)
  got <- distrib_cross3_y(f, y, list(sigma = 1.6, nu = 4))
  expect_named(got, c("sigma", "nu"))
  full <- distrib_cross3_y(d, y, list(mu = 0, sigma = 1.6, nu = 4))
  expect_identical(got$sigma, full$sigma)
  expect_identical(got$nu, full$nu)
})
