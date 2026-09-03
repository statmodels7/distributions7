test_that("expected_by_opg returns the outer product of the observed scores", {
  d <- pig1_distrib()
  th <- list(mu = c(2, 3, 4), sigma = c(0.5, 0.8, 1.1))
  y <- c(1, 4, 7)
  g <- distrib_gradient(d, y, th)
  out <- distributions7:::expected_by_opg(d, y, th)

  expect_identical(names(out), hess_names(d@params))
  expect_equal(out$mu_mu, -g$mu * g$mu)
  expect_equal(out$sigma_sigma, -g$sigma * g$sigma)
  expect_equal(out$mu_sigma, -g$mu * g$sigma)
  expect_true(all(vapply(out, length, integer(1)) == length(y)))
})

test_that("a scalar theta is recycled to the length of y", {
  d <- pig1_distrib()
  y <- c(0, 2, 5, 9)
  out <- distributions7:::expected_by_opg(d, y, list(mu = 3, sigma = 0.8))
  expect_true(all(vapply(out, length, integer(1)) == length(y)))
})

test_that("opg is no longer an alias of bartlett, and is the default", {
  # The two are readings of one identity and NOT the same computation: this is
  # what 0.44.0 changed, so the check is that they DIFFER per observation on a
  # family with no closed form, and that the default is the cheap one.
  d <- pig1_distrib()
  y <- c(1, 3, 6, 10)
  th <- list(mu = 4, sigma = 0.7)

  opg <- distrib_expected_hessian(d, y, th, approx = "opg")
  bar <- distrib_expected_hessian(d, y, th, approx = "bartlett")
  def <- distrib_expected_hessian(d, y, th)

  expect_equal(def, opg)
  expect_false(isTRUE(all.equal(opg$mu_mu, bar$mu_mu)))
  # bartlett does not depend on the observation, being an expectation; opg does
  expect_equal(length(unique(round(bar$mu_mu, 12))), 1L)
  expect_gt(length(unique(round(opg$mu_mu, 12))), 1L)
})

test_that("the outer product estimates what the identity computes", {
  # THE LICENCE FOR THE DEFAULT. The second Bartlett identity says the two
  # agree in expectation, so on a sample from the family the summed outer
  # product must approach the summed exact expectation. Checked on a family
  # whose expected information is closed form, so the reference is exact and
  # shares no arithmetic with the outer product.
  set.seed(11)
  d <- gaussian1_distrib()
  th <- list(mu = 1.5, sigma = 2)
  y <- distrib_rng(d, 20000L, th)

  opg <- distributions7:::expected_by_opg(d, y, th)
  ex <- distrib_expected_hessian(d, y, th)
  for (nm in names(ex)) {
    expect_equal(sum(opg[[nm]]) / length(y), sum(ex[[nm]]) / length(y),
                 tolerance = 0.05)
  }
})

test_that("a family with a closed form is untouched by approx", {
  for (fn in c("gaussian1_distrib", "gamma1_distrib", "poisson_distrib",
               "negbin2_distrib")) {
    d <- do.call(fn, list())
    th <- stats::setNames(as.list(rep(1, length(d@params))), d@params)
    y <- distrib_rng(d, 40L, th)
    a <- distrib_expected_hessian(d, y, th, approx = "opg")
    b <- distrib_expected_hessian(d, y, th, approx = "bartlett")
    expect_equal(a, b, info = fn)
  }
})

test_that("orders 3 and 4 route opg to bartlett", {
  # the outer product of scores is the SECOND-order identity and has no
  # counterpart above it
  d <- gaussian1_distrib()
  th <- list(mu = 0, sigma = 1)
  y <- distrib_rng(d, 30L, th)
  for (k in 3:4) {
    expect_equal(
      distributions7:::expected_derivative(d, y, th, order = k, approx = "opg"),
      distributions7:::expected_derivative(d, y, th, order = k,
                                           approx = "bartlett")
    )
  }
})

test_that("expected_hessian_exact is exported and names the six families", {
  expect_true(expected_hessian_exact(gaussian1_distrib()))
  expect_true(expected_hessian_exact(poisson_distrib()))
  for (fn in c("pig1_distrib", "pig2_distrib", "pseudohuber_distrib",
               "skewnormal1_distrib", "skewnormal2_distrib", "skewt_distrib")) {
    expect_false(expected_hessian_exact(do.call(fn, list())), info = fn)
  }
})

test_that("a fit reads its standard errors off the observed information where the expected one is not exact", {
  skip_on_cran()
  set.seed(3)
  d <- pig1_distrib()
  y <- distrib_rng(d, 300L, list(mu = 4, sigma = 0.6))
  fit <- fit_distrib(d, y)

  # the stored matrix IS the observed one, not the outer product
  V_obs <- vcov(fit, information = "observed")
  expect_equal(vcov(fit), V_obs)

  # and the expected route is reachable, different, and recomputed
  V_exp <- vcov(fit, information = "expected")
  expect_false(isTRUE(all.equal(V_obs, V_exp)))
  expect_identical(dimnames(V_exp), dimnames(V_obs))

  # a family WITH a closed form keeps the expected information
  set.seed(4)
  dg <- gaussian1_distrib()
  fg <- fit_distrib(dg, distrib_rng(dg, 300L, list(mu = 0, sigma = 1)))
  expect_equal(vcov(fg), vcov(fg, information = "expected"))
})

test_that("vcov's link scale is the parameter scale under the delta method", {
  skip_on_cran()
  set.seed(5)
  d <- gaussian1_distrib()
  fit <- fit_distrib(d, distrib_rng(d, 200L, list(mu = 1, sigma = 2)))
  for (info in c("observed", "expected")) {
    Ve <- vcov(fit, "link", information = info)
    Vt <- vcov(fit, "parameter", information = info)
    J <- diag(c(1, coef(fit)[["sigma"]]))
    expect_equal(unname(Vt), J %*% Ve %*% J, info = info)
  }
})
