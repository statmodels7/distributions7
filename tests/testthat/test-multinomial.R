# The multinomial: the first multivariate family here that is DISCRETE, so
# its support is a finite set of points and every expectation is an exact sum
# rather than a quadrature or a sample.

test_that("the support is enumerated and carries the whole mass", {
  d <- multinomial_distrib(3, size = 5)
  th <- as.list(stats::setNames(c(0.3, -0.2), d@params))
  supp <- mv_support(d, th)

  expect_identical(nrow(supp), as.integer(choose(5 + 3 - 1, 3 - 1)))
  expect_true(all(rowSums(supp) == 5))
  expect_true(all(supp >= 0))
  expect_equal(sum(distrib_pdf(d, supp, th)), 1, tolerance = 1e-12)

  # against dmultinom, written out one point at a time
  pr <- mv_location(d, th) / 5
  expect_equal(distrib_pdf(d, supp, th),
               apply(supp, 1L, function(r) stats::dmultinom(r, prob = pr)))

  # a point off the support has no mass
  expect_identical(distrib_pdf(d, matrix(c(2, 2, 2), 1), th), 0)
  expect_identical(distrib_pdf(d, matrix(c(2.5, 1.5, 1), 1), th), 0)
})


test_that("compositions enumerates every way to split an integer", {
  expect_identical(nrow(compositions(3, 2)), 4L)
  expect_identical(compositions(3, 2),
                   matrix(c(0L, 1L, 2L, 3L, 3L, 2L, 1L, 0L), ncol = 2))
  for (n in 0:5) for (k in 1:4) {
    m <- compositions(n, k)
    expect_identical(nrow(m), as.integer(choose(n + k - 1, k - 1)))
    expect_true(all(rowSums(m) == n))
    expect_identical(anyDuplicated(m), 0L)
  }
})


test_that("the base class refuses a support it cannot enumerate", {
  expect_error(mv_support(mvgaussian_distrib(2), list()),
               "does not enumerate a support")
  expect_error(mv_support(mvstudent_t_distrib(2), list()),
               "does not enumerate a support")
})


test_that("the analytical derivatives match one Richardson differentiation", {
  skip_if_not_installed("numDeriv")
  d <- multinomial_distrib(3, size = 5)
  th <- as.list(stats::setNames(c(0.3, -0.2), d@params))
  q0 <- unlist(th)
  set.seed(1)
  y <- distrib_rng(d, 40, th)
  llv <- function(q) {
    sum(distrib_pdf(d, y, as.list(stats::setNames(q, d@params)), log = TRUE))
  }

  g <- distrib_gradient(d, y, th)
  expect_equal(vapply(g, sum, numeric(1)), numDeriv::grad(llv, q0),
               tolerance = 1e-6, ignore_attr = TRUE)

  H <- numDeriv::hessian(llv, q0)
  h <- distrib_hessian(d, y, th)
  nm <- hess_names(d@params)
  pr <- hess_pairs(d@params)
  pos <- stats::setNames(seq_along(d@params), d@params)
  for (m in seq_along(nm)) {
    expect_equal(sum(h[[m]]), H[pos[[pr[[m]][1]]], pos[[pr[[m]][2]]]],
                 tolerance = 1e-5, label = nm[m])
  }
})


test_that("the expected information is the exact sum over the support", {
  # This is what a finite support buys, and it is the reason mv_support()
  # exists: the check is an equality rather than a comparison against noise.
  d <- multinomial_distrib(3, size = 5)
  for (th in list(as.list(stats::setNames(c(0.3, -0.2), d@params)),
                  as.list(stats::setNames(c(-1.1, 0.7), d@params)))) {
    supp <- mv_support(d, th)
    mass <- distrib_pdf(d, supp, th)
    hs <- distrib_hessian(d, supp, th)
    eh <- distrib_expected_hessian(d, supp[1, , drop = FALSE], th)
    for (nm in names(eh)) {
      expect_equal(eh[[nm]][1], sum(mass * hs[[nm]]), tolerance = 1e-12,
                   label = nm)
    }
  }
})


test_that("the moments and the marginals are the ones written down", {
  d <- multinomial_distrib(3, size = 5)
  th <- as.list(stats::setNames(c(0.3, -0.2), d@params))
  pr <- mv_location(d, th) / 5

  expect_equal(sum(pr), 1)
  expect_equal(mv_sigma(d, th), 5 * (diag(pr) - tcrossprod(pr)),
               tolerance = 1e-12)
  expect_lt(abs(det(mv_sigma(d, th))), 1e-12)

  for (j in 1:3) {
    m <- mv_marginal(d, th, which = j)
    expect_true(S7::S7_inherits(m$distrib, BinomialDistrib))
    expect_equal(distrib_pdf(m$distrib, 0:5, m$theta),
                 stats::dbinom(0:5, 5, pr[j]), label = as.character(j))
  }
  expect_error(mv_marginal(d, th, which = 1:2), "one coordinate at a time")
})


test_that("a fit recovers the probabilities", {
  d <- multinomial_distrib(3, size = 5)
  th <- as.list(stats::setNames(c(0.3, -0.2), d@params))
  set.seed(5)
  y <- distrib_rng(d, 3000, th)
  f <- fit_distrib(d, y)
  expect_true(f@converged)
  expect_equal(unname(coef(f)), unname(unlist(th)), tolerance = 0.1)
})


test_that("the validator passes, and its normalisation is exact", {
  d <- multinomial_distrib(3, size = 5)
  th <- as.list(stats::setNames(c(0.3, -0.2), d@params))
  set.seed(1)
  res <- check_distrib(d, th, verbose = FALSE)
  expect_true(all(res$status == "OK"),
              info = paste(res$check[res$status != "OK"], collapse = ", "))
  # the support being finite, the normalisation is a sum and not a sample
  expect_lt(res$statistic[res$check == "density integrates to 1"], 1e-14)
  # and a discrete family has no derivative in its response
  expect_false("response derivatives vs finite differences" %in% res$check)
  expect_true(has_mv_support(d))
  expect_false(has_mv_grad_y(d))

  # a normalisation wrong by a thousandth would sit inside the Monte Carlo
  # error of an importance-sampling check and is caught by the exact one
  good <- S7::method(distrib_pdf, MultinomialDistrib)
  on.exit(S7::method(distrib_pdf, MultinomialDistrib) <- good, add = TRUE)
  suppressMessages(
    S7::method(distrib_pdf, MultinomialDistrib) <- function(distrib, y, theta, log = FALSE) {
      out <- good(distrib, y, theta, log = TRUE) + log(1.001)
      if (log) out else exp(out)
    }
  )
  set.seed(1)
  bad <- check_distrib(d, th, verbose = FALSE)
  expect_equal(bad$status[bad$check == "density integrates to 1"], "FAIL")
})


test_that("the constructor validates its dimension and its size", {
  expect_error(multinomial_distrib(1, size = 3), "at least 2")
  expect_error(multinomial_distrib(3, size = 0), "positive integer")
  expect_error(multinomial_distrib(3, size = 2.5), "positive integer")
  expect_error(multinomial_distrib(3, size = 4, probs = "nonsense"),
               "parameters7 parameter")
})
