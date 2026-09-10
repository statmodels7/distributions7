# The centered parametrization reaches the direct one through a map that
# degenerates at the ceiling of the skewness. These pin the two things that
# repair cost: the map stays finite everywhere the link can reach, and a
# failure of it is reported in the centered family's own terms.

b_const <- function() sqrt(2 / pi)

# The expressions the repair replaced, transcribed here so the negative
# control does not read the package's own arithmetic back to itself.
old_map <- function(mu, sigma, g, s) {
  b <- b_const()
  cc <- s * (2 * (s * g) / (4 - pi))^(1 / 3)
  muz <- cc / sqrt(1 + cc^2)
  del <- muz / b
  om <- sigma / sqrt(1 - muz^2)
  c(mu = mu - om * muz, sigma = om, alpha = del / sqrt(1 - del^2))
}
old_d2 <- function(g, s) {
  b <- b_const()
  r <- s * (2 * (s * g) / (4 - pi))^(1 / 3)
  (b^2 + (b^2 - 1) * r^2) / b^2
}

test_that("the map is finite at every skewness the link can produce", {
  gm <- sn_max_skew()
  lk <- linkfunctions7::bounded_link(lwr = -gm, upr = gm)
  d <- skewnormal2_distrib()

  for (eta in c(20, 30, 36, 100, 800)) {
    for (sgn in c(1, -1)) {
      g <- linkfunctions7::linkinv(lk, sgn * eta)
      th <- list(mu = 0, sigma = 1, gamma1 = g)

      expect_true(all(is.finite(unlist(sn2_theta(th)))),
                  info = paste("eta =", sgn * eta))
      expect_true(all(is.finite(unlist(distrib_pdf(d, 0.5, th)))))
      expect_true(all(is.finite(unlist(distrib_gradient(d, 0.5, th)))))
      expect_true(all(is.finite(unlist(distrib_hessian(d, 0.5, th)))))
      expect_true(all(is.finite(unlist(distrib_deriv3(d, 0.5, th)))))
      expect_true(all(is.finite(unlist(distrib_deriv4(d, 0.5, th)))))
      expect_true(all(is.finite(unlist(
        distrib_gradient(d, 0.5, th, scale = "link")))))
    }
  }
})

test_that("the expressions the repair replaced fail exactly there", {
  # The control on the control: at the last representable skewness inside the
  # bound the old shape is not finite and the old 1 - delta^2 is negative,
  # while the quantity itself is of order 1e8 and perfectly representable.
  gm <- sn_max_skew()
  g <- gm * (1 - 2.22e-16)

  expect_false(is.finite(old_map(0, 1, g, 1)[["alpha"]]))
  expect_lt(old_d2(g, 1), 0)
  expect_true(is.finite(sn_one_minus_delta2(g, 1)))
  expect_gt(sn_one_minus_delta2(g, 1), 0)
  expect_true(is.finite(sn2_theta(list(mu = 0, sigma = 1, gamma1 = g))$alpha))
})

test_that("the map is unchanged where the old expressions were trustworthy", {
  worst <- 0
  for (g in c(-0.99, -0.9, -0.5, -1e-3, 0, 1e-3, 0.5, 0.9, 0.99)) {
    for (sg in c(0.1, 1, 100)) {
      s <- if (g >= 0) 1 else -1
      o <- old_map(2, sg, g, s)
      n <- unlist(sn_cp_to_dp(2, sg, g, s))
      worst <- max(worst, max(abs(o - n) / pmax(abs(o), 1e-12)))
    }
  }
  expect_lt(worst, 1e-12)
})

test_that("the shape agrees with the skewness it is built from", {
  # An independent route: skewnormal1's own skewness() shares no arithmetic
  # with the map, so the round trip checks the map rather than restating it.
  gm <- sn_max_skew()
  d1 <- skewnormal1_distrib()
  for (g in c(-0.99, -0.5, 1e-6, 0.5, 0.9, gm * (1 - 1e-9),
              gm * (1 - 2.22e-16))) {
    dp <- sn2_theta(list(mu = 0, sigma = 1, gamma1 = g))
    expect_equal(skewness(d1, dp), g, tolerance = 1e-12)
  }
})

test_that("a skewness the map cannot carry is reported as skewnormal2's own", {
  th <- list(mu = 0, sigma = 1, gamma1 = sn_max_skew())
  expect_error(sn2_theta(th), "skew normal2", fixed = TRUE)
  expect_error(sn2_theta(th), "gamma1", fixed = TRUE)
  # and never in the parent's vocabulary, which is what the repair is for
  err <- tryCatch(sn2_theta(th), error = conditionMessage)
  expect_false(grepl("skew normal1", err, fixed = TRUE))
  expect_false(grepl("alpha", err, fixed = TRUE))
})

test_that("a scale the map cannot carry is reported as skewnormal2's own", {
  th <- list(mu = 0, sigma = 1.1e308, gamma1 = 0.99)
  err <- tryCatch(sn2_theta(th), error = conditionMessage)
  expect_true(grepl("skew normal2", err, fixed = TRUE))
  expect_true(grepl("sigma", err, fixed = TRUE))
  expect_false(grepl("skew normal1", err, fixed = TRUE))
})

test_that("an ordinary triple passes the guard untouched", {
  th <- list(mu = 1, sigma = 2, gamma1 = 0.4)
  dp <- sn_cp_to_dp(1, 2, 0.4, 1)
  expect_null(sn2_reject_unmappable(dp, th))
  expect_identical(sn2_theta(th), dp)
})
