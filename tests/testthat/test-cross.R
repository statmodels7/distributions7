# distrib_cross_y(): the mixed second derivatives of the log-density, one in
# the response and one in each parameter. The independent reference throughout
# is the four-point mixed stencil applied directly to the log-density, which
# shares no code with either the closed forms or the fallback (the fallback
# differentiates distrib_grad_y, not the log-density).

mixed_stencil <- function(distrib, y, theta, param, hy = 1e-5, ht = 1e-5) {
  up <- theta
  dn <- theta
  up[[param]] <- theta[[param]] + ht
  dn[[param]] <- theta[[param]] - ht
  lp <- function(yy, th) distrib_pdf(distrib, yy, th, log = TRUE)
  (lp(y + hy, up) - lp(y - hy, up) - lp(y + hy, dn) + lp(y - hy, dn)) /
    (4 * hy * ht)
}

test_that("the gaussian closed form matches the mixed stencil", {
  d <- gaussian1_distrib()
  th <- list(mu = 0.3, sigma = 1.7)
  y <- c(-2, -0.5, 0.3, 1, 4)

  cr <- distrib_cross_y(d, y, th)
  expect_named(cr, c("mu", "sigma"))

  # ...and against the hand-written values, which the stencil then certifies.
  expect_equal(cr[["mu"]], rep(1 / 1.7^2, 5))
  expect_equal(cr[["sigma"]], 2 * (y - 0.3) / 1.7^3)

  for (p in c("mu", "sigma")) {
    expect_equal(cr[[p]], mixed_stencil(d, y, th, p),
      tolerance = 1e-5, label = paste("gaussian", p)
    )
  }
})

test_that("the Student t closed form matches the mixed stencil", {
  d <- student_t1_distrib()
  th <- list(mu = -0.2, sigma = 1.3, nu = 4.5)
  y <- c(-3, -0.2, 0.4, 2.5)

  cr <- distrib_cross_y(d, y, th)
  expect_named(cr, c("mu", "sigma", "nu"))
  for (p in c("mu", "sigma", "nu")) {
    expect_equal(cr[[p]], mixed_stencil(d, y, th, p),
      tolerance = 1e-4, label = paste("student t", p)
    )
  }
})

test_that("the Student t tends to the gaussian as nu grows", {
  y <- c(-1, 0.5, 2)
  tg <- distrib_cross_y(gaussian1_distrib(), y, list(mu = 0, sigma = 1.5))
  tt <- distrib_cross_y(
    student_t1_distrib(), y, list(mu = 0, sigma = 1.5, nu = 1e7)
  )
  expect_equal(tt[["mu"]], tg[["mu"]], tolerance = 1e-5)
  expect_equal(tt[["sigma"]], tg[["sigma"]], tolerance = 1e-5)
})

test_that("the fallback covers a family without a closed form", {
  d <- gamma2_distrib()
  th <- list(mu = 2, sigma2 = 0.8)
  y <- c(0.4, 1.5, 3, 6)

  cr <- distrib_cross_y(d, y, th)
  expect_named(cr, c("mu", "sigma2"))
  for (p in c("mu", "sigma2")) {
    expect_equal(cr[[p]], mixed_stencil(d, y, th, p),
      tolerance = 1e-4, label = paste("gamma", p)
    )
  }
})

test_that("the link scale multiplies each component by its own h'", {
  d <- gaussian1_distrib()
  th <- list(mu = 0.3, sigma = 1.7)
  y <- c(-1, 0.5, 2)

  cr <- distrib_cross_y(d, y, th)
  cl <- distrib_cross_y(d, y, th, scale = "link")

  # mu has the identity link, sigma the log link, so h'(eta) = sigma itself.
  expect_equal(cl[["mu"]], cr[["mu"]])
  expect_equal(cl[["sigma"]], cr[["sigma"]] * 1.7)
})

test_that("truncation leaves the mixed derivatives unchanged", {
  # l_T = l - log Z(theta): the constant has no y, so it drops from any
  # derivative involving the response. Both the delegation and the stencil on
  # the truncated log-density itself must say so.
  g <- gaussian1_distrib()
  d <- truncated(g, lower = -1, upper = 3)
  th <- list(mu = 0.3, sigma = 1.2)
  y <- c(-0.5, 0.8, 2.4)

  cr <- distrib_cross_y(d, y, th)
  expect_equal(cr, distrib_cross_y(g, y, th))
  for (p in c("mu", "sigma")) {
    expect_equal(cr[[p]], mixed_stencil(d, y, th, p),
      tolerance = 1e-5, label = paste("truncated", p)
    )
  }
})

test_that("fixed() keeps the components of the free parameters", {
  g <- gaussian1_distrib()
  d <- fixed(g, mu = 0.3)
  y <- c(-1, 0.5, 2)

  cr <- distrib_cross_y(d, y, list(sigma = 1.7))
  full <- distrib_cross_y(g, y, list(mu = 0.3, sigma = 1.7))
  expect_named(cr, "sigma")
  expect_equal(cr[["sigma"]], full[["sigma"]])

  # ...and on the link scale, which is where the joint estimation runs.
  crl <- distrib_cross_y(d, y, list(sigma = 1.7), scale = "link")
  expect_equal(crl[["sigma"]], full[["sigma"]] * 1.7)
})

test_that("a deliberately wrong closed form would be caught", {
  # The stencil tolerance must be able to see a 5% error, or the agreement
  # tests above prove nothing.
  d <- gaussian1_distrib()
  th <- list(mu = 0.3, sigma = 1.7)
  y <- c(-1, 0.5, 2)
  wrong <- 1.05 * distrib_cross_y(d, y, th)[["sigma"]]
  ref <- mixed_stencil(d, y, th, "sigma")
  expect_gt(max(abs(wrong - ref) / pmax(abs(ref), 1e-8)), 0.04)
})
