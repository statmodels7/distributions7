# The scalar C entry points of the fast route: score and (k, k) second
# derivative of the log-density in one parameter, each mirroring the
# family's own vector kernel expression by expression. identical(), not a
# tolerance: the consumer's twin test rests on these being the same numbers
# to the bit.

test_that("gaussian1's scalar entries mirror its kernels bit for bit", {
  set.seed(1)
  n <- 60
  y <- rnorm(n, 1, 2)
  mu <- runif(n, -1, 1)
  sg <- runif(n, 0.5, 3)
  th <- list(mu = mu, sigma = sg)
  g <- distrib_gradient(gaussian1_distrib(), y, th)
  h <- distrib_hessian(gaussian1_distrib(), y, th)
  p1 <- d7_scalar_probe("Gaussian1Distrib", 1L, y, cbind(mu, sg))
  p2 <- d7_scalar_probe("Gaussian1Distrib", 2L, y, cbind(mu, sg))
  expect_identical(p1$score, g[["mu"]])
  expect_identical(p1$curvature, h[["mu_mu"]])
  expect_identical(p2$score, g[["sigma"]])
  expect_identical(p2$curvature, h[["sigma_sigma"]])
})

test_that("gamma1's scalar entries mirror its kernels bit for bit", {
  set.seed(2)
  n <- 60
  y <- rgamma(n, shape = 2, scale = 1.5)
  mu <- runif(n, 1, 5)
  ph <- runif(n, 0.2, 1.5)
  th <- list(mu = mu, phi = ph)
  g <- distrib_gradient(gamma1_distrib(), y, th)
  h <- distrib_hessian(gamma1_distrib(), y, th)
  p1 <- d7_scalar_probe("Gamma1Distrib", 1L, y, cbind(mu, ph))
  p2 <- d7_scalar_probe("Gamma1Distrib", 2L, y, cbind(mu, ph))
  expect_identical(p1$score, g[["mu"]])
  expect_identical(p1$curvature, h[["mu_mu"]])
  expect_identical(p2$score, g[["phi"]])
  expect_identical(p2$curvature, h[["phi_phi"]])
})

test_that("an unknown family answers -1", {
  pr <- d7_scalar_probe("NoSuchDistrib", 1L, c(1, 2),
                        cbind(c(1, 1), c(1, 1)))
  expect_identical(pr$id, -1L)
})
