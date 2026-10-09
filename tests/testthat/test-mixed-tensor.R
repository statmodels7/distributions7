# The mixed fallbacks (cross_y, cross2_y, cross3_y, grad_y_hess, hess_y_hess)
# are one tensor stencil on the highest quantity a family implements itself,
# never a difference of a difference. Two references: a log-gamma built from
# transformation(), whose closed forms are elementary, and a gaussian defined
# by its density alone, against gaussian1_distrib()'s closed forms.

block_gap <- function(a, e) {
  ev <- unlist(e[names(e)])
  av <- unlist(a[names(e)])
  max(abs(av - ev)) / max(abs(ev))
}

test_that("a log-gamma from transformation() meets its closed forms", {
  # k = 1/sigma2; log p(b) = k log k - lgamma(k) + k b - k e^b, so every order
  # in b beyond the first is -k e^b, and dk/dsigma2 = -k^2. Measured on
  # 2026-10-09: 1.7e-11, 8.6e-9, 1.5e-7, 4.2e-11 and 1.5e-8; the previous
  # fallback gave 3.7e-7, 1.6e-4, 3.3e-3, 3.7e-3 and a factor of 7.2.
  lg <- fixed(transformation(gamma2_distrib(), log_transform()), mu = 1)
  s2 <- 0.2628267
  th <- list(sigma2 = s2)
  k <- 1 / s2
  b <- c(-1.2, -0.4, 0, 0.3, 0.8)
  rel <- function(a, e) max(abs(unlist(a) - e)) / max(abs(e))
  expect_lt(rel(distrib_cross_y(lg, b, th), -k^2 * (1 - exp(b))), 1e-9)
  expect_lt(rel(distrib_cross2_y(lg, b, th), k^2 * exp(b)), 1e-6)
  expect_lt(rel(distrib_cross3_y(lg, b, th), k^2 * exp(b)), 1e-5)
  expect_lt(rel(distrib_grad_y_hess(lg, b, th), 2 * k^3 * (1 - exp(b))), 1e-9)
  expect_lt(rel(distrib_hess_y_hess(lg, b, th), -2 * k^3 * exp(b)), 1e-6)
})

test_that("a density alone meets gaussian1's closed forms", {
  # measured on 2026-10-09 on the scale of each block: 3.0e-8, 1.1e-6,
  # 9.7e-7 and 1.0e-5; cross3_y is zero for a gaussian and is compared
  # absolutely. The previous fallback gave hess_y_hess off by more than 100.
  d0 <- density_only_distrib()
  d1 <- gaussian1_distrib()
  y <- c(-2.1, -0.3, 0.4, 1.7)
  th <- list(mu = 0.3, sigma = 1.4)
  expect_lt(block_gap(distrib_cross_y(d0, y, th), distrib_cross_y(d1, y, th)), 1e-6)
  expect_lt(block_gap(distrib_cross2_y(d0, y, th), distrib_cross2_y(d1, y, th)), 1e-4)
  expect_lt(block_gap(distrib_grad_y_hess(d0, y, th),
                      distrib_grad_y_hess(d1, y, th)), 1e-4)
  expect_lt(block_gap(distrib_hess_y_hess(d0, y, th),
                      distrib_hess_y_hess(d1, y, th)), 1e-3)
  c3 <- unlist(distrib_cross3_y(d0, y, th))
  expect_true(all(abs(c3 - unlist(distrib_cross3_y(d1, y, th))) < 1e-3))
})

test_that("the names and lengths are those of the generics", {
  d0 <- density_only_distrib()
  y <- c(0.1, 0.5, 2)
  th <- list(mu = 0.3, sigma = 1.4)
  expect_identical(names(distrib_cross_y(d0, y, th)), c("mu", "sigma"))
  expect_identical(names(distrib_hess_y_hess(d0, y, th)), hess_names(c("mu", "sigma")))
  expect_true(all(lengths(distrib_grad_y_hess(d0, y, th)) == 3L))
})
