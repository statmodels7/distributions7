# Three repairs to the multivariate families found by the measurements of
# 2026-10-17: the length of the expected derivatives on the link scale, the
# fifth derivative of the elliptical families at p = 3, and the zero response
# derivatives of the gaussian.

mv_repair_cases <- function() {
  list(
    mvgaussian1 = list(
      d = mvgaussian1_distrib(2),
      theta = list(mu1 = 0.5, mu2 = -0.3, sigma_log_L1 = 0.1,
                   sigma_log_L2 = -0.2, sigma_L2.1 = 0.4)),
    mvstudent_t1 = list(
      d = mvstudent_t1_distrib(2),
      theta = list(mu1 = 0.5, mu2 = -0.3, sigma_log_L1 = 0.1,
                   sigma_log_L2 = -0.2, sigma_L2.1 = 0.4, nu = 6)),
    dirichlet = list(
      d = dirichlet_distrib(3),
      theta = list(mean_alr1 = 0.3, mean_alr2 = -0.2, phi = 12)),
    multinomial = list(
      d = multinomial_distrib(3, size = 20),
      theta = list(probs_alr1 = 0.3, probs_alr2 = -0.2)))
}

test_that("expected derivatives on the link scale have one value per observation", {
  n <- 5
  for (nm in names(mv_repair_cases())) {
    cs <- mv_repair_cases()[[nm]]
    set.seed(1)
    y <- distrib_rng(cs$d, n, cs$theta)
    for (f in list(distrib_deriv3, distrib_deriv4)) {
      r <- f(cs$d, y, cs$theta, expected = TRUE, scale = "link")
      expect_true(all(lengths(r) == n), info = nm)
    }
  }
  # every link of the gaussian is the identity, so the two scales coincide
  cs <- mv_repair_cases()$mvgaussian1
  set.seed(1)
  y <- distrib_rng(cs$d, n, cs$theta)
  for (f in list(distrib_deriv3, distrib_deriv4)) {
    expect_equal(f(cs$d, y, cs$theta, expected = TRUE, scale = "link"),
                 f(cs$d, y, cs$theta, expected = TRUE, scale = "parameter"))
  }
})

test_that("the fifth derivative of the multivariate gaussian runs at p = 3", {
  skip_on_cran()
  d <- mvgaussian1_distrib(3)
  theta <- list(mu1 = 0.2, mu2 = -0.1, mu3 = 0.4, sigma_log_L1 = 0.1,
                sigma_log_L2 = -0.2, sigma_log_L3 = 0.05, sigma_L2.1 = 0.3,
                sigma_L3.1 = -0.2, sigma_L3.2 = 0.1)
  set.seed(1)
  y <- distrib_rng(d, 3, theta)
  d5 <- distrib_deriv5(d, y, theta)
  expect_length(d5, choose(9 + 4, 5))
  expect_true(all(lengths(d5) == 3))
  # the log-density is quadratic in the mean
  expect_equal(d5$mu1_mu1_mu2_mu3_sigma_L3.2, rep(0, 3), tolerance = 1e-6)
  # one component against a difference of the analytic fourth derivative
  h <- 1e-4
  tp <- theta; tp$sigma_L3.2 <- tp$sigma_L3.2 + h
  tm <- theta; tm$sigma_L3.2 <- tm$sigma_L3.2 - h
  k <- "sigma_log_L1_sigma_log_L2_sigma_L2.1_sigma_L3.1"
  num <- (distrib_deriv4(d, y, tp)[[k]] - distrib_deriv4(d, y, tm)[[k]]) / (2 * h)
  expect_equal(d5[[paste0(k, "_sigma_L3.2")]], num, tolerance = 1e-6)
})

test_that("the gaussian's higher response derivatives are zero arrays", {
  cs <- mv_repair_cases()$mvgaussian1
  d <- cs$d
  set.seed(1)
  y <- distrib_rng(d, 4, cs$theta)
  expect_identical(distrib_deriv3_y(d, y, cs$theta), array(0, c(2, 2, 2)))
  expect_identical(distrib_deriv4_y(d, y, cs$theta), array(0, c(2, 2, 2, 2)))
  for (sc in c("parameter", "link")) {
    c3 <- distrib_cross3_y(d, y, cs$theta, scale = sc)
    expect_named(c3, d@params)
    expect_true(all(vapply(c3, function(a) identical(a, array(0, c(2, 2, 2))),
                           logical(1))))
  }
  # fixed() delegates to the parent
  fx <- fixed(d, mu1 = 0.5, mu2 = -0.3)
  th <- cs$theta[c("sigma_log_L1", "sigma_log_L2", "sigma_L2.1")]
  expect_identical(distrib_deriv3_y(fx, y, th), array(0, c(2, 2, 2)))
  expect_identical(distrib_deriv4_y(fx, y, th), array(0, c(2, 2, 2, 2)))
  expect_named(distrib_cross3_y(fx, y, th), fx@params)
})
