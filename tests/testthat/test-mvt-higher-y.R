test_that("the multivariate t's third and fourth response derivatives are exact", {
  # One central difference of the order below, which is analytic, along each
  # response coordinate; the scale matrix is a log-Cholesky and an inverse AR(1).
  for (d in list(mvstudent_t1_distrib(3),
                 mvstudent_t2_distrib(3, omega = parameters7::ar1(3)))) {
    np <- length(d@params)
    th <- stats::setNames(as.list(c(0.5, -0.3, 0.2,
                                    seq(0.1, by = 0.13, length.out = np - 4),
                                    4.5)), d@params)
    set.seed(2)
    y <- distrib_rng(d, 4, th)
    d3 <- distrib_deriv3_y(d, y, th)
    d4 <- distrib_deriv4_y(d, y, th)
    expect_identical(dim(d3), c(3L, 3L, 3L, 4L))
    expect_identical(dim(d4), c(3L, 3L, 3L, 3L, 4L))
    h <- 1e-5
    for (k in 1:3) {
      yp <- y
      yp[, k] <- yp[, k] + h
      ym <- y
      ym[, k] <- ym[, k] - h
      n3 <- (distrib_hess_y(d, yp, th) - distrib_hess_y(d, ym, th)) / (2 * h)
      n4 <- (distrib_deriv3_y(d, yp, th) - distrib_deriv3_y(d, ym, th)) /
        (2 * h)
      expect_lt(max(abs(d3[, , k, ] - n3)), 1e-8)
      expect_lt(max(abs(d4[, , , k, ] - n4)), 1e-8)
    }
    # both are symmetric in every pair of response indices
    expect_lt(max(abs(d3 - aperm(d3, c(2, 1, 3, 4)))), 1e-14)
    expect_lt(max(abs(d3 - aperm(d3, c(3, 2, 1, 4)))), 1e-14)
    expect_lt(max(abs(d4 - aperm(d4, c(4, 2, 3, 1, 5)))), 1e-14)
  }
})

test_that("the multivariate t's mixed derivative is the third one moved in each parameter", {
  d <- mvstudent_t1_distrib(2)
  th <- list(mu1 = 0.5, mu2 = -0.3, sigma_log_L1 = 0.1, sigma_log_L2 = -0.2,
             sigma_L2.1 = 0.4, nu = 6)
  set.seed(1)
  y <- distrib_rng(d, 3, th)
  c3 <- distrib_cross3_y(d, y, th)
  expect_named(c3, d@params)
  h <- 1e-5
  for (a in d@params) {
    tp <- th
    tp[[a]] <- tp[[a]] + h
    tm <- th
    tm[[a]] <- tm[[a]] - h
    num <- (distrib_deriv3_y(d, y, tp) - distrib_deriv3_y(d, y, tm)) / (2 * h)
    expect_lt(max(abs(c3[[a]] - num)), 1e-8)
  }
  # a 5% error in one piece of the assembly is far outside that tolerance
  wrong <- c3$nu * 1.05
  num <- (distrib_deriv3_y(d, y, modifyList(th, list(nu = 6 + h))) -
            distrib_deriv3_y(d, y, modifyList(th, list(nu = 6 - h)))) / (2 * h)
  expect_gt(max(abs(wrong - num)), 1e-4)
})

test_that("fixed() passes the multivariate t's higher response derivatives on", {
  f <- fixed(mvstudent_t1_distrib(2), mu1 = 0, mu2 = 0)
  th <- list(sigma_log_L1 = 0.1, sigma_log_L2 = -0.2, sigma_L2.1 = 0.3, nu = 3)
  full <- c(list(mu1 = 0, mu2 = 0), th)
  y <- matrix(c(1, -0.5, 0.3, 2), 2, 2)
  expect_identical(distrib_deriv3_y(f, y, th),
                   distrib_deriv3_y(mvstudent_t1_distrib(2), y, full))
  expect_identical(distrib_deriv4_y(f, y, th),
                   distrib_deriv4_y(mvstudent_t1_distrib(2), y, full))
  c3 <- distrib_cross3_y(f, y, th)
  expect_named(c3, names(th))
  expect_identical(c3$nu,
                   distrib_cross3_y(mvstudent_t1_distrib(2), y, full)$nu)
})
