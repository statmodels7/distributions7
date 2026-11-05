test_that("expect_mc() simulates each distinct configuration once", {
  g <- gaussian1_distrib()
  mu <- c(0, 1, 0, 0, 2, 1)
  set.seed(11)
  e <- distrib_expected_hessian(g, rep(0, 6), list(mu = mu, sigma = 2),
                                method = expect_mc(nsim = 40))
  after <- stats::runif(1)
  # the reference: one simulation per distinct row, in order of appearance
  set.seed(11)
  ref <- lapply(c(0, 1, 2), function(m) {
    ys <- stats::rnorm(40, m, 2)
    h <- distrib_hessian(g, ys, list(mu = m, sigma = 2))
    vapply(h, base::mean, numeric(1))
  })
  expect_identical(stats::runif(1), after)
  ref <- do.call(rbind, ref)[match(mu, c(0, 1, 2)), ]
  for (nm in names(e)) expect_identical(e[[nm]], unname(ref[, nm]))
  expect_identical(e$sigma_sigma[c(1, 3, 4)], rep(e$sigma_sigma[1], 3))
})

test_that("expect_mc() with distinct rows consumes the stream row by row", {
  g <- gamma1_distrib()
  th <- list(mu = c(1, 2, 3), phi = c(0.5, 0.4, 0.5))
  set.seed(3)
  e <- distrib_expected_hessian(g, rep(1, 3), th, method = expect_mc(nsim = 25))
  set.seed(3)
  ref <- vapply(1:3, function(i) {
    t_i <- list(mu = th$mu[i], phi = th$phi[i])
    ys <- distrib_rng(g, 25, t_i)
    vapply(distrib_hessian(g, ys, t_i), base::mean, numeric(1))
  }, numeric(3))
  for (k in seq_along(e)) expect_identical(e[[k]], unname(ref[k, ]))
})

test_that("a binomial with a size per observation has the expected information by every method", {
  b <- binomial_distrib(size = c(5, 10, 20, 10))
  th <- list(mu = c(0.2, 0.5, 0.7, 0.5))
  closed <- distrib_expected_hessian(b, c(1, 2, 3, 4), th)$mu_mu
  int_h <- distrib_expected_hessian(b, c(1, 2, 3, 4), th, method = expect_integrate())$mu_mu
  int_s <- distrib_expected_hessian(b, c(1, 2, 3, 4), th,
                                    method = expect_integrate(identity = "score"))$mu_mu
  expect_equal(int_h, closed, tolerance = 1e-12)
  expect_equal(int_s, closed, tolerance = 1e-12)
  # a scalar mu with a vector size
  cl1 <- distrib_expected_hessian(b, c(1, 2, 3, 4), list(mu = 0.3))$mu_mu
  expect_equal(distrib_expected_hessian(b, c(1, 2, 3, 4), list(mu = 0.3),
                                        method = expect_integrate())$mu_mu,
               cl1, tolerance = 1e-12)
  # Monte Carlo: rows 2 and 4 coincide (size 10, mu 0.5) and share a draw
  set.seed(5)
  mc <- distrib_expected_hessian(b, c(1, 2, 3, 4), th, method = expect_mc(nsim = 30))$mu_mu
  set.seed(5)
  ref <- vapply(1:3, function(i) {
    s_i <- c(5, 10, 20)[i]; m_i <- c(0.2, 0.5, 0.7)[i]
    ys <- stats::rbinom(30, s_i, m_i)
    base::mean(-ys / m_i^2 - (s_i - ys) / (1 - m_i)^2)
  }, numeric(1))
  expect_equal(mc, ref[c(1, 2, 3, 2)], tolerance = 1e-14)
  expect_identical(mc[2], mc[4])
})

test_that("a kink root that is a double is found exactly", {
  kd <- kink_decomposition(laplace_distrib())
  set.seed(1)
  mus <- c(0, 1, -1, pi, 1e-3, 123.456, stats::rnorm(200, sd = 10))
  got <- vapply(mus, function(m)
    distributions7:::.kink_roots(kd, list(mu = m, sigma = 1.3), m - 40, m + 30), numeric(1))
  expect_identical(got, mus)
})
