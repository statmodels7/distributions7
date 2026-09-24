# Parameters that vary by observation, on the beta-binomial routes that sum
# over the support. Each summed the support against a vector of parameters,
# recycling one against the other: betabinom2's expected information came back
# with warnings and wrong numbers, and betabinom1's distribution function and
# quantile read past the end of the parameter vectors inside the kernel.

test_that("betabinom2's expected derivatives are per observation", {
  d <- betabinom2_distrib(size = 10)
  n <- 7
  th <- list(alpha = seq(1.5, 2.5, length.out = n),
             beta = seq(2.5, 3.5, length.out = n))
  y <- rep(0, n)
  for (f in list(function(t, yy) distrib_expected_hessian(d, yy, t),
                 function(t, yy) distrib_deriv3(d, yy, t, expected = TRUE),
                 function(t, yy) distrib_deriv4(d, yy, t, expected = TRUE))) {
    v <- expect_silent(f(th, y))
    for (i in seq_len(n)) {
      s <- f(list(alpha = th$alpha[i], beta = th$beta[i]), 0)
      for (k in names(s)) expect_identical(v[[k]][i], s[[k]])
    }
  }
})

test_that("betabinom1's distribution function and quantile are per observation", {
  d <- betabinom1_distrib(size = 10)
  th <- list(mu = c(0.2, 0.3, 0.45, 0.6, 0.7), sigma = c(0.1, 0.3, 0.2, 0.5, 0.25))
  q <- c(2, 3, 4, 5, 6)
  p <- c(0.1, 0.3, 0.5, 0.7, 0.9)
  cv <- expect_silent(distrib_cdf(d, q, th))
  qv <- expect_silent(distrib_quantile(d, p, th))
  for (i in seq_along(q)) {
    ti <- list(mu = th$mu[i], sigma = th$sigma[i])
    expect_equal(cv[i], distrib_cdf(d, q[i], ti), tolerance = 1e-14)
    expect_identical(qv[i], distrib_quantile(d, p[i], ti))
  }
  # a scalar quantile against vector parameters is recycled rather than read
  # once
  expect_length(distrib_cdf(d, 3, th), 5L)
})
