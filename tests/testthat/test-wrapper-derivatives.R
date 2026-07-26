# Third and fourth derivatives of the wrappers, which are assembled from the
# parent's by the set-partition machinery in R/wrapper_derivatives.R.

wrapper_deriv_cases <- function() {
  list(
    transformed = list(d = transformation(gaussian_distrib(), exp_transform()),
                       theta = list(mu = 0.5, sigma = 1.1), y = c(0.6, 1.4, 3.0)),
    zip = list(d = zero_inflated(poisson_distrib()),
               theta = list(mu = 3, zi = 0.25), y = c(0, 1, 4)),
    zinb = list(d = zero_inflated(negbin_distrib()),
                theta = list(mu = 4, theta = 1.5, zi = 0.3), y = c(0, 2, 6)),
    hurdle = list(d = zero_adjusted(poisson_distrib()),
                  theta = list(mu = 3, za = 0.4), y = c(0, 1, 5)),
    zaga = list(d = zero_adjusted(gamma_distrib()),
                theta = list(mu = 3, sigma2 = 2, za = 0.3), y = c(0, 1.2, 5)),
    ztp = list(d = truncated(poisson_distrib(), lower = 1),
               theta = list(mu = 2.5), y = c(1, 3, 6)),
    tnorm = list(d = truncated(gaussian_distrib(), -1, 2),
                 theta = list(mu = 0.5, sigma = 1.5), y = c(-0.5, 0.3, 1.6)),
    tgam = list(d = truncated(gamma_distrib(), 0.5, 8),
                theta = list(mu = 3, sigma2 = 2), y = c(0.8, 2, 6))
  )
}

test_that("wrapper third and fourth derivatives match finite differences", {
  for (nm in names(wrapper_deriv_cases())) {
    case <- wrapper_deriv_cases()[[nm]]
    d <- case$d; th <- case$theta; y <- case$y

    a3 <- distrib_deriv3(d, y, th)
    n3 <- numerical_deriv3(d, y, th)
    expect_setequal(names(a3), names(n3))
    for (k in names(n3)) {
      expect_equal(a3[[k]], n3[[k]], tolerance = 1e-5,
                   label = sprintf("%s deriv3 %s", nm, k))
    }

    a4 <- distrib_deriv4(d, y, th)
    n4 <- numerical_deriv4(d, y, th)
    expect_setequal(names(a4), names(n4))
    for (k in names(n4)) {
      expect_equal(a4[[k]], n4[[k]], tolerance = 1e-4,
                   label = sprintf("%s deriv4 %s", nm, k))
    }
  }
})

test_that("a transformed distribution's higher derivatives are exactly the parent's", {
  # The Jacobian does not depend on theta, so nothing is approximated here.
  g <- gaussian_distrib()
  td <- transformation(g, exp_transform())
  th <- list(mu = 0.5, sigma = 1.1)
  y <- c(0.6, 1.4, 3.0)

  expect_equal(distrib_deriv3(td, y, th), distrib_deriv3(g, log(y), th))
  expect_equal(distrib_deriv4(td, y, th), distrib_deriv4(g, log(y), th))
})

test_that("zero-inflation leaves the parent's derivatives untouched above zero", {
  d <- zero_inflated(poisson_distrib())
  th <- list(mu = 3, zi = 0.25)
  y <- c(2, 5, 9)

  parent3 <- distrib_deriv3(poisson_distrib(), y, list(mu = 3))
  parent4 <- distrib_deriv4(poisson_distrib(), y, list(mu = 3))

  # [[ ]] rather than $: on a list, $ does partial matching, so a wrong key can
  # quietly resolve to a longer one -- or to NULL when the prefix is ambiguous.
  expect_equal(distrib_deriv3(d, y, th)[["mu_mu_mu"]], parent3[["mu_mu_mu"]])
  expect_equal(distrib_deriv4(d, y, th)[["mu_mu_mu_mu"]], parent4[["mu_mu_mu_mu"]])

  # ... and every component mixing theta with zi vanishes above zero
  d3 <- distrib_deriv3(d, y, th)
  for (k in c("mu_mu_zi", "mu_zi_zi")) expect_equal(d3[[k]], rep(0, 3), label = k)
})

test_that("the pure zi derivatives at zero match the closed form", {
  # At y = 0 the log-likelihood is log(f0 + zi(1 - f0)), so the k-th derivative
  # in zi is (-1)^{k-1} (k-1)! ((1-f0)/L0)^k.
  d <- zero_inflated(poisson_distrib())
  mu <- 3; zi <- 0.25
  th <- list(mu = mu, zi = zi)
  f0 <- dpois(0, mu)
  L0 <- zi + (1 - zi) * f0
  r <- (1 - f0) / L0

  expect_equal(distrib_deriv3(d, 0, th)[["zi_zi_zi"]], 2 * r^3)
  expect_equal(distrib_deriv4(d, 0, th)[["zi_zi_zi_zi"]], -6 * r^4)
})

test_that("the hurdle likelihood separates at every order", {
  d <- zero_adjusted(negbin_distrib())
  th <- list(mu = 4, theta = 1.5, za = 0.35)
  y <- c(0, 2, 7)

  for (dk in list(distrib_deriv3(d, y, th), distrib_deriv4(d, y, th))) {
    mixed <- names(dk)[grepl("za", names(dk)) & !grepl("^(za_)+za$", names(dk))]
    expect_gt(length(mixed), 0)
    for (k in mixed) expect_equal(dk[[k]], rep(0, 3), label = paste("mixed", k))
  }
})

test_that("the zero-adjusted continuous derivatives switch off at the atom", {
  d <- zero_adjusted(gamma_distrib())
  th <- list(mu = 3, sigma2 = 2, za = 0.3)
  y <- c(0, 1.2, 5)
  parent <- distrib_deriv3(gamma_distrib(), y, list(mu = 3, sigma2 = 2))

  got <- distrib_deriv3(d, y, th)[["mu_mu_mu"]]
  expect_equal(got[1], 0)                          # the atom carries no theta information
  expect_equal(got[-1], parent[["mu_mu_mu"]][-1])  # elsewhere it is the parent's
})

test_that("the truncated third derivative agrees with one obtained from the cdf", {
  # An independent route: l_T = l - log Z, and log Z is a function of theta that
  # can be differentiated numerically straight from the parent's cdf, touching
  # none of the partition machinery.
  #
  # The bounds must be ASYMMETRIC about mu. With mu at the centre of the interval
  # the truncated law is symmetric, every odd derivative of log Z vanishes, and
  # the comparison degenerates into 0 == 0 -- a test that passes without
  # exercising anything.
  s <- 1.5
  for (b in list(c(-1, 3, 0.5), c(0, 4, 1.0), c(-2, 1, 0.3))) {
    lo <- b[1]; up <- b[2]; mu <- b[3]
    d <- truncated(gaussian_distrib(), lo, up)
    th <- list(mu = mu, sigma = s)
    y <- c(lo + 0.2, (lo + up) / 2, up - 0.3)

    logZ <- function(m) log(pnorm((up - m) / s) - pnorm((lo - m) / s))
    h <- 1e-3
    d3_logZ <- (logZ(mu + 2 * h) - 2 * logZ(mu + h) +
                2 * logZ(mu - h) - logZ(mu - 2 * h)) / (2 * h^3)

    expect_gt(abs(d3_logZ), 1e-3)    # the test would be vacuous otherwise
    expect_equal(distrib_deriv3(d, y, th)[["mu_mu_mu"]],
                 distrib_deriv3(gaussian_distrib(), y, th)[["mu_mu_mu"]] - d3_logZ,
                 tolerance = 1e-5,
                 label = sprintf("truncated Gaussian on (%g, %g)", lo, up))
  }
})

test_that("expected higher derivatives are still available on the wrappers", {
  for (nm in c("zip", "ztp", "zaga")) {
    case <- wrapper_deriv_cases()[[nm]]
    e <- distrib_deriv3(case$d, 0, case$theta, expected = TRUE)
    expect_equal(names(e), deriv_names(case$d@params, 3))
    expect_true(all(is.finite(unlist(e))), label = paste(nm, "expected deriv3 finite"))
  }
})

test_that("check_distrib passes on the wrappers at all four orders", {
  skip_on_cran()
  for (nm in names(wrapper_deriv_cases())) {
    set.seed(77)
    case <- wrapper_deriv_cases()[[nm]]
    out <- check_distrib(case$d, theta = case$theta, n = 30, nsim = 3e4,
                         orders = 1:4, verbose = FALSE)
    expect_equal(out$status, rep("OK", nrow(out)),
                 label = paste(nm, ":", paste(out$check[out$status == "FAIL"], collapse = ", ")))
  }
})

test_that("order_indices lines up with the names deriv_names gives", {
  # The pairing is what attaches a component to its multi-index; a mismatch
  # would label "mu_sigma" with the index (sigma, sigma) and be invisible.
  for (params in list(c("mu", "sigma"), c("mu", "theta", "zi"), "mu")) {
    for (k in 2:4) {
      idxs <- distributions7:::order_indices(params, k)
      nms <- deriv_names(params, k)
      expect_equal(length(idxs), length(nms))
      rebuilt <- vapply(idxs, function(i) paste(i, collapse = "_"), character(1))
      expect_equal(rebuilt, nms, label = sprintf("order %d, %d params", k, length(params)))
    }
  }
})
