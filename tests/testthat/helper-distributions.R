# Shared helpers for the test suite.

# List of all implemented distributions with a valid interior theta.
all_distrib_cases <- function() {
  list(
    gaussian  = list(d = gaussian_distrib(),  theta = list(mu = 1.5, sigma = 2.0)),
    lognormal = list(d = lognormal_distrib(), theta = list(mu = 0.5, sigma2 = 1.3)),
    gamma     = list(d = gamma_distrib(),     theta = list(mu = 3.0, sigma2 = 2.0)),
    cauchy    = list(d = cauchy_distrib(),    theta = list(mu = 0.5, sigma = 1.4)),
    logistic  = list(d = logistic_distrib(),  theta = list(mu = 0.5, sigma = 1.4)),
    invgauss  = list(d = invgauss_distrib(),  theta = list(mu = 2.0, phi = 0.7)),
    beta      = list(d = beta_distrib(),      theta = list(mu = 0.4, phi = 6.0)),
    student_t = list(d = student_t_distrib(), theta = list(mu = 0.5, sigma = 1.3, nu = 6.0)),
    poisson   = list(d = poisson_distrib(),   theta = list(mu = 4.0)),
    bernoulli = list(d = bernoulli_distrib(), theta = list(mu = 0.35)),
    binomial  = list(d = binomial_distrib(size = 10), theta = list(mu = 0.35)),
    negbin    = list(d = negbin_distrib(),    theta = list(mu = 4.0, theta = 1.7))
  )
}

# Central finite-difference gradient of the log-density w.r.t. each parameter.
fd_gradient_ref <- function(d, y, theta, h = 1e-5) {
  out <- lapply(seq_along(theta), function(j) {
    tp <- tm <- theta
    hh <- h * max(1, abs(theta[[j]]))
    tp[[j]] <- theta[[j]] + hh
    tm[[j]] <- theta[[j]] - hh
    (distrib_pdf(d, y, tp, log = TRUE) - distrib_pdf(d, y, tm, log = TRUE)) / (2 * hh)
  })
  names(out) <- names(theta)
  out
}

# Central finite-difference Hessian of the log-density, named like the
# analytic output ("mu_mu", "mu_sigma", ...).
fd_hessian_ref <- function(d, y, theta, h = 1e-4) {
  p <- length(theta)
  out <- list()
  for (j in 1:p) {
    for (k in j:p) {
      hj <- h * max(1, abs(theta[[j]]))
      hk <- h * max(1, abs(theta[[k]]))
      shift <- function(a, b) {
        t2 <- theta
        t2[[j]] <- t2[[j]] + a * hj
        t2[[k]] <- t2[[k]] + b * hk
        t2
      }
      val <- (distrib_pdf(d, y, shift(1, 1), log = TRUE) -
                distrib_pdf(d, y, shift(1, -1), log = TRUE) -
                distrib_pdf(d, y, shift(-1, 1), log = TRUE) +
                distrib_pdf(d, y, shift(-1, -1), log = TRUE)) / (4 * hj * hk)
      out[[paste0(names(theta)[j], "_", names(theta)[k])]] <- val
    }
  }
  out
}
