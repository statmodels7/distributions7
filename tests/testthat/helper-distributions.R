# Shared helpers for the test suite.

# List of all implemented distributions with a valid interior theta.
all_distrib_cases <- function() {
  list(
    gaussian  = list(d = gaussian1_distrib(),  theta = list(mu = 1.5, sigma = 2.0)),
    lognormal = list(d = lognormal1_distrib(), theta = list(mu = 0.5, sigma2 = 1.3)),
    gamma     = list(d = gamma2_distrib(),     theta = list(mu = 3.0, sigma2 = 2.0)),
    cauchy    = list(d = cauchy_distrib(),    theta = list(mu = 0.5, sigma = 1.4)),
    logistic  = list(d = logistic_distrib(),  theta = list(mu = 0.5, sigma = 1.4)),
    invgauss  = list(d = invgauss1_distrib(),  theta = list(mu = 2.0, phi = 0.7)),
    beta      = list(d = beta1_distrib(),      theta = list(mu = 0.4, phi = 6.0)),
    student_t = list(d = student_t1_distrib(), theta = list(mu = 0.5, sigma = 1.3, nu = 6.0)),
    poisson   = list(d = poisson_distrib(),   theta = list(mu = 4.0)),
    bernoulli = list(d = bernoulli_distrib(), theta = list(mu = 0.35)),
    binomial  = list(d = binomial_distrib(size = 10), theta = list(mu = 0.35)),
    negbin    = list(d = negbin2_distrib(),    theta = list(mu = 4.0, theta = 1.7))
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


# What a fit actually did, as one line. A test that asserts convergence and
# fails prints only FALSE, which says nothing about whether the run stopped a
# hair above its tolerance or never approached it; passing this as `info` puts
# the number in the report, on whichever platform the failure happens.
fit_report <- function(fit, distrib, y) {
  sc <- tryCatch(
    max(abs(vapply(
      distrib_gradient(distrib, y, as.list(coef(fit)), scale = "link"),
      sum, numeric(1)
    ))) / n_obs(distrib, y),
    error = function(e) NA_real_
  )
  sprintf("method %s, %d iterations, criterion '%s', note '%s', score/n %.3e",
          fit@method, fit@iterations, fit@criterion, fit@note, sc)
}

# A discrete family that APPROXIMATES its expected information: the mass, the
# score, the Hessian and the draws of pig1_distrib(), and no expected method of
# its own, so the base class's strategies apply. Every shipped family now
# computes its expected information exactly, and the approximation machinery
# (approx =, the refusals that read expected_hessian_exact(), the outer
# product) is tested on this one.
PigBare <- S7::new_class("PigBare", parent = discrete_distrib, package = NULL)
S7::method(distrib_pdf, PigBare) <- function(distrib, y, theta, log = FALSE, ...) {
  distrib_pdf(pig1_distrib(), y, theta, log = log)
}
S7::method(distrib_gradient, PigBare) <- function(distrib, y, theta,
                                                  scale = c("parameter", "link"), ...) {
  distrib_gradient(pig1_distrib(), y, theta)
}
S7::method(distrib_hessian, PigBare) <- function(distrib, y, theta,
                                                 scale = c("parameter", "link"), ...) {
  distrib_hessian(pig1_distrib(), y, theta)
}
S7::method(distrib_rng, PigBare) <- function(distrib, n, theta, ...) {
  distrib_rng(pig1_distrib(), n, theta)
}
pig_bare_distrib <- function() {
  p <- pig1_distrib()
  PigBare(distrib_name = "pig bare", dimension = "univariate",
          bounds = c(0, Inf), params = c("mu", "sigma"), n_params = 2,
          params_bounds = p@params_bounds, link_params = p@link_params,
          params_smooth = c(mu = TRUE, sigma = TRUE),
          params_interpretation = p@params_interpretation)
}
