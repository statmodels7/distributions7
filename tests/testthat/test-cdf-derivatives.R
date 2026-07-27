# Derivatives of the distribution function with respect to the parameters.
#
# Each family is checked against the route the implementation does NOT use, so
# that no comparison is tautological:
#
#   continuous  implemented by differencing the cdf   -> checked against the
#                                                        partial-expectation integral
#   discrete    implemented by the partial sum        -> checked against finite
#                                                        differences of the cdf
#   closed form (location-scale) -> checked against both

cdf_deriv_cases <- function() {
  list(
    gaussian = list(d = gaussian_distrib(), th = list(mu = 1.2, sigma = 1.7), q = c(0, 1.2, 3)),
    logistic = list(d = logistic_distrib(), th = list(mu = 0.5, sigma = 1.4), q = c(-1, 0.5, 2)),
    gamma    = list(d = gamma_distrib(),    th = list(mu = 3, sigma2 = 2),    q = c(1, 3, 7)),
    beta     = list(d = beta_distrib(),     th = list(mu = 0.4, phi = 6),     q = c(0.2, 0.5, 0.8)),
    poisson  = list(d = poisson_distrib(),  th = list(mu = 4),                q = c(1, 4, 9)),
    negbin   = list(d = negbin_distrib(),   th = list(mu = 4, theta = 1.7),   q = c(1, 4, 9)),
    binomial = list(d = binomial_distrib(size = 10), th = list(mu = 0.35),    q = c(2, 4, 7))
  )
}

# dF/dtheta_i = int_{-inf}^{q} f * score_i, computed here rather than taken
# from the package.
partial_score <- function(d, q, th, p) {
  vapply(q, function(qq) {
    if (S7::S7_inherits(d, discrete_distrib)) {
      grid <- seq(d@bounds[1], qq)
      sum(distrib_pdf(d, grid, th) * distrib_gradient(d, grid, th)[[p]])
    } else {
      stats::integrate(function(y) distrib_pdf(d, y, th) * distrib_gradient(d, y, th)[[p]],
                       d@bounds[1], qq, rel.tol = 1e-10)$value
    }
  }, numeric(1))
}

fd_cdf <- function(d, q, th, p, h = 1e-5) {
  j <- match(p, d@params)
  hh <- h * max(1, abs(th[[j]]))
  tp <- tm <- th
  tp[[j]] <- th[[j]] + hh
  tm[[j]] <- th[[j]] - hh
  (distrib_cdf(d, q, tp) - distrib_cdf(d, q, tm)) / (2 * hh)
}

test_that("continuous cdf gradients match the partial expectation of the score", {
  for (nm in c("gaussian", "logistic", "gamma", "beta")) {
    case <- cdf_deriv_cases()[[nm]]
    d <- case$d; th <- case$th; q <- case$q
    a <- distrib_grad_cdf(d, q, th, log = FALSE)
    for (p in d@params) {
      expect_equal(a[[p]], partial_score(d, q, th, p), tolerance = 1e-6,
                   label = sprintf("%s dF/d%s", nm, p))
    }
  }
})

test_that("discrete cdf gradients match finite differences of the cdf", {
  # The implementation sums f * score over the support, so differencing the cdf
  # is the independent route here.
  for (nm in c("poisson", "negbin", "binomial")) {
    case <- cdf_deriv_cases()[[nm]]
    d <- case$d; th <- case$th; q <- case$q
    a <- distrib_grad_cdf(d, q, th, log = FALSE)
    for (p in d@params) {
      expect_equal(a[[p]], fd_cdf(d, q, th, p), tolerance = 1e-5,
                   label = sprintf("%s dF/d%s", nm, p))
    }
  }
})

test_that("the location-scale closed forms are the density and its first moment", {
  # dF/dmu = -f(y),  dF/dsigma = -z f(y)
  for (nm in c("gaussian", "logistic", "cauchy", "laplace")) {
    d <- switch(nm, gaussian = gaussian_distrib(), logistic = logistic_distrib(),
                cauchy = cauchy_distrib(), laplace = laplace_distrib())
    th <- stats::setNames(list(1.2, 1.7), d@params)
    q <- c(-1, 0.4, 1.2, 3.5)
    z <- (q - th[[1]]) / th[[2]]
    f <- distrib_pdf(d, q, th)

    a <- distrib_grad_cdf(d, q, th, log = FALSE)
    expect_equal(a[[1]], -f, label = paste(nm, "dF/dmu"))
    expect_equal(a[[2]], -z * f, label = paste(nm, "dF/dscale"))

    # ... and they agree with the partial expectation, which knows none of this
    expect_equal(a[[1]], partial_score(d, q, th, d@params[1]), tolerance = 1e-6,
                 label = paste(nm, "against the partial score"))
  }
})

test_that("the log scale and the upper tail are what they claim", {
  for (nm in names(cdf_deriv_cases())) {
    case <- cdf_deriv_cases()[[nm]]
    d <- case$d; th <- case$th; q <- case$q
    Fq <- distrib_cdf(d, q, th)
    raw <- distrib_grad_cdf(d, q, th, log = FALSE)

    lo <- distrib_grad_cdf(d, q, th, lower.tail = TRUE, log = TRUE)
    up <- distrib_grad_cdf(d, q, th, lower.tail = FALSE, log = TRUE)
    for (p in d@params) {
      expect_equal(lo[[p]], raw[[p]] / Fq, label = sprintf("%s log lower %s", nm, p))
      expect_equal(up[[p]], -raw[[p]] / (1 - Fq), label = sprintf("%s log upper %s", nm, p))
    }
  }
})

test_that("cdf hessians match finite differences of the log-cdf", {
  for (nm in names(cdf_deriv_cases())) {
    case <- cdf_deriv_cases()[[nm]]
    d <- case$d; th <- case$th; q <- case$q
    a <- distrib_hess_cdf(d, q, th)
    for (k in hess_names(d@params)) {
      pr <- distributions7:::hess_pairs(d@params)[[k]]
      i <- match(pr[1], d@params); j <- match(pr[2], d@params)
      hi <- 1e-4 * max(1, abs(th[[i]])); hj <- 1e-4 * max(1, abs(th[[j]]))
      sh <- function(u, v) {
        t2 <- th; t2[[i]] <- t2[[i]] + u * hi; t2[[j]] <- t2[[j]] + v * hj; t2
      }
      num <- (log(distrib_cdf(d, q, sh(1, 1))) - log(distrib_cdf(d, q, sh(1, -1))) -
              log(distrib_cdf(d, q, sh(-1, 1))) + log(distrib_cdf(d, q, sh(-1, -1)))) /
             (4 * hi * hj)
      expect_equal(a[[k]], num, tolerance = 1e-4, label = sprintf("%s d2 log F %s", nm, k))
    }
  }
})

test_that("a censored likelihood gets the right score", {
  # The motivating use: an observation known only to exceed c contributes
  # log(1 - F(c)), so its score is distrib_grad_cdf(lower.tail = FALSE).
  set.seed(1)
  d <- gaussian_distrib()
  th <- list(mu = 1, sigma = 2)
  y <- distrib_rng(d, 200, th)
  cens <- y > 2.5
  y[cens] <- 2.5
  expect_true(any(cens) && any(!cens))

  loglik <- function(mu, sg) {
    t2 <- list(mu = mu, sigma = sg)
    sum(distrib_pdf(d, y[!cens], t2, log = TRUE)) +
      sum(log(distrib_cdf(d, y[cens], t2, lower.tail = FALSE)))
  }
  score <- function(p) {
    sum(distrib_gradient(d, y[!cens], th)[[p]]) +
      sum(distrib_grad_cdf(d, y[cens], th, lower.tail = FALSE, log = TRUE)[[p]])
  }
  h <- 1e-5
  expect_equal(score("mu"), (loglik(1 + h, 2) - loglik(1 - h, 2)) / (2 * h),
               tolerance = 1e-6)
  expect_equal(score("sigma"), (loglik(1, 2 + h) - loglik(1, 2 - h)) / (2 * h),
               tolerance = 1e-6)
})

test_that("interval censoring assembles from the unlogged derivatives", {
  # log(F(b) - F(a)) has score (dF(b) - dF(a)) / (F(b) - F(a)).
  d <- gamma_distrib()
  th <- list(mu = 3, sigma2 = 2)
  a <- 1.5; b <- 4.0

  ga <- distrib_grad_cdf(d, a, th, log = FALSE)
  gb <- distrib_grad_cdf(d, b, th, log = FALSE)
  P <- distrib_cdf(d, b, th) - distrib_cdf(d, a, th)

  ll <- function(mu, s2) {
    t2 <- list(mu = mu, sigma2 = s2)
    log(distrib_cdf(d, b, t2) - distrib_cdf(d, a, t2))
  }
  h <- 1e-5
  expect_equal((gb$mu - ga$mu) / P, (ll(3 + h, 2) - ll(3 - h, 2)) / (2 * h),
               tolerance = 1e-6)
  expect_equal((gb$sigma2 - ga$sigma2) / P, (ll(3, 2 + h) - ll(3, 2 - h)) / (2 * h),
               tolerance = 1e-6)
})

test_that("cdf derivatives are correct with vectorized parameters", {
  d <- gaussian_distrib()
  mu <- c(0, 1, 2)
  th <- list(mu = mu, sigma = 1.5)
  q <- c(0.5, 1.0, 1.5)

  a <- distrib_grad_cdf(d, q, th, log = FALSE)
  for (i in 1:3) {
    one <- distrib_grad_cdf(d, q[i], list(mu = mu[i], sigma = 1.5), log = FALSE)
    expect_equal(a$mu[i], one$mu, tolerance = 1e-10)
    expect_equal(a$sigma[i], one$sigma, tolerance = 1e-10)
  }
})
