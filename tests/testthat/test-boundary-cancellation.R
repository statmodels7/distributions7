## Derivatives at a boundary limit, where the direct form cancels.
##
## Several families carry a parameter whose derivative vanishes as the family
## tends to a simpler one -- the gamma to a normal as its dispersion goes to
## zero, the negative binomial to the Poisson, the beta-binomial to the
## binomial, the Student t to the gaussian.  Each derivative is then written
## as a difference of terms that agree to leading order, and the direct form
## loses every digit it has before the values a link can produce.  The
## rewrites live in `src/psi_diff.h` and in the families' own kernels.
##
## Every case below asserts the same two things, and the SECOND is what keeps
## the first honest: the shipped derivative tracks the leading behaviour of
## the quantity it computes, and the form it replaced does NOT -- so reverting
## to the direct expression fails the test rather than passing it silently.
## The limits and thresholds are transcribed from the expansions and from the
## measurements, and share no arithmetic with the implementations.

test_that("the gamma's dispersion score survives a vanishing dispersion", {
  d <- gamma1_distrib()
  ## At y = mu the data term log(z) - z + 1 is exactly zero and the score is
  ##   -[log(s) - psi(s)]/phi^2,  s = 1/phi,
  ## with log(s) - psi(s) = 1/(2s) + 1/(12 s^2) - ...
  asym <- function(phi) {
    s <- 1 / phi
    -(1 / (2 * s) + 1 / (12 * s * s)) / phi^2
  }
  direct <- function(phi) -(log(1 / phi) - digamma(1 / phi)) / phi^2
  for (phi in c(1e-8, 1e-10, 1e-12, 1e-14)) {
    got <- distrib_gradient(d, 3, list(mu = 3, phi = phi))[["phi"]]
    expect_true(is.finite(got))
    expect_equal(got, asym(phi), tolerance = 1e-10)
  }
  ## and the direct form does not: 0.2 per cent out at phi = 1e-12, and
  ## EXACTLY ZERO at 1e-14, where the value is -5e13.
  expect_gt(abs(direct(1e-12) - asym(1e-12)) / abs(asym(1e-12)), 1e-4)
  expect_identical(direct(1e-14), 0)

  ## the orders above it stay finite and keep their signs
  for (phi in c(1e-6, 1e-10, 1e-14)) {
    th <- list(mu = 3, phi = phi)
    expect_gt(distrib_hessian(d, 3, th)[["phi_phi"]], 0)
    expect_lt(distrib_deriv3(d, 3, th)[["phi_phi_phi"]], 0)
    expect_gt(distrib_deriv4(d, 3, th)[["phi_phi_phi_phi"]], 0)
  }
})

test_that("the negative binomial's dispersion score survives the Poisson limit", {
  d <- negbin2_distrib()
  y <- 3
  mu <- 4
  ## dl/dtheta = [y - (y - mu)^2]/(2 theta^2) + O(theta^-3)
  asym <- function(th) (y - (y - mu)^2) / (2 * th^2)
  direct <- function(th) {
    (digamma(y + th) - digamma(th)) + log(th / (th + mu)) + (mu - y) / (th + mu)
  }
  for (th in c(1e5, 1e6, 1e7, 1e9)) {
    got <- distrib_gradient(d, y, list(mu = mu, theta = th))[["theta"]]
    expect_true(is.finite(got))
    expect_equal(got, asym(th), tolerance = 1e-4)
  }
  ## the direct form is 20 per cent out at theta = 1e7 and has changed SIGN
  ## by 1e8, where the value is +1e-16 and it reads -7e-16.
  expect_gt(abs(direct(1e7) - asym(1e7)) / abs(asym(1e7)), 0.1)
  expect_lt(direct(1e8) * asym(1e8), 0)

  ## and it stays finite at the largest theta a log link can produce
  for (gen in list(distrib_gradient, distrib_hessian)) {
    v <- gen(d, y, list(mu = mu, theta = 1e308))
    expect_true(all(vapply(v, is.finite, logical(1))))
  }
})

test_that("psi_shift_diff matches the exact recurrence at every order", {
  ## For an INTEGER shift the difference is a one-signed sum,
  ##   psi^(n)(x+k) - psi^(n)(x) = (-1)^n n! sum_{j<k} 1/(x+j)^(n+1),
  ## which cannot cancel and shares no arithmetic with the asymptotic series.
  exact <- function(n, k, x) {
    if (k == 0) return(0)
    (-1)^n * factorial(n) * sum(1 / (x + 0:(k - 1))^(n + 1))
  }
  direct <- function(n, k, x) {
    if (n == 0L) digamma(x + k) - digamma(x) else psigamma(x + k, n) - psigamma(x, n)
  }
  for (n in 0:3) {
    for (x in c(1e2, 1e6, 1e10, 1e12)) {
      ex <- exact(n, 3, x)
      expect_equal(psi_shift_diff(n, 3, x), ex, tolerance = 1e-14)
    }
    ## and the direct form does not, from about x = 1e8 upward
    expect_gt(abs(direct(n, 3, 1e12) - exact(n, 3, 1e12)) / abs(exact(n, 3, 1e12)),
              1e-5)
  }
  ## a large shift, where the recurrence is dear and the series is not
  for (n in c(0L, 1L)) {
    expect_equal(psi_shift_diff(n, 500, 1e9), exact(n, 500, 1e9),
                 tolerance = 1e-14)
  }
  ## zero at a zero shift, by either branch
  expect_identical(psi_shift_diff(0L, 0, 3), 0)
  expect_identical(psi_shift_diff(2L, 0, 1e6), 0)
  ## and it is vectorized in both arguments
  expect_equal(psi_shift_diff(1L, c(1, 2, 3), 1e6),
               vapply(1:3, function(k) exact(1L, k, 1e6), numeric(1)),
               tolerance = 1e-14)
})

test_that("NB1's dispersion derivatives survive the Poisson limit", {
  ## Here the limit is theta -> 0, where the size r = mu/theta runs away.
  d <- negbin1_distrib()
  y <- 3
  mu <- 4
  lim_mu <- y / mu - 1
  lim_th <- ((y - mu)^2 - y) / (2 * mu)
  ## In the mean the approach is O(theta), which is a stronger statement
  ## than any fixed tolerance: measured 1.2e-04, 1.2e-06, 1.2e-08, 1.4e-10.
  ##
  ## In the dispersion it is O(theta) only down to about 1e-6 and then meets
  ## a floor that RISES again -- 8e-08, 1.2e-07, 7.6e-06, 1.9e-03 at 1e-6 to
  ## 1e-12 -- because the score is P * (-mu/theta^2) + Q, two terms of size
  ## 1e+10 at theta = 1e-10 summing to 0.25.  No accuracy in P closes that;
  ## it would need the two combined symbolically, as the two logarithms
  ## inside P already are.  What is asserted is what holds: five digits over
  ## the whole range, where the direct form is out by a factor of 400 at
  ## theta = 1e-8 and of the wrong sign at 1e-10.
  for (th in c(1e-4, 1e-6, 1e-8, 1e-10)) {
    g <- distrib_gradient(d, y, list(mu = mu, theta = th))
    expect_true(is.finite(g[["mu"]]) && is.finite(g[["theta"]]))
    expect_lt(abs(g[["mu"]] - lim_mu) / abs(lim_mu), 20 * th)
    expect_lt(abs(g[["theta"]] - lim_th) / abs(lim_th), 1e-5)
  }
  ## Written directly the score in theta is a factor of 400 out at
  ## theta = 1e-8 and of the WRONG SIGN at 1e-10.
  direct_th <- function(th) {
    r <- mu / th
    om <- 1 + th
    P <- digamma(y + r) - digamma(r) - log(om)
    Q <- -r / om + y / th - y / om
    P * (-mu / th^2) + Q
  }
  expect_gt(abs(direct_th(1e-8) - lim_th) / abs(lim_th), 10)
  expect_lt(direct_th(1e-10) * lim_th, 0)

  ## the observed Hessian keeps its sign and its Poisson limit, -1/mu
  ## in the mean, at every theta above
  for (th in c(1e-4, 1e-8, 1e-12)) {
    h <- distrib_hessian(d, y, list(mu = mu, theta = th))
    expect_true(all(vapply(h, is.finite, logical(1))))
    expect_lt(abs(h[["mu_mu"]] + y / mu^2) / (y / mu^2), 20 * th)
  }
})

test_that("the beta-binomial's score survives the binomial limit", {
  ## betabinom1 is (mu, sigma) with alpha = mu/sigma and beta = (1-mu)/sigma,
  ## so the concentration is S = 1/sigma and the binomial limit is sigma -> 0,
  ## where the score is the binomial's, (y - n mu)/(mu (1 - mu)).
  d <- betabinom1_distrib(size = 10)
  y <- 3
  n <- 10
  p <- 0.5
  lim <- (y - n * p) / (p * (1 - p))
  expect_equal(lim, -8)
  ## the approach is O(sigma), which is a stronger statement than any fixed
  ## tolerance: measured 9e-06, 9e-08, 1e-09, 4e-13 at the four below.
  for (sg in c(1e-6, 1e-8, 1e-10, 1e-12)) {
    got <- distrib_gradient(d, y, list(mu = p, sigma = sg))[["mu"]]
    expect_true(is.finite(got))
    expect_lt(abs(got - lim) / abs(lim), 20 * sg)
  }
  ## the chain from the shapes written directly divides the digamma
  ## differences by sigma, so their loss is amplified: 4e-08 out at
  ## sigma = 1e-8 and 3.6e-04 at 1e-12.
  direct <- function(sg) {
    S <- 1 / sg
    A <- p * S
    B <- (1 - p) * S
    dS <- digamma(n + S) - digamma(S)
    ((digamma(y + A) - digamma(A) - dS) -
       (digamma(n - y + B) - digamma(B) - dS)) / sg
  }
  expect_gt(abs(direct(1e-12) - lim) / abs(lim), 1e-5)
})

test_that("the Student t's expected information in nu keeps its sign", {
  d <- student_t1_distrib()
  ## E[l_nu_nu] = -7/(2 nu^4) + 13/nu^5 - ...  : strictly negative, and the
  ## direct difference of trigammas reads POSITIVE from nu = 3.2e5.
  asym <- function(nu) -3.5 / nu^4
  for (nu in c(1e5, 1e6, 1e8, 1e10)) {
    e <- distrib_expected_hessian(
      d, 0, list(mu = 0, sigma = 1, nu = nu)
    )[["nu_nu"]]
    expect_true(is.finite(e))
    expect_lt(e, 0)
    expect_equal(e, asym(nu), tolerance = 1e-3)
  }
  direct <- function(nu) {
    (trigamma((nu + 1) / 2) - trigamma(nu / 2)) / 4 +
      (nu + 5) / (2 * nu * (nu + 1) * (nu + 3))
  }
  expect_gt(direct(1e6), 0)  # the sign the rewrite exists to keep

  ## and nothing overflows at the nu the log link can produce
  g <- distrib_gradient(d, 0.5, list(mu = 0, sigma = 1, nu = 1e300))
  expect_true(all(vapply(g, is.finite, logical(1))))
})

test_that("the multivariate t's nu derivatives survive the gaussian limit", {
  d <- mvstudent_t1_distrib(2)
  y <- matrix(c(0.4, -0.2), 1)
  nms <- d@params
  inu <- length(nms)
  th <- as.list(c(0, 0, log(1), log(1), 0, 8))
  names(th) <- nms
  for (nu in c(1e5, 1e7, 1e9)) {
    th[[inu]] <- nu
    g <- distrib_gradient(d, y, th)[[nms[inu]]]
    h <- distrib_hessian(d, y, th)[[paste0(nms[inu], "_", nms[inu])]]
    expect_true(is.finite(g) && is.finite(h))
    expect_lt(h, 0)
  }
})

test_that("the Student t's third derivatives survive its own chart", {
  ## Every component divided by D^3 with D = nu sigma^2 + r^2, and D^3
  ## overflows at D = 5.6e+102 while the log link reaches 1.8e+308: at a nu
  ## the family's own chart can produce the whole surface came back NaN,
  ## which is where statmodels7's exact outer gradient reads. 0.31.0 made
  ## the score and the observed Hessian finite there and did not reach
  ## orders three and four.
  d <- student_t1_distrib()
  y <- c(-2.3, -0.4, 0, 0.7, 3.1)
  m <- 0.2
  s <- 1
  for (nu in c(1e50, 1e150, 1e300, .Machine$double.xmax)) {
    th <- list(mu = rep(m, length(y)), sigma = rep(s, length(y)),
               nu = rep(nu, length(y)))
    v <- distrib_deriv3(d, y, th)
    expect_true(all(vapply(v, function(q) all(is.finite(q)), logical(1))),
                info = paste("nu =", nu))
  }

  ## and it converges on its own closed limit, which shares no arithmetic
  ## with the kernel: nu_nu_nu -> 1.5 (1 + 2 z^2 - z^4) / nu^4
  z2 <- ((m - y) / s)^2
  lim <- sum(1.5 * (1 + 2 * z2 - z2^2))
  prev <- Inf
  for (nu in c(1e6, 1e8, 1e10)) {
    th <- list(mu = rep(m, length(y)), sigma = rep(s, length(y)),
               nu = rep(nu, length(y)))
    got <- sum(distrib_deriv3(d, y, th)[["nu_nu_nu"]])
    e <- abs(got - lim / nu^4) / abs(lim / nu^4)
    expect_lt(e, 1e-3)
    expect_lt(e, prev)          # and it converges, as 1/nu
    prev <- e
  }

  ## the form it replaces cannot: D^3 is already infinite there
  Dcube <- function(nu, s, r) (nu * s^2 + r^2)^3
  expect_true(is.finite(Dcube(1e50, 1, 1)))
  expect_false(is.finite(Dcube(1e150, 1, 1)))
})
