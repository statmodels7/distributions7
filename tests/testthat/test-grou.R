# Generalized Ratio-of-Uniforms sampling: the default RNG for continuous
# distributions that provide neither a native generator nor a quantile function.

# Distributions defined by nothing but their density. The classes and their
# methods are created once, at the top level, so that re-running the helpers does
# not emit "Overwriting method" messages.
id_link <- linkfunctions7::identity_link()
lg_link <- linkfunctions7::log_link()

GrouNorm <- S7::new_class("GrouNorm", parent = continuous_distrib, package = NULL)
S7::method(distrib_pdf, GrouNorm) <- function(distrib, y, theta, log = FALSE) {
  stats::dnorm(y, theta[[1]], theta[[2]], log = log)
}

GrouGamma <- S7::new_class("GrouGamma", parent = continuous_distrib, package = NULL)
S7::method(distrib_pdf, GrouGamma) <- function(distrib, y, theta, log = FALSE) {
  stats::dgamma(y, theta[[1]], theta[[2]], log = log)
}

GrouCauchy <- S7::new_class("GrouCauchy", parent = continuous_distrib, package = NULL)
S7::method(distrib_pdf, GrouCauchy) <- function(distrib, y, theta, log = FALSE) {
  stats::dcauchy(y, theta[[1]], theta[[2]], log = log)
}

GrouBeta <- S7::new_class("GrouBeta", parent = continuous_distrib, package = NULL)
S7::method(distrib_pdf, GrouBeta) <- function(distrib, y, theta, log = FALSE) {
  stats::dbeta(y, theta[[1]], theta[[2]], log = log)
}

GrouExp <- S7::new_class("GrouExp", parent = continuous_distrib, package = NULL)
S7::method(distrib_pdf, GrouExp) <- function(distrib, y, theta, log = FALSE) {
  stats::dexp(y, theta[[1]], log = log)
}
S7::method(distrib_quantile, GrouExp) <- function(distrib, p, theta, lower.tail = TRUE, log.p = FALSE) {
  stats::qexp(p, theta[[1]])
}

bare_instance <- function(cls, nm, bounds, params, pbounds, links) {
  cls(
    distrib_name = nm, dimension = "univariate", bounds = bounds,
    params = params, params_interpretation = stats::setNames(params, params),
    n_params = length(params), params_bounds = pbounds, link_params = links
  )
}

bare_norm <- function() {
  bare_instance(GrouNorm, "grou norm", c(-Inf, Inf), c("mu", "sigma"),
    list(mu = c(-Inf, Inf), sigma = c(0, Inf)), list(mu = id_link, sigma = lg_link))
}

bare_gamma <- function() {
  bare_instance(GrouGamma, "grou gamma", c(0, Inf), c("a", "r"),
    list(a = c(0, Inf), r = c(0, Inf)), list(a = lg_link, r = lg_link))
}

bare_cauchy <- function() {
  bare_instance(GrouCauchy, "grou cauchy", c(-Inf, Inf), c("mu", "sigma"),
    list(mu = c(-Inf, Inf), sigma = c(0, Inf)), list(mu = id_link, sigma = lg_link))
}

bare_beta <- function() {
  bare_instance(GrouBeta, "grou beta", c(0, 1), c("a", "b"),
    list(a = c(0, Inf), b = c(0, Inf)), list(a = lg_link, b = lg_link))
}

test_that("GRoU draws follow the target density on unbounded, half-line and bounded supports", {
  set.seed(101)
  cases <- list(
    gaussian = list(d = bare_norm(), th = list(mu = 1.5, sigma = 2),
      p = function(q) stats::pnorm(q, 1.5, 2)),
    gamma = list(d = bare_gamma(), th = list(a = 4.5, r = 1.5),
      p = function(q) stats::pgamma(q, 4.5, 1.5)),
    beta = list(d = bare_beta(), th = list(a = 2.5, b = 4),
      p = function(q) stats::pbeta(q, 2.5, 4))
  )

  for (nm in names(cases)) {
    y <- rng_grou(cases[[nm]]$d, 20000, cases[[nm]]$th)
    expect_length(y, 20000)
    ks <- suppressWarnings(stats::ks.test(y, cases[[nm]]$p))
    expect_gt(ks$p.value, 0.001)
  }
})

test_that("GRoU copes with heavy tails, where r = 2 keeps the region bounded", {
  set.seed(102)
  y <- rng_grou(bare_cauchy(), 20000, list(mu = 0.5, sigma = 1.4))
  ks <- suppressWarnings(stats::ks.test(y, function(q) stats::pcauchy(q, 0.5, 1.4)))
  expect_gt(ks$p.value, 0.001)
})

test_that("recentring makes GRoU independent of where the distribution sits", {
  # Without the shift to the mode the bounding rectangle degenerates and the
  # acceptance rate collapses to zero.
  set.seed(103)
  y <- rng_grou(bare_norm(), 20000, list(mu = 1000, sigma = 3))
  ks <- suppressWarnings(stats::ks.test(y, function(q) stats::pnorm(q, 1000, 3)))
  expect_gt(ks$p.value, 0.001)
})

test_that("find_pdf_anchor locates the mode of a density far from the origin", {
  # Regression: a golden-section search on the compactified scale was off by
  # hundreds of units here, because dy/dt grows like y^2.
  expect_equal(find_pdf_anchor(bare_norm(), list(mu = 1000, sigma = 3)), 1000,
    tolerance = 1e-6)
  expect_equal(find_pdf_anchor(bare_norm(), list(mu = -1e5, sigma = 10)), -1e5,
    tolerance = 1e-6)
  # gamma(a, r) has its mode at (a - 1) / r
  expect_equal(find_pdf_anchor(bare_gamma(), list(a = 4.5, r = 1.5)), 3.5 / 1.5,
    tolerance = 1e-6)
})

test_that("the default RNG method picks GRoU only when inversion is not analytic", {
  # No quantile method: GRoU, and much faster than one uniroot per draw
  d <- bare_gamma()
  th <- list(a = 4.5, r = 1.5)
  set.seed(104)
  y <- distrib_rng(d, 5000, th)
  expect_length(y, 5000)
  expect_gt(suppressWarnings(stats::ks.test(y, function(q) stats::pgamma(q, 4.5, 1.5)))$p.value, 0.001)

  # With an analytic quantile, inverse transform sampling is used verbatim
  qd <- bare_instance(GrouExp, "grou exp", c(0, Inf), "rate",
    list(rate = c(0, Inf)), list(rate = lg_link))
  set.seed(105)
  a <- distrib_rng(qd, 5, list(rate = 2))
  set.seed(105)
  expect_equal(a, stats::qexp(stats::runif(5), 2))
})

test_that("distributions with a native RNG are unaffected", {
  set.seed(106)
  a <- distrib_rng(gaussian_distrib(), 5, list(mu = 0, sigma = 1))
  set.seed(106)
  expect_equal(a, stats::rnorm(5))
})

test_that("a density diverging at one edge is reparameterised, not refused", {
  # f(y) ~ |y - a|^(alpha - 1) near the edge: sampling |Y - a|^(1/lambda) with
  # lambda * alpha > 1 has a bounded density, and mapping back is exact.
  set.seed(107)
  for (shape in c(0.1, 0.4, 0.85)) {
    y <- rng_grou(bare_gamma(), 20000, list(a = shape, r = 1.5))
    expect_length(y, 20000)
    ks <- suppressWarnings(stats::ks.test(y, function(q) stats::pgamma(q, shape, 1.5)))
    expect_gt(ks$p.value, 0.001)
  }

  # the divergence may sit at the upper edge instead
  set.seed(108)
  y <- rng_grou(bare_beta(), 20000, list(a = 2, b = 0.3))
  ks <- suppressWarnings(stats::ks.test(y, function(q) stats::pbeta(q, 2, 0.3)))
  expect_gt(ks$p.value, 0.001)

  # and at the lower edge of a bounded support
  set.seed(109)
  y <- rng_grou(bare_beta(), 20000, list(a = 0.3, b = 2))
  ks <- suppressWarnings(stats::ks.test(y, function(q) stats::pbeta(q, 0.3, 2)))
  expect_gt(ks$p.value, 0.001)
})

test_that("the reparameterisation does not care where the divergent edge sits", {
  # The power is applied to |Y - edge|, the distance from the singular edge, not
  # to Y itself; since a divergence can only occur at a finite edge that distance
  # is non-negative wherever the edge lies. So the positivity restriction that
  # would apply to a Box-Cox transform of the variable does not arise, and
  # supports on the negative half-line work like any other.
  shifted <- function(nm, pdf, bounds) {
    cls <- S7::new_class(nm, parent = continuous_distrib, package = NULL)
    S7::method(distrib_pdf, cls) <- pdf
    cls(distrib_name = nm, dimension = "univariate", bounds = bounds,
        params = "a", params_interpretation = c(a = "shape"), n_params = 1,
        params_bounds = list(a = c(0, Inf)), link_params = list(a = lg_link))
  }

  # lower edge at -5
  d1 <- shifted("ShiftGamma", function(distrib, y, theta, log = FALSE) {
    stats::dgamma(y + 5, theta[[1]], 1, log = log)
  }, c(-5, Inf))
  set.seed(120)
  y <- rng_grou(d1, 20000, list(a = 0.4))
  expect_true(all(y > -5))
  expect_gt(suppressWarnings(stats::ks.test(y, function(q) stats::pgamma(q + 5, 0.4, 1)))$p.value, 0.001)

  # support entirely negative
  d2 <- shifted("NegBeta", function(distrib, y, theta, log = FALSE) {
    v <- stats::dbeta((y + 10) / 8, theta[[1]], 3, log = TRUE) - log(8)
    if (log) v else exp(v)
  }, c(-10, -2))
  set.seed(121)
  y <- rng_grou(d2, 20000, list(a = 0.35))
  expect_true(all(y > -10 & y < -2))
  expect_gt(suppressWarnings(stats::ks.test(y,
    function(q) stats::pbeta((q + 10) / 8, 0.35, 3)))$p.value, 0.001)

  # upper edge, on the negative half-line
  d3 <- shifted("ReflGamma", function(distrib, y, theta, log = FALSE) {
    stats::dgamma(-3 - y, theta[[1]], 1, log = log)
  }, c(-Inf, -3))
  set.seed(122)
  y <- rng_grou(d3, 20000, list(a = 0.3))
  expect_true(all(y < -3))
  expect_gt(suppressWarnings(stats::ks.test(y,
    function(q) stats::pgamma(-3 - q, 0.3, 1, lower.tail = FALSE)))$p.value, 0.001)
})

test_that("a density diverging at both edges is handled by a two-sided map", {
  # No single power straightens both edges, so a map behaving like a different
  # power at each end is used instead.
  d <- bare_beta()
  for (ab in list(c(0.5, 0.5), c(0.2, 0.3), c(0.05, 0.7), c(0.99, 0.99))) {
    set.seed(110)
    y <- rng_grou(d, 20000, list(a = ab[1], b = ab[2]))
    expect_length(y, 20000)
    expect_true(all(y >= 0 & y <= 1))
    ks <- suppressWarnings(stats::ks.test(y, function(q) stats::pbeta(q, ab[1], ab[2])))
    expect_gt(ks$p.value, 0.001,
      label = sprintf("beta(%.2f, %.2f) KS p", ab[1], ab[2]))
  }
})

test_that("mass that cannot be resolved lands on the edge instead of being lost", {
  # Beta(0.9, 0.1) puts about 2.5% of its mass within one unit in the last place
  # of 1. No sampler working in double precision can place it any more finely
  # than at 1 itself, but it must not be dropped: rejecting those draws would
  # quietly redistribute that mass over the rest of the distribution.
  d <- bare_beta()
  set.seed(113)
  y <- rng_grou(d, 1e5, list(a = 0.9, b = 0.1))
  target <- stats::pbeta(1 - .Machine$double.eps / 2, 0.9, 0.1, lower.tail = FALSE)
  expect_equal(mean(y == 1), target, tolerance = 0.1)
  expect_equal(max(y), 1)
})

test_that("the discrete fallback inverts the step cdf exactly and in one pass", {
  # Nothing is solved numerically here: the cumulative mass function is a step
  # function, so inversion is exact and a single binary search serves the whole
  # sample.
  DiscP <- S7::new_class("DiscPois", parent = discrete_distrib, package = NULL)
  S7::method(distrib_pdf, DiscP) <- function(distrib, y, theta, log = FALSE) {
    stats::dpois(y, theta[[1]], log = log)
  }
  d <- bare_instance(DiscP, "grou pois", c(0, Inf), "mu",
    list(mu = c(0, Inf)), list(mu = lg_link))

  p <- c(0, 0.001, 0.1, 0.5, 0.9, 0.999)
  for (mu in c(0.5, 4, 250)) {
    expect_equal(distrib_quantile(d, p, list(mu = mu)), stats::qpois(p, mu),
      label = paste("mu =", mu))
  }

  # the draws follow the pmf (a continuous KS test would be invalid on ties)
  set.seed(111)
  y <- distrib_rng(d, 5e4, list(mu = 4))
  obs <- as.numeric(table(factor(pmin(y, 15), levels = 0:15)))
  expected <- 5e4 * c(stats::dpois(0:14, 4), stats::ppois(14, 4, lower.tail = FALSE))
  chisq <- sum((obs - expected)^2 / expected)
  expect_gt(stats::pchisq(chisq, df = 15, lower.tail = FALSE), 0.001)

  # per-observation parameters still work
  set.seed(112)
  y <- distrib_rng(d, 4000, list(mu = rep(c(2, 60), each = 2000)))
  expect_equal(mean(y[1:2000]), 2, tolerance = 0.2)
  expect_equal(mean(y[2001:4000]), 60, tolerance = 2)
})

test_that("vector-valued theta is handled group by group", {
  set.seed(108)
  d <- bare_norm()
  y <- distrib_rng(d, 4000, list(mu = rep(c(0, 50), each = 2000), sigma = 1))
  expect_length(y, 4000)
  expect_equal(mean(y[1:2000]), 0, tolerance = 0.15)
  expect_equal(mean(y[2001:4000]), 50, tolerance = 0.15)
})

test_that("rng_grou validates its arguments", {
  expect_error(rng_grou(bare_norm(), 10, list(mu = c(0, 1), sigma = 1)), "scalar")
  expect_length(rng_grou(bare_norm(), 0, list(mu = 0, sigma = 1)), 0)
})
