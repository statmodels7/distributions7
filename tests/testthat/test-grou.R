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

test_that("an unbounded density is refused and falls back to inversion", {
  d <- bare_gamma()
  th <- list(a = 0.4, r = 1) # density diverges at zero
  expect_error(rng_grou(d, 10, th), "bounded density")
  set.seed(107)
  expect_warning(y <- distrib_rng(d, 100, th), "inverse transform")
  expect_length(y, 100)
  expect_gt(suppressWarnings(stats::ks.test(y, function(q) stats::pgamma(q, 0.4, 1)))$p.value, 0.001)
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
