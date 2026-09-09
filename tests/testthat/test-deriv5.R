# The fifth-order derivative surface.
#
# No family writes the fifth order out. Every one reaches numerical_deriv5(),
# ONE central difference of the analytic fourth -- not a difference of a
# difference, which is the rule the rest of this package's numerical surface
# obeys and which numerical_deriv4() deliberately does not (it differences the
# analytic Hessian twice, which it can afford because a second difference of an
# exact quantity is still only two orders removed from one).
#
# The reference here is Richardson extrapolation on the analytic fourth. That
# is the rule for a quantity which itself contains a difference: a plain
# central difference would be the same arithmetic twice. numDeriv may be named
# in tests, where it is a Suggests dependency; R/ uses numericals7 throughout
# and carries numDeriv_grad() so that numDeriv stays a suggestion.

skip_if_not_installed("numDeriv")


# Richardson on the analytic fourth, component by component. The component of
# a non-decreasing five-tuple is the four-tuple left when its last index is
# dropped, differentiated in the parameter that index names.
ref_deriv5 <- function(d, y, theta, scale = "parameter") {
  params <- d@params
  idx <- deriv_indices(params, 5)
  nms <- deriv_names(params, 5)
  links <- d@link_params
  x0 <- vapply(params, function(nm) {
    if (scale == "link") linkfunctions7::linkfun(links[[nm]], theta[[nm]]) else theta[[nm]]
  }, numeric(1))

  out <- vector("list", length(nms))
  names(out) <- nms
  for (t in seq_along(nms)) {
    m <- idx[[t]][5]
    key4 <- paste(params[idx[[t]][1:4]], collapse = "_")
    f <- function(v) {
      x <- x0
      x[m] <- v
      th2 <- theta
      for (i in seq_along(params)) {
        th2[[params[i]]] <- if (scale == "link") {
          linkfunctions7::linkinv(links[[params[i]]], x[i])
        } else {
          x[i]
        }
      }
      distrib_deriv4(d, y, th2, scale = scale)[[key4]]
    }
    # the component is a vector over observations against one parameter value,
    # so this is a jacobian and not a gradient
    out[[t]] <- as.vector(numDeriv::jacobian(f, x0[m]))
  }
  out
}

# PER COMPONENT against its own scale. A global denominator would hide a defect
# in a small component entirely, and several families have components that are
# identically zero -- a gaussian's log-density is quadratic in mu, so every
# order-5 component carrying three or more mu's vanishes. Such a component has
# no denominator of its own and is judged against the family's scale.
rel5 <- function(a, b) {
  glob <- max(vapply(b, function(v) max(abs(v)), numeric(1)), 1e-300)
  max(vapply(names(b), function(k) {
    den <- max(abs(b[[k]]))
    max(abs(a[[k]] - b[[k]])) / if (den > glob * 1e-12) den else glob
  }, numeric(1)))
}

d5_cases <- function() {
  list(
    gaussian1  = list(d = gaussian1_distrib(),  th = list(mu = 0.3, sigma = 1.2),
                      y = c(-1.3, -0.2, 0.4, 1.1)),
    gamma1     = list(d = gamma1_distrib(),     th = list(mu = 2.0, phi = 0.7),
                      y = c(0.4, 0.9, 1.6, 2.7)),
    poisson    = list(d = poisson_distrib(),    th = list(mu = 2.5),
                      y = c(0, 1, 2, 4)),
    negbin2    = list(d = negbin2_distrib(),    th = list(mu = 3.0, theta = 2.2),
                      y = c(0, 1, 3, 5)),
    student_t1 = list(d = student_t1_distrib(), th = list(mu = 0.2, sigma = 1.1, nu = 6),
                      y = c(-1.3, -0.2, 0.4, 1.1)),
    beta1      = list(d = beta1_distrib(),      th = list(mu = 0.4, phi = 6),
                      y = c(0.15, 0.31, 0.5, 0.68)),
    weibull1   = list(d = weibull1_distrib(),   th = list(mu = 1.5, sigma = 2.0),
                      y = c(0.4, 0.9, 1.6, 2.7))
  )
}


test_that("the fifth order matches Richardson on the analytic fourth", {
  for (nm in names(d5_cases())) {
    cs <- d5_cases()[[nm]]
    for (sc in c("parameter", "link")) {
      got <- distrib_deriv5(cs$d, cs$y, cs$th, scale = sc)
      ref <- ref_deriv5(cs$d, cs$y, cs$th, sc)
      expect_lt(rel5(got, ref), 1e-6,
                label = sprintf("%s deriv5 on the %s scale", nm, sc))
    }
  }
})


test_that("the check can fail", {
  # The negative control. It must go into a component that is NOT identically
  # zero: mu_mu_mu_mu_sigma of a gaussian is, so injecting there would prove
  # nothing at all.
  d <- gaussian1_distrib()
  th <- list(mu = 0.3, sigma = 1.2)
  y <- c(-1.3, -0.2, 0.4, 1.1)

  good <- distrib_deriv5(d, y, th)
  expect_identical(good[["mu_mu_mu_mu_sigma"]], rep(0, length(y)))
  expect_gt(max(abs(good[["mu_mu_sigma_sigma_sigma"]])), 1)

  bad <- good
  bad[["mu_mu_sigma_sigma_sigma"]] <- bad[["mu_mu_sigma_sigma_sigma"]] * 1.05
  ref <- ref_deriv5(d, y, th)
  expect_lt(rel5(good, ref), 1e-6)
  expect_gt(rel5(bad, ref), 1e-3)
})


test_that("the components are the order-5 enumeration and nothing else", {
  y <- c(-1.3, -0.2, 0.4, 1.1)
  g <- distrib_deriv5(gaussian1_distrib(), y, list(mu = 0, sigma = 1))
  expect_identical(names(g), unname(deriv_names(c("mu", "sigma"), 5)))
  expect_length(g, 6L)
  expect_true(all(vapply(g, length, integer(1)) == length(y)))

  s <- distrib_deriv5(student_t1_distrib(), y, list(mu = 0, sigma = 1, nu = 6))
  expect_length(s, 21L)
})


test_that("a higher-accuracy rule agrees, which is what check_distrib asks", {
  # Two rules from numericals7, five nodes against three, sharing only their
  # centre. This is the reference check_distrib() uses, R/ not being allowed
  # to name numDeriv.
  for (nm in names(d5_cases())) {
    cs <- d5_cases()[[nm]]
    a2 <- distrib_deriv5(cs$d, cs$y, cs$th)
    a4 <- numerical_deriv5(cs$d, cs$y, cs$th, accuracy = 4L)
    expect_lt(rel5(a2, a4), 1e-6, label = nm)
  }
})


test_that("the stencil comes from numericals7 and not from a written-out quotient", {
  # The nodes, the weights and the step are the library's, so raising the
  # accuracy is an argument rather than new arithmetic. Reconstructing the
  # default rule from numericals7 by hand must reproduce the method exactly.
  d <- gaussian1_distrib()
  th <- list(mu = 0.3, sigma = 1.2)
  y <- c(-1.3, -0.2, 0.4, 1.1)

  s <- numericals7::fd_offsets(1L, accuracy = 2L)$central
  w <- numericals7::fd_weights(s, 1L)
  h <- numericals7::fd_step(th$sigma, 1L, accuracy = 2L, bounds = c(0, Inf))

  acc <- 0
  for (j in which(w != 0)) {
    th2 <- th
    th2$sigma <- th$sigma + s[j] * h
    acc <- acc + w[j] * distrib_deriv4(d, y, th2)[["sigma_sigma_sigma_sigma"]]
  }
  expect_equal(distrib_deriv5(d, y, th)[["sigma_sigma_sigma_sigma_sigma"]],
               acc / h, tolerance = 1e-14)
})


test_that("has_exact_deriv4 reads the owner, and a wrapper asks its parent", {
  expect_true(has_exact_deriv4(gaussian1_distrib()))
  expect_true(has_exact_deriv4(gamma1_distrib()))
  expect_true(has_exact_deriv4(weibull3_distrib()))       # a reparametrization

  # wrappers of an analytic family are analytic
  expect_true(has_exact_deriv4(fixed(gaussian1_distrib(), mu = 0)))
  expect_true(has_exact_deriv4(zero_inflated(poisson_distrib())))
  expect_true(has_exact_deriv4(zero_adjusted(gamma1_distrib())))
  expect_true(has_exact_deriv4(truncated(gaussian1_distrib(), lower = -2)))
  expect_true(has_exact_deriv4(folded(gaussian1_distrib())))
  expect_true(has_exact_deriv4(transformation(gaussian1_distrib(), exp_transform())))
})


test_that("a density-only family answers, and says the answer is unchecked", {
  # The promise is that a distribution needs only distrib_pdf(), and it is not
  # withdrawn at the fifth order. What is withdrawn is the claim that the
  # answer means anything: there the fourth is itself a difference, so the
  # fifth is a difference of a difference.
  Bare <- S7::new_class("BareGaussD5", parent = continuous_distrib)
  S7::method(distrib_pdf, Bare) <- function(distrib, y, theta, log = FALSE, ...) {
    stats::dnorm(y, theta$mu, theta$sigma, log = log)
  }
  bare <- Bare(distrib_name = "bare", dimension = "univariate",
               bounds = c(-Inf, Inf), params = c("mu", "sigma"), n_params = 2,
               params_bounds = list(mu = c(-Inf, Inf), sigma = c(0, Inf)),
               link_params = list(mu = linkfunctions7::identity_link(),
                                  sigma = linkfunctions7::log_link()),
               params_smooth = c(mu = TRUE, sigma = TRUE),
               params_interpretation = c(mu = "mean", sigma = "standard deviation"))

  expect_false(has_exact_deriv4(bare))
  expect_false(has_exact_deriv4(truncated(bare, lower = -3)))
  expect_false(has_exact_deriv4(fixed(bare, mu = 0)))
  expect_false(has_exact_deriv4(folded(bare)))

  # it answers rather than raising
  got <- distrib_deriv5(bare, c(-1, 0, 1), list(mu = 0.3, sigma = 1.2))
  expect_length(got, 6L)
  expect_true(all(vapply(got, function(v) all(is.finite(v)), logical(1))))
})
