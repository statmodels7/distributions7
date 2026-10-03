# The third and fourth derivatives of the expected information, by one stencil
# on the analytic second derivative, against symbolic derivatives (base R D())
# of the information written out.

d34_symbolic <- function(expr, vars) {
  for (v in vars) expr <- D(expr, v)
  expr
}

# E[l_pp] on the link scale is minus the information I; the reference
# differentiates -I symbolically in the free coordinates.
d34_cases <- function() list(
  poisson = list(d = poisson_distrib(), eta = c(mu = 1.1), y = 3,
                 I = quote(exp(a)), map = c(mu = "a"),
                 theta = function(e) list(mu = exp(e[["mu"]]))),
  gaussian1 = list(d = gaussian1_distrib(), eta = c(mu = 0.3, sigma = -0.4),
                   y = 0.5, I = quote(exp(-2 * b)), map = c(mu = "a", sigma = "b"),
                   theta = function(e) list(mu = e[["mu"]],
                                            sigma = exp(e[["sigma"]]))),
  negbin2 = list(d = negbin2_distrib(), eta = c(mu = 1.1, theta = 0.7), y = 3,
                 I = quote(exp(b) * exp(a) / (exp(b) + exp(a))),
                 map = c(mu = "a", theta = "b"),
                 theta = function(e) list(mu = exp(e[["mu"]]),
                                          theta = exp(e[["theta"]]))),
  beta1 = list(d = beta1_distrib(), eta = c(mu = 0.4, phi = 3.0), y = 0.6,
               I = quote(exp(b)^2 * (trigamma(exp(a) / (1 + exp(a)) * exp(b)) +
                                       trigamma(1 / (1 + exp(a)) * exp(b))) *
                           (exp(a) / (1 + exp(a)) / (1 + exp(a)))^2),
               map = c(mu = "a", phi = "b"),
               theta = function(e) list(mu = plogis(e[["mu"]]),
                                        phi = exp(e[["phi"]]))))

test_that("the names enumerate a pair and a sorted tuple", {
  P <- c("mu", "sigma")
  expect_length(d3expected_names(P), 3 * 4)
  expect_length(d4expected_names(P), 3 * 5)
  expect_identical(d3expected_key(P, 2, 1, 2, 1, 2), "mu_sigma_mu_sigma_sigma")
  expect_true(d3expected_key(P, 1, 1, 2, 1, 2) %in% d3expected_names(P))
  expect_true(d4expected_key(P, 1, 1, 2, 2, 2, 1) %in% d4expected_names(P))
})

test_that("orders three and four on the link scale match the symbolic derivatives", {
  for (nm in names(d34_cases())) {
    cs <- d34_cases()[[nm]]
    P <- names(cs$eta)
    i <- which(P == names(cs$map)[1])
    env <- as.list(stats::setNames(cs$eta, cs$map[P]))
    th <- cs$theta(cs$eta)
    a3 <- distrib_d3expected_hessian(cs$d, cs$y, th, scale = "link")
    a4 <- distrib_d4expected_hessian(cs$d, cs$y, th, scale = "link")
    expect_identical(names(a3), d3expected_names(P))
    expect_identical(names(a4), d4expected_names(P))
    for (tu in deriv_indices(P, 3L)) {
      want <- -eval(d34_symbolic(cs$I, cs$map[P[tu]]), env)
      got <- a3[[d3expected_key(P, i, i, tu[1], tu[2], tu[3])]]
      expect_lt(abs(got - want) / max(1, abs(want)), 1e-8, label = paste(nm, "order 3"))
    }
    for (tu in deriv_indices(P, 4L)) {
      want <- -eval(d34_symbolic(cs$I, cs$map[P[tu]]), env)
      got <- a4[[d4expected_key(P, i, i, tu[1], tu[2], tu[3], tu[4])]]
      expect_lt(abs(got - want) / max(1, abs(want)), 1e-5, label = paste(nm, "order 4"))
    }
  }
})

test_that("on the parameter scale the gaussian's mean entry is -1/sigma^2", {
  d <- gaussian1_distrib()
  s <- 1.7
  a3 <- distrib_d3expected_hessian(d, 0, list(mu = 0.2, sigma = s))
  a4 <- distrib_d4expected_hessian(d, 0, list(mu = 0.2, sigma = s))
  P <- d@params
  # d^k/ds^k (-s^-2) = -(-2)(-3)...(-(k+1)) s^-(k+2)
  expect_equal(a3[[d3expected_key(P, 1, 1, 2, 2, 2)]], 24 / s^5, tolerance = 1e-8)
  expect_equal(a4[[d4expected_key(P, 1, 1, 2, 2, 2, 2)]], -120 / s^6, tolerance = 1e-5)
  expect_equal(a3[[d3expected_key(P, 1, 1, 1, 2, 2)]], 0, tolerance = 1e-8)
})

test_that("a family without an analytic second derivative signals an error", {
  expect_error(distrib_d3expected_hessian(laplace_distrib(), 0,
                                          list(mu = 0, sigma = 1)),
               "no analytic second derivative")
  expect_error(distrib_d4expected_hessian(
    zero_inflated(poisson_distrib()), 1, list(mu = 2, zi = 0.2)),
    "no analytic second derivative")
})

test_that("the kernel's information and its derivative match the generics", {
  for (cs in list(
    list(d = poisson_distrib(), p = "mu", th = list(mu = 2.3), y = 3, e = 0.4),
    list(d = negbin2_distrib(), p = "mu", th = list(mu = 2.3, theta = 1.7),
         y = 3, e = 0.4),
    list(d = beta1_distrib(), p = "mu", th = list(mu = 0.4, phi = 20), y = 0.6,
         e = -0.3),
    list(d = gaussian1_distrib(), p = "sigma", th = list(mu = 0.2, sigma = 1),
         y = 0.5, e = -0.2))) {
    k <- distrib_kernel(cs$d, cs$p)
    th <- cs$th
    th[[cs$p]] <- linkfunctions7::linkinv(cs$d@link_params[[cs$p]], cs$e)
    key <- paste(cs$p, cs$p, sep = "_")
    want <- -distrib_expected_hessian(cs$d, cs$y, th, scale = "link")[[key]]
    expect_equal(k$information(cs$y, cs$th, cs$e), want, tolerance = 1e-12)
    expect_equal(k$dinformation(cs$y, cs$th, cs$e),
                 numDeriv::grad(function(e) k$information(cs$y, cs$th, e), cs$e),
                 tolerance = 1e-8)
  }
})
