# THE EXPECTED POLYGAMMA DIFFERENCES OF THE TWO NEGATIVE BINOMIALS.
#
# Both families need E[psi^(n)(Y + a)] - psi^(n)(a) over their own support, at
# a shift that is a count. The two terms agree to leading order wherever the
# family tends to a Poisson -- negbin2 as theta runs away, negbin1 as theta
# goes to zero and the size r = mu/theta runs away -- so computing the
# expectation and subtracting afterwards loses the answer. For an INTEGER
# shift the difference has an exact recurrence whose terms carry one sign,
#
#   psi^(n)(x + k) - psi^(n)(x) = (-1)^n n! sum_{j<k} 1/(x + j)^(n+1),
#
# and the kernels accumulate it beside the mass instead.
#
# Every reference below is that recurrence summed in R against dnbinom, which
# shares no arithmetic with the kernels; each test carries the form the kernel
# replaced as a negative control, so the tolerance cannot be met by a route
# that has quietly gone back to subtracting at the end.

supporto <- function(mu, th)
  0:(ceiling(qnbinom(1 - 1e-14, size = th, mu = mu)) + 60)

# E[psi^(n)(Y + theta)] - psi^(n)(theta) for negbin2, at orders 3 and 4
rif_nb2 <- function(mu, th, nd) {
  k <- supporto(mu, th)
  w <- dnbinom(k, size = th, mu = mu)
  w <- w / sum(w)
  cf <- (if (nd %% 2 == 0) 1 else -1) * factorial(nd)
  sum(w * c(0, cumsum(cf / (th + 0:(max(k) - 1))^(nd + 1))))
}

# what the kernels return, carried back to the bare difference. ORDER 1 IS NOT
# HERE: the expected hessian returns E[l_theta_theta] itself, as one sum, and
# reconstructing the difference from it would put back the very cancellation
# that form exists to remove -- see the third block below, which reads it as
# returned. Orders 3 and 4 ADD a closed term back rather than subtracting one,
# which recovers the large quantity and cancels nothing.
got_nb2 <- function(mu, th, nd) {
  s <- th + mu
  if (nd == 2) {
    unlist(negbin_deriv3_expected_cpp(0, mu, th, 1L))[["theta_theta_theta"]] +
      mu * (2 * th + mu) / (th * th * s * s)
  } else {
    unlist(negbin_deriv4_expected_cpp(0, mu, th, 1L))[["theta_theta_theta_theta"]] +
      2 * mu * (th * s - (2 * th + mu)^2) / (th^3 * s^3)
  }
}

# the route this replaced, transcribed: the mass seeded as a difference of two
# logarithms, one trigamma per term, and psi'(theta) taken off at the end
nb2_vecchio <- function(mu, th) {
  ratio <- mu / (th + mu)
  lratio <- log(mu) - log(th + mu)
  kmax <- min(100 + mu + 20 * sqrt(mu * (1 + mu / th)) +
                40 * (mu + th) / th, 2e9)
  s <- 0; cum <- 0
  lpk <- th * (log(th) - log(th + mu))
  logscale <- lpk <= -640
  pk <- exp(lpk); k <- 0
  while (k <= kmax) {
    s <- s + trigamma(k + th) * pk
    cum <- cum + pk
    if (cum >= 1 - 1e-12 && k >= 100) { k <- k + 1; break }
    if (logscale) {
      lpk <- lpk + lratio + log((k + th) / (k + 1))
      pk <- exp(lpk)
      if (lpk > -640) logscale <- FALSE
    } else {
      pk <- pk * (k + th) / (k + 1) * ratio
    }
    k <- k + 1
  }
  if (cum < 1) s <- s + trigamma(k + th) * (1 - cum)
  s - trigamma(th)
}

rel <- function(a, b) abs(a - b) / abs(b)


test_that("negbin2's expected theta derivatives hold into the Poisson limit", {
  # the first three cells are where the family is effectively a Poisson and
  # the two polygammas agree to eleven digits or more
  celle <- list(c(0.1, 501200), c(1, 501200), c(0.1, 158500),
                c(10, 10000), c(5, 100), c(1000, 2), c(0.5, 0.5))
  for (nd in 2:3) {
    for (cc in celle) {
      B <- rif_nb2(cc[[1]], cc[[2]], nd)
      expect_lt(rel(got_nb2(cc[[1]], cc[[2]], nd), B), 1e-10,
                label = sprintf("order %d at mu = %g, theta = %g",
                                nd + 1, cc[[1]], cc[[2]]))
    }
  }
})


test_that("the polygamma route it replaces fails where this one holds", {
  # read on the SUMMAND, E[psi'(Y+theta)] - psi'(theta), which is what the old
  # spelling got wrong; both halves are asserted, since without the second the
  # check would pass for a kernel that had gone back to the old form under a
  # loose tolerance
  rif1 <- function(mu, th) {
    k <- supporto(mu, th)
    w <- dnbinom(k, size = th, mu = mu)
    w <- w / sum(w)
    -sum(w * c(0, cumsum(1 / (th + 0:(max(k) - 1))^2)))
  }
  for (cc in list(c(0.1, 501200), c(1, 501200), c(0.1, 158500))) {
    expect_gt(rel(nb2_vecchio(cc[[1]], cc[[2]]), rif1(cc[[1]], cc[[2]])), 1e-7,
              label = sprintf("the old form at mu = %g, theta = %g",
                              cc[[1]], cc[[2]]))
  }
  # and away from the limit it is fine, so what separates the two is the
  # cancellation and not a difference of formula
  for (cc in list(c(5, 100), c(1000, 2))) {
    expect_lt(rel(nb2_vecchio(cc[[1]], cc[[2]]), rif1(cc[[1]], cc[[2]])), 1e-10,
              label = sprintf("the old form at mu = %g, theta = %g",
                              cc[[1]], cc[[2]]))
  }
})


test_that("the expected information stays positive into the Poisson limit", {
  # E[l_theta_theta] is returned as ONE sum, the two pieces of the composition
  # having merged, because each is of order mu/theta^2 while their sum is of
  # order mu^2/(2 theta^4). The asymptote is DERIVED and not fitted: expanding
  # both pieces in 1/theta, the theta^-3 term contributes mu^2/theta^4 and the
  # theta^-4 term -3mu^2/(2 theta^4).
  info <- function(mu, th)
    -unlist(negbin_expected_hessian_cpp(0, mu, th, 1L))[["theta_theta"]]
  # the composition this replaced, in R
  vecchia <- function(mu, th) {
    k <- supporto(mu, th)
    w <- dnbinom(k, size = th, mu = mu)
    w <- w / sum(w)
    -(-sum(w * c(0, cumsum(1 / (th + 0:(max(k) - 1))^2))) + mu / (th * (th + mu)))
  }

  # an expected information cannot be negative, at any point of the range
  for (mu in c(1, 4, 100)) {
    for (e in 1:8) {
      th <- 10^e
      expect_gt(info(mu, th), 0,
                label = sprintf("the information at mu = %g, theta = 1e%d",
                                mu, e))
    }
  }
  # and it converges onto the derived asymptote, monotonically in theta. The
  # gap is read over 1e4 to 1e7 and NOT further: the asymptote is the leading
  # term alone, so its own next-order truncation is of relative size mu/theta
  # and at mu = 100, theta = 1e8 that is 1e-6, which is what the gap measures
  # there rather than anything about the kernel.
  for (mu in c(4, 100)) {
    d <- vapply(4:7, function(e)
      abs(info(mu, 10^e) - mu^2 / (2 * 10^(4 * e))) / (mu^2 / (2 * 10^(4 * e))),
      numeric(1))
    expect_true(all(diff(d) < 0),
                label = sprintf("the gap to the asymptote falls at mu = %g", mu))
    expect_true(all(head(d, -1) / d[-1] > 3),
                label = sprintf("it falls by a decade per decade at mu = %g", mu))
    expect_lt(d[[length(d)]], 1e-4,
              label = sprintf("the gap at mu = %g, theta = 1e7", mu))
  }
  # the negative control: the composition really does lose the sign there, so
  # the check above cannot be met by a kernel that has gone back to it
  expect_lt(vecchia(4, 1e7), 0)
  expect_lt(vecchia(100, 1e8), 0)
  # and away from the limit the two agree
  for (cc in list(c(4, 1e2), c(100, 1e3), c(5, 10))) {
    expect_lt(rel(info(cc[[1]], cc[[2]]), vecchia(cc[[1]], cc[[2]])), 1e-9,
              label = sprintf("both routes at mu = %g, theta = %g",
                              cc[[1]], cc[[2]]))
  }
})


test_that("negbin1's expected information holds at both ends of its range", {
  # mu_mu is E[psi'(Y + r) - psi'(r)]/theta^2 and carries no cancellation of
  # its own, so a tolerance on it measures the series and not the composition
  rif_nb1 <- function(mu, th) {
    r <- mu / th
    p <- 1 / (1 + th)
    k <- 0:(ceiling(qnbinom(1 - 1e-14, size = r, prob = p)) + 60)
    w <- dnbinom(k, size = r, prob = p)
    w <- w / sum(w)
    -sum(w * c(0, cumsum(1 / (r + 0:(max(k) - 1))^2))) / th^2
  }
  d <- negbin1_distrib()
  # r from 6e5 down to 0.001: the first is the Poisson limit the difference
  # exists for, the last is where psi_T_rest took two trigamma calls a term
  for (cc in list(c(3000, 0.005), c(500, 2), c(5, 0.5), c(3000, 0.5),
                  c(5, 100), c(0.1, 100))) {
    got <- distrib_expected_hessian(
      d, 0, list(mu = cc[[1]], theta = cc[[2]]))[["mu_mu"]]
    expect_lt(rel(got, rif_nb1(cc[[1]], cc[[2]])), 1e-8,
              label = sprintf("mu = %g, theta = %g", cc[[1]], cc[[2]]))
  }
})
