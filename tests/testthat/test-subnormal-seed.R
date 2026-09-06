# THE MASS AT ZERO AS A SUBNORMAL, in the two families whose expected
# quantities are summed against their own pmf recurrence.
#
# Those series carry P(Y = 0) in log scale until it is representable and then
# switch to the multiplicative form. What decides when to switch is the
# EXPONENT, and a rule reading exp(lpk) instead gets it wrong over a narrow
# band: between the smallest denormal and the smallest normal double the
# exponential returns a number that is not zero -- so such a rule leaves log
# scale -- and carries almost no significand -- so the recurrence it seeds is
# wrong from its first step. At the far edge of the band the mass at zero is
# 4.94e-324, which is one bit.
#
# The band is reached by ordinary models. For negbin1, P(Y = 0) is
# exp(-(mu/theta) log1p(theta)), so at theta = 0.5 it is mu between about 873
# and 919; measured there, mu_theta came back -2.73e-02 where it is
# 1.2089e-04. For negbin2 it is exp(theta log(theta/(theta+mu))), which a
# large mean with a mild overdispersion reaches; measured at mu = 1.702e5,
# theta = 100, E[trigamma(Y + theta)] was 3.7e-02 out.
#
# Every reference below is a sum taken in R over the family's own mass and
# shares no arithmetic with the kernels.

mass_at_zero <- function(family, mu, th) {
  if (family == "negbin1") exp(-(mu / th) * log1p(th))
  else exp(th * (log(th) - log(th + mu)))
}
is_denormal <- function(x) x > 0 & x < .Machine$double.xmin


test_that("negbin1's expected information holds across the band", {
  # read on mu_mu, which is the series over theta^2 and carries no
  # cancellation; mu_theta and theta_theta are differences of terms far
  # larger than themselves, so a tolerance on those would measure the
  # cancellation rather than the series
  rif <- function(mu, th) {
    r <- mu / th
    p <- 1 / (1 + th)
    k <- 0:(ceiling(qnbinom(1 - 1e-14, size = r, prob = p)) + 60)
    w <- dnbinom(k, size = r, prob = p)
    sum((trigamma(k + r) - trigamma(r)) * w) / sum(w) / th^2
  }
  d <- negbin1_distrib()
  th <- 0.5

  # inside the band, and that these points really are inside it
  for (mu in c(880, 900, 918.832)) {
    expect_true(is_denormal(mass_at_zero("negbin1", mu, th)),
                label = sprintf("mu = %g is in the band", mu))
    got <- distrib_expected_hessian(d, 0, list(mu = mu, theta = th))[["mu_mu"]]
    expect_equal(got, rif(mu, th), tolerance = 1e-9,
                 label = sprintf("mu = %g", mu))
  }

  # and either side: an ordinary mass at zero, and one that underflows to
  # exactly zero, where the log scale is never left
  for (mu in c(500, 860, 950, 2000)) {
    got <- distrib_expected_hessian(d, 0, list(mu = mu, theta = th))[["mu_mu"]]
    expect_equal(got, rif(mu, th), tolerance = 1e-9,
                 label = sprintf("mu = %g", mu))
  }
})


test_that("negbin2's expected quantities hold across the band", {
  # both helpers that carry the recurrence are exercised: nb_E_ltt through the
  # expected information, nb_E_psigamma_diff through the expected third and
  # fourth derivatives
  griglia <- function(mu, th)
    0:(ceiling(qnbinom(1 - 1e-14, size = th, mu = mu)) + 60)
  # ⚠️ Read on E[l_theta_theta] AS RETURNED, and not on E[trigamma(Y+theta)]
  # reconstructed from it. Since 0.48.0 the expected hessian returns the two
  # pieces already merged into one sum, so recovering the polygamma
  # expectation means adding trigamma(theta) and taking mu/(theta(theta+mu))
  # off again -- and at mu = 1.702e5, theta = 100 that second term is 2.0e+02
  # times the returned value, so the reconstruction loses two digits of its
  # own. Measured, the kernel agrees with the reference below at 3.96e-09 or
  # better at every point here while the reconstruction is 3.35e-08 out: it
  # was the reconstruction that failed, not the seed this block is about.
  #
  # ⚠️ The tolerance below is 1e-7 where it used to be 1e-9, and that is the
  # PRICE of the merge rather than a widening for convenience. These five
  # points are mu >> theta, which is the one regime where merging costs
  # accuracy: measured against the same reference, the composition it replaces
  # reads 3.45e-12 to 2.28e-11 here and the merged sum 6.98e-10 to 3.96e-09,
  # a factor of 58 to 265. It is paid for in the other direction, where the
  # composition is 2.47e+07 out at mu = 0.1, theta = 1e6 and loses the sign
  # of an expected information; see test-negbin-polygamma-recurrence.R.
  rif_trig <- function(mu, th) {
    k <- griglia(mu, th)
    w <- dnbinom(k, size = th, mu = mu)
    w <- w / sum(w)
    j <- 0:(max(k) - 1)
    sum(w * c(0, cumsum((th * (2 * j - mu) + j^2) /
                          (th * (th + mu) * (th + j)^2))))
  }
  got_trig <- function(mu, th)
    unlist(negbin_expected_hessian_cpp(0, mu, th, 1L))[["theta_theta"]]
  rif_d <- function(fn, mu, th) {
    k <- griglia(mu, th)
    w <- dnbinom(k, size = th, mu = mu)
    w <- w / sum(w)
    o <- fn(k, rep(mu, length(k)), rep(th, length(k)), 1L)
    vapply(o, function(v) sum(v * w), numeric(1))
  }

  th <- 100
  # the band, and one point either side of it
  for (mu in c(1.702e5, 1.54e5, 1.338e5, 1.096e5, 4.024e4)) {
    expect_equal(got_trig(mu, th), rif_trig(mu, th), tolerance = 1e-7,
                 label = sprintf("E[l_theta_theta] at mu = %g", mu))
    expect_equal(unlist(negbin_deriv3_expected_cpp(0, mu, th, 1L)),
                 rif_d(negbin_deriv3_cpp, mu, th), tolerance = 1e-7,
                 label = sprintf("E[d3] at mu = %g", mu))
    expect_equal(unlist(negbin_deriv4_expected_cpp(0, mu, th, 1L)),
                 rif_d(negbin_deriv4_cpp, mu, th), tolerance = 1e-7,
                 label = sprintf("E[d4] at mu = %g", mu))
  }
  # the widest point really is in the band, so the test cannot pass by
  # never reaching it
  expect_true(is_denormal(mass_at_zero("negbin2", 1.702e5, th)))
  expect_false(is_denormal(mass_at_zero("negbin2", 4.024e4, th)))
})
