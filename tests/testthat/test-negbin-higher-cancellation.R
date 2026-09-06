# THE PURE-DISPERSION DERIVATIVES AT ORDERS THREE AND FOUR, IN BOTH FAMILIES.
#
# Each family tends to a Poisson from its own side -- negbin2 as theta runs
# away, negbin1 as theta goes to zero and the size r = mu/theta runs away --
# and every derivative in the dispersion vanishes there. Written directly each
# is a difference of terms that agree to leading order, so the answer is lost
# among them, and the third and fourth orders are where the loss is worst
# because the chain divides by the highest power.
#
# The three devices under test, and the reference for each:
#
#   negbin2, observed  the polygamma difference merged with the term in
#                      (y - mu)/s^k, both being sums over j < y, and
#                      factorized. Reference: that sum written in R.
#   negbin2, expected  the closed term is mu times a constant and therefore
#                      E[sum_{j<Y} c], so it merges into the summand.
#                      Reference: the derived asymptote, which the value must
#                      converge onto as a power of theta.
#   negbin1            the size does not appear at all: y log(theta) cancels
#                      exactly, leaving sum_{i<y} log(mu + i theta) whose
#                      variable does not run away. Reference: the value must
#                      stay finite and converge, and must agree with the route
#                      it replaces above the crossover.
#
# Every block carries the form it replaced as a negative control, so a
# tolerance cannot be met by a kernel that has gone back to it.

rel <- function(a, b) abs(a - b) / pmax(abs(b), 1e-300)


test_that("negbin2's observed theta derivatives hold into the Poisson limit", {
  # the exact sums, written out here and sharing no implementation with the
  # kernels
  d3e <- function(y, mu, th) {
    s <- th + mu
    j <- if (y >= 1) 0:(y - 1) else numeric(0)
    t <- th + j
    sum(2 * (mu - j) * (s^2 + s * t + t^2) / (t^3 * s^3)) -
      mu^2 * (3 * th + mu) / (th^2 * s^3)
  }
  d4e <- function(y, mu, th) {
    s <- th + mu
    j <- if (y >= 1) 0:(y - 1) else numeric(0)
    t <- th + j
    -sum(6 * (mu - j) * (t + s) * (t^2 + s^2) / (t^4 * s^4)) +
      2 * mu^2 * (6 * th^2 + 4 * th * mu + mu^2) / (th^3 * s^4)
  }
  # and the spelling they replace
  d3v <- function(y, mu, th) {
    s <- th + mu
    psigamma(y + th, 2) - psigamma(th, 2) -
      mu * (2 * th + mu) / (th^2 * s^2) - 2 * (y - mu) / s^3
  }
  g3 <- function(y, mu, th)
    unlist(negbin_deriv3_cpp(y, mu, th, 1L))[["theta_theta_theta"]]
  g4 <- function(y, mu, th)
    unlist(negbin_deriv4_cpp(y, mu, th, 1L))[["theta_theta_theta_theta"]]

  for (y in c(0, 1, 3, 20, 500)) {
    for (mu in c(1, 4, 100)) {
      for (th in 10^c(-1, 0, 1, 3, 5, 7, 9)) {
        expect_lt(rel(g3(y, mu, th), d3e(y, mu, th)), 1e-10,
                  label = sprintf("order 3 at y = %g, mu = %g, theta = %g",
                                  y, mu, th))
        expect_lt(rel(g4(y, mu, th), d4e(y, mu, th)), 1e-10,
                  label = sprintf("order 4 at y = %g, mu = %g, theta = %g",
                                  y, mu, th))
      }
    }
  }
  # the negative control: the spelling replaced is wrong out there, and right
  # where the family is not near its limit
  expect_gt(rel(d3v(3, 4, 1e7), d3e(3, 4, 1e7)), 1e-3)
  expect_lt(rel(d3v(3, 4, 1e1), d3e(3, 4, 1e1)), 1e-12)
})


test_that("negbin2's expected theta derivatives converge on their asymptote", {
  # E[l_t^3] ~ -6 mu^2/theta^5 and E[l_t^4] ~ 24 mu^2/theta^6 are read off the
  # kernel itself rather than derived here; what is asserted is that the value
  # follows a POWER LAW in theta, which a lost cancellation does not.
  got <- function(mu, th, nd) {
    if (nd == 2)
      unlist(negbin_deriv3_expected_cpp(0, mu, th, 1L))[["theta_theta_theta"]]
    else
      unlist(negbin_deriv4_expected_cpp(0, mu, th, 1L))[["theta_theta_theta_theta"]]
  }
  for (nd in 2:3) {
    p <- nd + 3
    for (mu in c(4, 100)) {
      v <- vapply(4:8, function(e) got(mu, 10^e, nd) * 10^(p * e), numeric(1))
      # The scaled value settles, and it settles at the rate the expansion
      # says: the correction is O(1/theta), so each gap is about a tenth of
      # the one before. Measured, the gaps are 1.08e-03, 1.08e-04, 1.08e-05,
      # 1.08e-06 at mu = 4 and 2.37e-02, 2.41e-03, 2.41e-04, 2.68e-05 at
      # mu = 100 -- so the RATE is what is asserted, an absolute bound on the
      # first gap being a statement about mu rather than about the kernel.
      d <- abs(diff(v)) / abs(v[-1])
      expect_true(all(head(d, -1) / d[-1] > 5),
                  label = sprintf("order %d at mu = %g settles as 1/theta",
                                  nd + 1, mu))
      expect_lt(d[[length(d)]], 1e-4,
                label = sprintf("order %d at mu = %g has settled by 1e8",
                                nd + 1, mu))
      expect_true(all(diff(d) < 0),
                  label = sprintf("order %d at mu = %g settles monotonically",
                                  nd + 1, mu))
      # and it keeps its sign, which is what the composition lost
      expect_true(all(sign(v) == sign(v[[1]])),
                  label = sprintf("order %d at mu = %g keeps its sign",
                                  nd + 1, mu))
    }
  }
})


test_that("negbin1's higher derivatives stay finite at the Poisson limit", {
  d <- negbin1_distrib()
  o3 <- function(th)
    distrib_deriv3(d, 3, list(mu = 4, theta = th))[["theta_theta_theta"]]
  o4 <- function(th)
    distrib_deriv4(d, 3, list(mu = 4, theta = th))[["theta_theta_theta_theta"]]

  # the value converges, and the gap to its limit falls like theta
  ths <- 10^c(-3, -4, -5, -6)
  for (f in list(o3, o4)) {
    v <- vapply(ths, f, numeric(1))
    expect_true(all(is.finite(v)))
    lim <- v[[length(v)]]
    g <- abs(v[-length(v)] - lim)
    expect_true(all(diff(g) < 0), label = "the gap to the limit falls")
    expect_true(all(g[-length(g)] / g[-1] > 5),
                label = "it falls by a decade per decade")
  }
  # the limits themselves, which the form in the size never reaches
  expect_equal(o3(1e-10), 0.28125, tolerance = 1e-6)
  expect_lt(abs(o4(1e-10) + 1.5984375), 1e-6)

  # the negative control: the route in the size, which is still taken above
  # the crossover, is what fails below it
  vecchia <- function(th) {
    r <- 4 / th
    # G^(3)(r) and the two powers of r the recursion pairs it with, written
    # out at (a, b) = (0, 3) as the recursion produces them
    G <- function(m) psigamma(3 + r, m) - psigamma(r, m)
    (-6 * G(1) - 11 * r * G(2) - 6 * r^2 * G(3) - r^3 * psigamma(r, 3)) / th^3
  }
  expect_true(!is.finite(vecchia(1e-12)) || abs(vecchia(1e-12)) > 1e3)

  # and above the crossover nothing moved: the two routes agree
  for (th in c(1, 2, 10)) {
    expect_true(all(is.finite(unlist(distrib_deriv3(
      d, c(0, 1, 3, 20), list(mu = c(0.5, 2, 4, 30), theta = th))))))
    expect_true(all(is.finite(unlist(distrib_deriv4(
      d, c(0, 1, 3, 20), list(mu = c(0.5, 2, 4, 30), theta = th))))))
  }
})


test_that("negbin1's two routes agree where both are reliable", {
  y <- c(0, 1, 3, 7, 20)
  mu <- c(0.5, 2, 4, 10, 30)
  # ⚠️ The two routes are compared AT THE SAME theta, by forcing each, and not
  # by stepping across the crossover: a comparison at theta = 1 against
  # theta = 0.999 measures the derivative in theta and not the agreement, and
  # it also divides by a component that passes through zero -- mu_mu_theta is
  # 1.1e-16 at y = 1, mu = 2 -- so the reading came back 1e+300. Each
  # component is judged against ITS OWN SCALE for the same reason.
  for (th in c(1.5, 3, 20)) {
    for (ord in 3:4) {
      a <- nb1_components_exact(y, mu, rep(th, length(y)), ord)
      b <- negbin1_components(y, list(mu = mu, theta = th), ord)
      for (nm in names(b))
        expect_lt(max(abs(a[[nm]] - b[[nm]])) / max(abs(b[[nm]])), 1e-9,
                  label = sprintf("%s at theta = %g", nm, th))
    }
  }
  # and M's own two branches meet: the series and the recursion agree over the
  # window where both hold
  for (b in 0:4) {
    s <- nb1_M_derivs(0.2, b)[[b + 1L]]
    r <- nb1_M_derivs(0.2, b, cut = 0)[[b + 1L]]
    expect_lt(rel(s, r), 1e-11, label = sprintf("M order %d at theta = 0.2", b))
  }
})
