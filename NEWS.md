# distributions7 0.39.0

* Every moment method returns one value per parameter setting. 147 of the 168
  forced the shape with `moment_const(theta, k, 0)`, which recycles to the
  length the family's `k` parameters imply; eighteen were written without it
  and answered with the length of whichever components entered the formula, so
  a quantity that does not read the location came back of length 1 when only
  the location varied:

  ```r
  skewness(weibull1_distrib(), list(mu = c(0.1, 1, 100), sigma = 2))
  #> 0.6311107          # length 1, where three settings were asked for
  ```

* Nothing warned. A caller binding moments to the rows of a data frame got one
  number recycled down the column for a Laplace and the right numbers for a
  Gaussian, the two differing only in whether the component being varied
  happens to enter the value.

* The eighteen are the elastic net's `mean`, `variance` and `skewness`;
  `variance` for the Gumbel, both Laplace charts and the pseudo-Huber;
  `kurtosis` for the pseudo-Huber; `mean` for both Poisson-inverse-Gaussian
  charts; `variance`, `skewness` and `kurtosis` for the skew normal in its
  direct chart and for the skew t; and `skewness` and `kurtosis` for the
  Weibull in its scale chart. `kurtosis` for the skew normal's centered chart
  and the Weibull's two in its mean chart follow through their delegations, so
  eighteen edits close 21 (family, generic) pairs.

* No value moved. `moment_const(theta, k, 0)` is a vector of zeros, so every
  one of the 160 moments a scalar theta produces is what it was, the largest
  absolute change over the set being 0.

* Found by an instrument that varies EVERY parameter in turn. A survey that
  varies the first alone reports eighteen pairs and cannot see the three `mean`
  methods that read the first parameter but no later one, which is how the
  figure recorded when the defect was first noticed came to be three short.

* `test-moments-recycling.R` is the regression test, and its sweep cannot be
  satisfied by repairing one family: it asks the question of every parameter of
  every univariate family. Beside it are the closed forms transcribed by hand,
  so a fix that moved a value would fail rather than pass quietly.

# distributions7 0.38.0

* The von Mises distribution function is a series and no longer a quadrature.
  Both parametrizations used the base class's fallback -- one numerical
  integration of the density per observation -- which made them the dearest
  families in the package by three orders of magnitude: measured, a residual
  over a million observations cost 136 seconds against 0.13 to 0.47 for every
  other family, and the whole of it was inside `distrib_cdf`.

* What replaces it is the Fourier expansion of the density integrated term by
  term,
  \deqn{F(x) = \frac{x + \pi}{2\pi} + \frac{1}{\pi I_0(\kappa)}
    \sum_{j\ge1} \frac{I_j(\kappa)}{j}
    \big[\sin(j(x - \mu)) + \sin(j(\pi + \mu))\big],}
  whose second sine is the lower limit of integration and is what makes it
  the distribution function of the family AS WRITTEN, on \eqn{[-\pi, \pi)}
  with the location inside it rather than a variable wrapped around the
  circle. Only the ratios \eqn{I_j/I_0} are needed, and
  `numericals7::bessel_i_ratios()` gives them by a recurrence whose loop runs
  over the series index.

* Measured against the quadrature it replaces, at concentrations from 0.05 to
  200 and three locations: agreement between 1e-15 and 5e-15, which is the
  quadrature's own accuracy, and **61x to 71x faster**. A hundred thousand
  points now take 0.76 seconds where the quadrature took some forty. The
  quantile, which root-finds on the distribution function, and every
  cdf-derivative fallback of the family are faster by the same factor.

* HOW MANY TERMS is measured rather than assumed. Against the same series at
  four times the length, machine precision is reached at 10 terms at
  \eqn{\kappa = 0.5}, 26 at 10, 90 at 100, 242 at 1000 and 404 at 3000 --
  always under \eqn{8.5\sqrt{\kappa} + 10}, which is the rule used and
  which a test pins at four concentrations. The sum is accumulated over
  blocks of observations: the natural expression forms an \eqn{n \times m}
  matrix, which at a hundred thousand points and a concentration of a hundred
  is already hundreds of megabytes.

* Both von Mises families leave the list of those whose distribution function
  is not available in closed form. What remains on it are the six with the
  genuine obstruction, where the derivative of an incomplete gamma or beta in
  its shape is hypergeometric.

# distributions7 0.37.1

* `fit_distrib()` names `crit_grad()` as its stopping rule rather than taking
  the optimizer's default, which optimizers7 0.6.0 widened into a disjunction
  that also stops on a stalled objective. The restart loop reads `converged` as
  the signal to try another start and to fall back to BFGS, so a rule reporting
  success at a stall turns a multi-start search into a single-start one that
  keeps the stall. Measured on `folded(gaussian1_distrib())` at n = 3000 with
  `mu = 1.2` and `sigma = 2`, both runs reporting convergence: under the wider
  rule Fisher scoring stops after 141 iterations at a score of 0.57, with
  `mu = 0.103`, `sigma = 3.572` and a log-likelihood of -5142.97, against 14
  iterations at 2.3e-07, `mu = 1.269`, `sigma = 1.970` and -4729.38 under this
  one. The tolerance is still optimizers7's: the rule is named here, the
  constant is not.

# distributions7 0.37.0

* The Student t's THIRD derivatives survive the `nu` its own chart can
  produce. Every component divided by \eqn{D^3} with
  \eqn{D = \nu\sigma^2 + r^2}, and \eqn{D^3} overflows at
  \eqn{D = 5.6\times10^{102}} while the log link reaches
  \eqn{1.8\times10^{308}}: over `nu` at `sigma` = 1, three of the ten were
  non-finite at 1e150, four at 1e300 and **eight at `double.xmax`**, ten of
  ten on the link scale.

  ⚠️ 0.31.0 made the score and the observed Hessian finite there by the same
  rewrite and did not reach orders three and four -- the shape this file
  records as *when a defect is a shape of mistake, grep for the shape* --
  and orders three and four are exactly what statmodels7's exact outer
  gradient reads, so a fit whose `nu` had gone to its clamp could not be
  certified.

  The substitution is the one that removes \eqn{D}: with \eqn{z = r/\sigma},
  \eqn{u = z^2/\nu}, \eqn{t = 1/(1+u) = \nu\sigma^2/D} and
  \eqn{a = 1 + 1/\nu}, every \eqn{D^{-k}} carries a \eqn{(\nu\sigma^2)^k}
  the numerator supplies, and what is left is a bounded function of
  \eqn{(z, u, t)} over a power of \eqn{\sigma}. **Nine of the ten are exact
  algebra**, with no series and so no crossover to calibrate. The two
  carrying a real cancellation are written out:

  ```
  sigma^3 l_sss = 2nu - 2(1+nu) t^3 (1 - 3u)            -> -2
                = -2 + 2 z^2 a (6 + 3u + u^2) t^3       exactly,
      from  1 - t^3(1 - 3u) = u(6 + 3u + u^2) t^3
  ```

  and the rational part of `l_nununu`, whose bracket
  \eqn{-4(1 + 2t^3 - 3t^2)} factorizes as \eqn{-4(t-1)^2(2t+1)} with
  \eqn{t - 1 = -ut}, leaving
  \eqn{-4z^4t^2(2t+1)/\nu^4 - 8t^3/\nu^3}. Its polygamma pair needs the only
  series here: from the duplication twice differentiated,
  \eqn{\psi''(z+\tfrac12) = 8\psi''(2z) - \psi''(z)}, so

  ```
  psi''((nu+1)/2) - psi''(nu/2) = 8/nu^3 + 12/nu^4 - 20/nu^6 + 84/nu^8
  ```

  the two Bernoulli expansions cancelling term by term at \eqn{\nu^{-2}}.
  `t_T3rest()` carries the \eqn{8/\nu^3} inside, as `t_S()` carries its own
  \eqn{2/\nu^2}, because the only consumer pairs it with a term that cancels
  precisely that.

  Validated four ways, three of them independent of the rewrite:
  against the form it replaces where that still holds (7.7e-16 to 2.2e-15
  up to `nu` = 30); against a central difference of the ANALYTIC Hessian
  over all ten components and `nu` from 2.5 to 1e4 (5.5e-10 to 1.4e-07);
  against the closed limit \eqn{1.5(1 + 2z^2 - z^4)/\nu^4}, which the
  kernel converges onto as \eqn{1/\nu} -- 2.4e-04, 2.4e-05, 2.4e-07,
  2.4e-09, **2.4e-11** at `nu` from 1e5 to 1e12; and by being finite on the
  parameter scale at every `nu` to `double.xmax`.

  ⚠️ My own asymptote was wrong twice before it was right, and the code was
  right both times: the first draft dropped the \eqn{24z^2/\nu^4} that
  \eqn{8\nu^{-3}(1-t^3)} contributes, and the second was a mis-summed hand
  arithmetic. Computing the coefficient rather than writing it out settles
  it in one line.

* ⚠️ **What this does NOT repair: the chain to the LINK scale.** On the
  parameter scale the third derivatives are now finite at every `nu`; on
  the link scale `nu_nu_nu` is still non-finite at 1e150 and
  `mu_nu_nu`, `sigma_nu_nu`, `nu_nu_nu` at 1e300 and above, because
  `to_link_scale()` forms \eqn{(h')^k} with \eqn{h' = \nu} against a
  component of order \eqn{\nu^{-k}}. That is the item already recorded for
  the gaussian's second order at `sigma` = 1.3e154, and closing it means
  restructuring `bell_partial()` to multiply the parameter-scale component
  in first -- the hottest shared path in the package. The FOURTH
  derivatives are untouched here and cede from 1e150 as before.

# distributions7 0.36.0

* The NB1's derivatives in its dispersion no longer cancel, and the header
  of 0.34.0 gains an R twin, `psi_shift_diff()`, for the orders that live
  there.

  This is the same Poisson limit the negative binomial of 0.32.0 reaches,
  from the other side: NB1's size is `r = mu/theta`, so it is `theta -> 0`
  that makes the family tend to the Poisson, and `psi(y+r) - psi(r)` and
  `psi'(y+r) - psi'(r)` both vanish there while the chain rule divides them
  by `theta^2` and `theta^4`.  The amplification is what makes it the worst
  of the set found so far.  Measured against the Poisson limit, at
  `mu` = 4 and `y` = 3, where the score in `theta` is -0.25:

  ```
  theta      1e-6        1e-8            1e-10
  direct   -0.2538     -99.3        +889005
  now      -0.2500      -0.2500        -0.2500
  ```

  a factor of 400 out at 1e-8 and of the wrong sign at 1e-10; the second
  derivative reads 6.1e+03 where it should be of order one, and the fourth
  1.2e+36.

  The two logarithms of the score combine exactly, as the beta-binomial's
  did: with `y/r = y theta/mu`,

  ```
  log1p(y/r) - log1p(theta) = log1p( theta (y - mu) / (mu (1 + theta)) )
  ```

  so none of the leading behavior is formed and then subtracted, and the
  trigamma pair goes through `psi_T_rest()` with its own leading term
  `(y/r)/(y+r)` written out.  The expected information sums the DIFFERENCE
  term by term rather than taking `psi'(r)` off an expectation of
  `psi'(Y+r)`, the two agreeing to leading order.

  `psi_shift_diff(n, k, x)` gives `psi^(n)(x+k) - psi^(n)(x)` at any order
  from one expansion, each power differenced as
  `x^-p (exp(-p log1p(k/x)) - 1)` through `expm1` and `log1p`.  Validated
  against the exact recurrence `(-1)^n n! sum_{j<k} (x+j)^(-n-1)`, which
  holds because the shift is a count and which shares no arithmetic with
  the series: **1.1e-16 to 4.4e-16 at every order from 0 to 3 and every
  argument to 1e12**, where the direct difference is 2.6e-08 out at 1e8 and
  1.1e-03 at 1e12.  At a shift of 500, where the recurrence is dear and the
  series is not, it is exact to the bit.  The expectation is likewise 0 to
  2.1e-16 against that recurrence weighted by the mass, against 1.2e-08 for
  the difference of sums.

  The score, the observed Hessian and their higher orders reproduce the
  forms they replace to 2e-16 wherever those forms still hold, and every
  order agrees with `numDeriv` on the log-mass.

* ⚠️ **Two things this does NOT repair, both measured and both left
  standing rather than described away.**

  The expected information's `theta` block is a cancellation of some
  thirteen digits among three terms of size `mu/theta^2`, which no accuracy
  in the expectation can close: with `E[P_r]` exact to 1.6e-16 the matrix
  is still indefinite from `theta` = 1e-6, its determinant reading -52.6
  there and 6.6e+04 at 1e-7 while the `mu` block stays exact
  (0.2499997 against 0.25).  It is the same shape already recorded for the
  negative binomial and the multivariate t: those expectations cancel one
  order deeper than the score, and closing them needs the composition
  written out symbolically rather than a better summand.

  The score in `theta` itself reaches the Poisson limit to five digits over
  the whole range and no further: its approach is `O(theta)` down to about
  1e-6 and then meets a floor that RISES again, 8e-08, 1.2e-07, 7.6e-06 and
  1.9e-03 at 1e-6 to 1e-12.  The chain is `P (-mu/theta^2) + Q`, two terms
  of size 1e+10 at `theta` = 1e-10 summing to 0.25, so a `P` exact to the
  last bit still leaves 1e-06; closing it means combining the two
  symbolically, as the two logarithms inside `P` already are.  Against a
  factor of 400 and a wrong sign, five digits is the improvement, and the
  test asserts exactly that rather than more.

  The third and fourth derivatives carry a SECOND cancellation, in
  `negbin1_components()`, where terms of size `r^j G^(a+j)(r)` -- of order
  8e+06 at `theta` = 5e-4 -- sum to a value of order one.  `psi_shift_diff()`
  makes each `G` exact and leaves that sum where it was: old and new agree
  to the printed digit down to `theta` = 0.05 and diverge from each other
  below it, with neither trustworthy.  ⚠️ The reference is no help there
  and says so: a central difference of the analytic Hessian at
  `theta` = 5e-4 returns 2.98023224e-01, which is exactly `1e7 * 2^-25`,
  i.e. one ulp of the quantity being differenced.

# distributions7 0.35.0

* The gamma's derivatives in its dispersion no longer cancel, and the four
  quantities the repair needs join the shared `src/psi_diff.h`.

  As the dispersion goes to zero the shape `s = 1/phi` grows and the family
  tends to a normal, so every derivative in `phi` is a polygamma minus its
  own leading asymptote:

  ```
  f1 = log(s) + 1 - psi(s) + log(z) - z      f2 =  1/s   - psi'(s)
  f3 = -1/s^2 - psi''(s)                     f4 =  2/s^3 - psi'''(s)
  ```

  each of which loses its digits as `s` grows.  The score is the clearest
  case: at `y = mu` the data term is exactly zero and the score is
  `-[log(s) - psi(s)]/phi^2`, whose direct form is 1.3e-09 out at
  `phi` = 1e-6, 1.8e-05 at 1e-10, 0.2 per cent at 1e-12 and reads
  **exactly zero** at 1e-14, where the value is -5e+13.  Rewritten it tracks
  the asymptote to between 0 and 2.3e-16 over that whole range.

  ⚠️ Unlike the negative binomial and the beta-binomial, whose limits an
  ordinary fit reaches, this one is prophylactic and is written as such: a
  gamma fit at a dispersion of 1e-4 reaches `s` = 1.0e+04, where the loss is
  1e-11, and reaching 1e+08 would take a coefficient of variation of 1e-4.
  That is a degenerate fit, which is precisely when its derivatives should
  not be noise.

  `f1` splits into two cancelling pairs rather than one, the second being
  `1 + log(z) - z = log1p(w) - w` with `w = z - 1`, which is the `psi_Ew()`
  the negative binomial already carries.  The crossover is 50 for all four,
  measured: the series sits within 5.4e-15, 7.5e-14, 3.7e-13 and 1.2e-12
  there and the direct forms still have every digit.

  The scalar C entry point of the score-driven fast route moved with the
  kernel, the two being held to `identical()` by their twin test.

  Validated against the form it replaces where that still holds (0 to
  7.7e-16 on both components), against `numDeriv` on the log-density at
  every order (2.5e-12 to 4.4e-10), and against a difference of the Hessian
  at the third order (3.7e-10); and the orders above the score come back as
  exact powers of the dispersion, finite and correctly signed, to `s` = 1e15.

* `tests/testthat/test-boundary-cancellation.R` pins all five families
  repaired since 0.31.0 -- the gamma, the negative binomial, the
  beta-binomial, the Student t and the multivariate t -- at the boundary
  each one tends to.  Every case asserts two things, and the second is what
  keeps the first honest: the shipped derivative tracks the limit, and the
  form it replaced does **not**, so reverting to a direct expression fails
  the test rather than passing it silently.

  ⚠️ Writing it corrected the record of 0.34.0.  The measurements quoted
  there are the compiled `betabinom1`'s, which the release did repair -- its
  score converges on the binomial's `(y - n mu)/(mu(1-mu))` where the direct
  chain is 4.0e-08 out at `sigma` = 1e-8 and 3.6e-04 at 1e-12.  The
  parametrization by the shapes, `betabinom2`, computes the same differences
  in R and was **not** touched; it cedes later, by 2.3 per cent at a
  concentration of 1e+14, and is recorded here rather than left to be found.

# distributions7 0.34.0

* The beta-binomial's shape derivatives no longer cancel, and the three
  quantities the repair needs are now shared with the negative binomial in
  `src/psi_diff.h` rather than written twice.

  As the concentration `S = A + B` grows the family tends to the binomial and
  every derivative in the shapes vanishes, so each was a difference of
  digammas at arguments a whole SIZE apart.  Measured, `dl/dA` is wrong by
  5.6e-05 at `S` = 1e6, **exactly zero** at 1e9 where the value is 1.9e-17,
  and 1.9e+08 out at 1e12 -- while a fit reaches there without trying: at a
  true concentration of 3000 it reports 1.7e+08, and on binomial data
  3.1e+09.

  ⚠️ The log-mass has carried the exact form since 0.20.0, summing
  `log(A+j)` over the support rather than differencing two `lbeta`; these are
  the derivatives OF THAT SUM and had not followed it.  A repair applied
  where a defect was found and not to the quantities derived from it is the
  shape this file records elsewhere as *when a defect is a shape of mistake,
  grep for the shape*.

  With `psi(x+k) - psi(x) = psi_A_rest(k,x) + log1p(k/x)` the two logarithms
  combine into a single `log1p` of a small quantity,
  `log1p(y/A) - log1p(n/S) = log1p((y S - A n)/(A(S+n)))`, and the trigamma
  pair the same way through `psi_T_rest()`.  Nothing added is `O(n)`: the
  series are `O(1)` in the size, where the log-mass's own sums are not.

  Validated against `numDeriv` on the log-mass, which shares no arithmetic
  (4.2e-10 to 5.8e-09 on the gradient, 1.3e-12 to 3.5e-11 on the Hessian);
  by the score keeping its sign and decaying as `1/S` -- -5.4993e-02,
  -5.5714e-05, -5.5713e-08, exactly a thousandfold per decade; and by the
  negative binomial being unmoved to 3e-16 across the move onto the shared
  header.

  ⚠️ A residual reported rather than explained: that `1/S` law holds cleanly
  to `S` = 2e9 and then departs, by 1.7 per cent at 2e12 and a factor of 2.6
  at 2e15.  The values stay finite and correctly signed, and the departure
  has not been chased.

* `src/psi_diff.h` carries `psi_A_rest()`, `psi_T_rest()` and `psi_Ew()`,
  each taking the SHIFT -- a count, a size or a dimension, and therefore an
  integer -- and returning the difference with its own leading behavior
  subtracted, so a caller pairing it with that behavior cancels
  symbolically.  The header states which families reach which boundary and
  what each direct form was measured to cost.

# distributions7 0.33.0

* The multivariate Student t's score and observed Hessian in `nu` no longer
  cancel.  The structure is the univariate family's with the dimension `p` in
  place of 1, and it fails the same way: as `nu` grows the family tends to the
  multivariate gaussian, every derivative in `nu` vanishes, and each was
  written as a difference of terms agreeing to leading order --
  `psi((nu+p)/2) - psi(nu/2)` is `p/nu` and so is what it is subtracted from.
  Measured on the score the family returns, the direct form is wrong by
  4.9e-05 at `nu` = 1e6, by 0.397 at 1e8 and by 838 at 1e10, and it CHANGES
  SIGN: at `nu` = 1e9 it read +8.2e-15 where the trend of the values below it
  gives -2.2e-16.

  **The repair needs no series**, because `p` is an integer dimension, so the
  shift between the two arguments is a whole number of steps of the recurrence
  `psi(x+1) = psi(x) + 1/x`.  For even `p` that gives sums whose terms all
  carry ONE SIGN, and nothing cancels:

  ```
  A_p(nu)                                   = -sum_{j<p/2} 4j / (nu (nu+2j))
  [psi'((nu+p)/2) - psi'(nu/2)]/2 + p/nu^2  = sum_{j<p/2} 8j(nu+j) / (nu^2 (nu+2j)^2)
  ```

  Both are exactly zero at `p` = 2, which is what each identity gives there
  and what the direct forms return as noise at 1e-16.  For odd `p` the shift
  is a half-integer and the recurrence carries the quantity onto the
  univariate `A_1(nu)`, which keeps a series above a measured crossover -- the
  same expansion `student_t.cpp` carries, and the one place in the package
  where that series exists twice.  The remaining pair of the score,
  `(nu+p) q/(nu(nu+q))` and `log1p(q/nu)`, is `D(u) + (p/nu) u/(1+u)` with
  `u = q/nu` and the same `D(u) = u/(1+u) - log1p(u)`.

  Validated three ways: the helpers against the direct forms where those still
  have their digits (1e-16 to 6e-12 over `p` = 2 to 6 and `nu` = 2.5 to 400,
  even and odd alike); the package's own gradient and Hessian against
  `numDeriv` on the log-density, which shares no arithmetic (4.4e-10 to
  2.6e-07, which is numDeriv's own accuracy on a fifteen-parameter function);
  and the score in `nu` decaying cleanly and keeping its sign to the edge of
  the chart -- -6.5e-08, -8.29e-11, -8.37e-15, -8.37e-19, -8.37e-25,
  -8.37e-201 and exactly 0 at 1.79e308.

* ⚠️ The EXPECTED information in `nu` is not repaired here, for the same
  structural reason as the negative binomial's in 0.32.0: it cancels one order
  DEEPER than the other two, its three terms agreeing at `nu^-2` and again at
  `nu^-3`, so it is a derivation of its own rather than a transcription.

# distributions7 0.32.0

* The negative binomial's score and observed Hessian in the dispersion no
  longer cancel.  As `theta` grows the family tends to the Poisson and every
  derivative in it vanishes, so each was written as a sum of terms that cancel
  to leading order -- and the score's four terms cancel PAIRWISE:
  `psi(y+th) - psi(th)` is `y/th`, `log(th/(th+mu))` is `-mu/th` and
  `(mu-y)/(th+mu)` is `(mu-y)/th`, and the three sum to zero.  The
  cancellation is therefore of order `theta`, not of order `theta/y` as the
  digamma difference alone suggests: measured, the direct form is wrong by
  **1.0e-03 at theta = 1e6, by 4.4 at 1e7, and it CHANGES SIGN at 1e8**.

  And a fit reaches there routinely.  On 2000 counts with `mu = 4` drawn at a
  true `theta` of 100, `fit_distrib()` reports **1.6e+07**; on Poisson counts
  it reports 2.3e+05.  Where such a fit stops is therefore decided by which
  wrong value happens to cross the tolerance rather than by the likelihood.

  Each cancellation is now performed symbolically.  With `a = theta`,
  `b = theta + y` and `c = theta + mu`:

  ```
  dl/dtheta   = [psi(b) - psi(a) - log1p(y/a)] + [log1p(w) - w],  w = (y-mu)/c
  d2l/dtheta2 = (y-mu)^2/(b c^2) + [psi'(b) - psi'(a) + y/(a b)]
  ```

  The Hessian's first quotient is an EXACT identity --
  `-y/(ab) + mu/(ac) + (y-mu)/c^2 = (y-mu)^2/(b c^2)` -- so its three leading
  terms need no series and no crossover at all; only the trigamma remainder
  does.  The score's two brackets are `y/(2ab) + ...` and `-w^2/2 + ...`, each
  with a series below its own measured crossover.

  The two derivations check each other: to leading order the score is
  `[y - (y-mu)^2]/(2 th^2)` and the Hessian `[(y-mu)^2 - y]/th^3`, which is
  its derivative.

  Validated four ways: every one of the five components agrees with the form
  it replaces where that form still has its digits (0 to 1.2e-15 at
  `theta` = 0.5, 3 and 30); `numDeriv` on the log-mass, which shares no
  arithmetic, gives 1.2e-09 on the gradient and 4.0e-11 on the Hessian; the
  score now tracks its Poisson-limit asymptote to 1.3e-04 at `theta` = 1e5
  and exactly at 1e9 and 1e12, where before it was 4.4x wrong at 1e7; and
  every component is finite at `theta` of 1e15, 1e100 and 1.79e308, the value
  the log link clamps to -- `(2a+y)/(ab)` is written `2/b + y/(ab)` so the
  product `a*b`, which overflows past 1.3e154, is never formed.

* ⚠️ **The EXPECTED information in `theta` is NOT repaired**, and it is what
  `iwls()` reads by default.  Measured, at `theta` = 1e6 it returns
  **-1.7e-16**, and an expected information cannot be negative.  Its leading
  order needs one term MORE of the observed Hessian than the rewrite above
  carries, the `theta^-3` term vanishing under expectation, so it is a
  derivation of its own rather than a transcription of these two.

* ⚠️ And `stats::dnbinom` itself loses the Poisson limit: at `y` = 2 and
  `theta` = 1e10 the difference from `dpois` reads **-4.1e-08** where the
  value is +1.0e-10.  It is not repaired and the reason is that the
  consequence is far smaller than for the score: the log-mass VALUE agrees
  with `dpois` to twelve digits, and what is corrupted is only its difference
  from the limit, which nothing in the package computes -- where the score IS
  that difference.  Repairing it would mean writing a log-mass of our own in
  place of R's, which is a different decision from rewriting a derivative.

# distributions7 0.31.0

* The Student t's derivatives in `nu` no longer cancel, and nothing in the
  family overflows at the `nu` its own chart can produce.  Two defects, both
  reachable from an ordinary fit and both silent.

  **The cancellation.**  The score, the observed `nu_nu` and the expected
  `nu_nu` were written as differences of digamma or trigamma at arguments
  half a unit apart, which agree to leading order: the expected information
  **lost its sign from `nu` = 3.2e5**, reading +2.2e-23 where the value is
  -3.5e-24, and on the link scale the factor `nu^2` then made it read
  **-500 = -n/2**.  A scoring step on a negative information walks the fit
  the wrong way.  Each now takes an asymptotic branch above a measured
  crossover, from the duplication `psi(2z) = [psi(z) + psi(z+1/2)]/2 + log 2`
  and its derivative:

  ```
  A(nu)      = psi((nu+1)/2) - psi(nu/2) - 1/nu
             = 1/(2 nu^2) - 1/(4 nu^4) + 1/(2 nu^6) - ...        (nu >= 200)
  S(nu)      = psi'((nu+1)/2) - psi'(nu/2) + 2/nu^2
             = -2/nu^3 + 2/nu^5 - 6/nu^7 - ...                   (nu >= 100)
  E[l_nu_nu] = -7/(2 nu^4) + 13/nu^5 - 79/(2 nu^6) + 119/nu^7    (nu >= 1000)
  ```

  and the score's two remaining terms, `((nu+1)res^2)/(nu den)` and
  `log(1 + res^2/(nu sigma^2))`, which are both `u + O(u^2)`, are replaced by
  the one function they leave behind, `D(u) = u/(1+u) - log1p(u)`, with a
  series of its own below `u` = 1e-3.  The crossovers are where the two
  routes agree best, measured; the series reproduce the direct forms to six
  significant digits at `nu` = 1e3.

  **The overflow.**  The kernels formed `nu sigma^2`, which is `Inf` at the
  `nu` the log link clamps to (1.8e308), and `(nu+1) res / den` is then
  `Inf/Inf`: **the whole score came back `NaN`** on a fit that had
  legitimately run `nu` towards its boundary.  Every expression that divided
  a numerator growing with `nu` by a denominator growing with `nu` is written
  in the ratio instead, with `z^2 = res^2/sigma^2` formed BEFORE dividing by
  `nu` so the product is never taken.  All fifteen components are finite at
  `nu` = 1e6, 1e100 and 1.79e308.

  Validated against four references: the direct forms where those still have
  their digits (1e-16 on every component at `nu` = 2.5, 5, 30, and 7e-14 on
  the two deliberately improved ones at 300), the sign of the expected
  information over 24 values of `nu` from 3.16 to 1e12 (no violations against
  three before), `numDeriv` on the log-density (3.9e-10 on the gradient and
  6.5e-11 on the Hessian, which is numDeriv's own accuracy), and the gaussian
  limit, which the log-density approaches as 1/nu.

# distributions7 0.30.0

* The generalized Pareto's third and fourth derivatives cost **17 ms at
  n = 20000 where they cost 660**, which was the largest single number in
  the derivative census and 97.9 per cent of it sat in one R function.
  Two things were wrong with it and neither was parallelism. Its two
  branches were each evaluated over the WHOLE vector and subset afterwards,
  so a sample straddling the cut paid for both in full; and the near-zero
  branch raised two elementwise powers per term, which are algebraically
  one -- with `u = xi z`, `xi^(k-b) z^(k+1) = u^(k-b) z^(b+1)`, so
  `z^(b+1)/sigma^a` leaves the loop and what remains is a POLYNOMIAL IN u
  with scalar coefficients. That is a scalar recursion of forty-one steps
  an element, so it is compiled (`gpd_poly_cpp`, Horner from the highest
  power down, which sums a decaying series smallest-first where the loop
  summed largest-first). Checked against the previous implementation over
  266 components spanning both branches and every order: 4.4e-16 relative
  to each component's own scale.

* **153 of the 159 exported kernels take a `threads` count**, against 60
  before. The six without are the two jet twins, which exist only as the
  tests' independent reference, and the pseudo-Huber's four (below), so
  every kernel is either threaded or refused for a stated reason.
  Measured at n = 40000 on eight threads: chisq 5.8x on the score and 6.7x
  at fourth order, gengamma 6.4x, gpd 4.4x, lognormal 4.2x, weibull 4.4x,
  logistic and gumbel 4.0x, skewnormal 4.2x, cauchy 2.6x, the binomial pair
  2.1x to 2.5x, and every one of them identical at any count.

* The Poisson-inverse Gaussian's kernels are among them, 3.9x and 4.8x on
  eight threads, and the explicit route still agrees with the mechanical
  jet transcription the tests keep beside it to 1.6e-14. Their jet twins
  take no count, existing only as that reference.

* `pseudohuber_distrib()` is deliberately NOT converted: its kernels call
  R's `bessel_k`, which can raise a warning, and a warning from a worker
  thread ends the session. It is the same refusal `numericals7`'s
  `log_bessel_k()` carries, and for the same reason.

* ⚠️ The conversion surfaced the trap `d7_par.h` warns of, in three
  families at once: `chisq`, `exponential` and `geometric` hoisted the
  parameter out of the loop and wrote it inside (`if (!mu_is_scalar) m =
  mu[i];`), which is shared state once the iterations are split. It shows
  ONLY where the parameter varies by observation -- with a scalar every
  thread writes the same value and the answer comes out right by accident
  -- so the twin test added with them uses a parameter per observation.

* `kMinTiny`, a fourth cost class for bodies of about four nanoseconds an
  observation. The geometric's score is two divisions and does not break
  even until about 150000: measured 0.79x at 40000, 1.29x at 200000, 1.58x
  at 1000000, where `kMinCheap` was measured for bodies twice as dear.

# distributions7 0.29.0

* Every `distrib_pdf()` method takes `...`. The generic is
  `function(distrib, y, theta, ...)` and 43 of 45 methods did not absorb
  what it may be handed, so any caller passing an argument the family does
  not read broke it -- which is what happened the moment the fitting layer
  began passing a thread count. The derivative generics' methods had carried
  `...` all along; this brings the density surface into line with them.

* `vonmises1_distrib()` and `vonmises2_distrib()` carry the count down to
  `numericals7::log_bessel_i()`, which is where this family spends its time:
  profiled at **80.8 per cent** of a fit whose concentration is modelled, the
  concentration then being a vector rather than one number. Measured end to
  end at n = 8000 with both parameters smoothed, 5.7 s against 2.0, and the
  coefficients and the log-likelihood are identical to the bit.

# distributions7 0.28.0

* A parallel kernel is reproducible again, and the cross-count twins ask
  for `identical()` rather than the tolerance of 1e-13 that 0.27.4 settled
  on. That release read the last-bit differences out of R's polygamma path
  as the runtime's and unbindable; re-measured, the reading was wrong on
  both counts. They are not deterministic -- `gamma1`'s third derivative at
  `phi = 1/19` returned six distinct results over six identical calls at
  one thread count, and `negbin2`'s returned five to ten -- so a fit's
  answer moved between two runs whenever a shape landed near one of the
  arguments where `psigamma` diverges (measured at x = 19 and x = 40, 1.3
  and 0.8 ulp; `bessel_k`, and `pgamma` and `pbeta` on the log scale,
  behave the same way). And they are bindable: the worker of
  `d7::par_for()` now installs the calling thread's floating-point
  environment before running its chunk, which makes the parallel branch
  reproduce the sequential value exactly at every argument probed. It costs
  one call per chunk and the measured gains are unchanged (`gamma1` 6.7x ->
  7.3x at eight threads, `beta1`'s third derivative 5.9x -> 7.2x).

* `threads` says how many, not merely whether. `d7::par_for()` passes the
  count to `parallelFor()`, whose `resolveValue()` prefers an explicit
  positive value to `RCPP_PARALLEL_NUM_THREADS`, so a fit that sized the
  pool through `numericals7::local_threads()` is unaffected. Every other
  caller was running on all of the machine's cores whatever it asked for:
  measured on 24 cores, `threads = 2` gave 13.9x, the same as
  `threads = 24`, and now gives 2.03x against 4.00x at four and 7.79x at
  eight.

* The comment in `d7_par.h` states which Rmath routines a body may call as
  a measured list rather than as a family name. `digamma` and `trigamma`
  are thread-stable and `psigamma` at higher orders is not, so "the digamma
  family" was never the right unit; and the routines that can raise a
  warning -- the p/q family, `lchoose`, the Bessel functions -- are
  excluded for a different reason, a warning from a worker thread killing
  the process. `betabinom.cpp` records that its `lchoose` calls are
  admissible only because the support guard keeps a non-integer argument
  from ever reaching them.

# distributions7 0.27.4

* The cross-count twins of the parallel kernels compare at a tolerance of
  1e-13 instead of `identical()`. The Windows CI runner's R runtime
  returns per-thread last bits from its own polygamma path: one ulp,
  deterministic, the same value at the same index across three independent
  binaries -- with the sequential branch routed through the worker's own
  function, and again with that function noinline -- which is the opposite
  of a race signature and not something package code can bind. The
  decomposition guarantee (no reduction is ever split) stands, and the
  tolerance still fails a split reduction or a data race by ten orders.

# distributions7 0.27.3

* The worker's loop in `d7::par_for()` is marked noinline, so the
  sequential branch and the parallel one execute the single compiled copy
  they both call. Routing the sequential branch through the same source
  function (0.27.2) had not been enough: the compiler inlined it at each
  call site and optimized the two copies apart, and the Windows CI runner
  went on producing one-ulp differences between one and two threads in the
  negbin kernels. Bit-identity across counts has to be a property of the
  binary, not of the source.

# distributions7 0.27.2

* `d7::par_for()`'s sequential branch runs through the worker over the whole
  range instead of writing a loop of its own, so both branches execute the
  same compiled function: the bit-identity across thread counts becomes a
  property of the code rather than of the optimizer. The Windows CI runner
  had reported last-bit differences between one and two threads in the two
  negbin kernels whose per-element arithmetic wraps `R::psigamma` calls,
  which no reordering of the decomposition can explain and which does not
  reproduce on this machine's compiler.

# distributions7 0.27.1

* The negative binomial's expected series helpers (`nb_E_trigamma`,
  `nb_E_psigamma`) no longer call `qnbinom` or `dnbinom` inside the
  parallel bodies: the quantile's search reaches `pbeta`, whose warning
  path calls into the R API, and a warning raised from a worker thread
  trips R's C-stack check and killed the test process on four of the five
  CI platforms. The series now stops on its own accumulated mass at the
  point the quantile located, with the underflowing head carried in log
  scale; the switch back to the multiplicative recurrence waits until the
  mass is comfortably normal, since seeding it at a subnormal was measured
  to carry a 2.5x error to the mode at `mu = theta = 1e4`. Values agree
  with the previous sizing to 1e-10 or better across nine regimes, and the
  rule in `d7_par.h` now names Rmath's p/q functions as off limits inside
  a worker.

# distributions7 0.27.0

* The per-observation derivative kernels of the transcendental compiled
  families -- poisson, negbin2, negbin1, beta1, student_t1, invgauss1,
  betabinom1 and gamma2, forty-five kernels over fourteen files -- run
  their loops through the same d7 driver gaussian1 and gamma1 already use,
  at the transcendental threshold (128), with `threads` on their
  derivative methods. The loop bodies are untouched, so the results are
  bit-identical at any count, which the suite asserts kernel by kernel
  with `identical()`; every file was read before conversion for the
  shared-buffer shape that would have made a mechanical pass a data race
  (none carried one).

# distributions7 0.26.0

* Scalar C entry points for the fast route of a score-driven filter
  (piano_parallel.txt, section 2a), registered with `R_RegisterCCallable`:
  `d7_scalar_id` keyed by the family's S7 class name (gaussian1 and
  gamma1; an unknown name answers -1 and the consumer keeps its R
  callbacks) and `d7_score_curv`, the score and the (k, k) second
  derivative of the log-density in one parameter at one observation on the
  parameter scale, mirroring the family's own vector kernels expression by
  expression. A twin test holds them against `distrib_gradient()` and
  `distrib_hessian()` with `identical()`. The remaining compiled families
  take the same few lines each when a measurement names them; the sixteen
  in vectorized R have no C body to point at and stay on the callbacks.

# distributions7 0.25.0

* `fit_distrib(threads = numericals7::n_threads())` accepts the toolkit's
  thread policy. The count travels as an argument down to the family's
  compiled per-observation kernels; the process-level RcppParallel setting
  is sized at the fit's entry and restored on exit. At the default,
  `n_threads(1)`, the code takes exactly the sequential path.
* The per-observation derivative kernels of `gaussian1_distrib()` (all four
  orders) and `gamma1_distrib()` (all four orders, observed and expected)
  run their loops through one RcppParallel driver, decomposed over the
  elements of the output: observation i's derivatives are computed and
  written in full by one thread, so no reduction is split and the result is
  bit-identical at any thread count, which a test asserts with
  `identical()`. Their derivative methods take a `threads` argument
  (default 1) through the generics' dots. Below a measured internal
  threshold a kernel stays sequential whatever the count says; the
  remaining compiled families take the same one-line conversion when a
  measurement names them.

# distributions7 0.24.0

* `distrib_dexpected_hessian()` is the derivative of the expected information
  in the parameters, one component per pair `(a, b)` and differentiating
  parameter `c`. It exists for a marginal criterion whose penalized matrix
  carries the EXPECTED information: that matrix enters through its
  determinant, so its gradient asks for `dK/dbeta`, which is `-l'''` with the
  observed information and `-dE[l'']/deta` with the expected one -- and the
  two are different objects, because differentiating an expectation moves the
  measure as well as the integrand. The missing piece,
  `E[l_ab l_c]`, is a mixed moment no Bartlett identity isolates: the third
  ties the SYMMETRIZED sum, not the single term. The components are symmetric
  in `(a, b)` and NOT in `c`, so they are keyed by `dexpected_names()` rather
  than by the sorted triples `deriv_names()` uses at order three.

* The default method is ONE central difference of the family's own expected
  information, which is a single stencil on an analytic quantity wherever that
  information is a written-out formula -- the licence the skew t already has
  for its degrees of freedom, and not the nested differencing the package
  forbids. Validated against the gaussian's hand-written components (7.2e-11)
  and against the identity `E[l_abc] + E[l_ab l_c]` computed by quadrature on
  a beta, every one of whose six components is non-zero (8.5e-10).

* It REFUSES where the expected information is itself approximated, and the
  reason is cost rather than accuracy: measured at 100 observations the six
  families that approximate it cost 1880 to 147300 ms against a median of
  0.183 ms for the thirty-four that do not, so 2p of those calls per
  evaluation is not a slower route but an unusable one.

* `has_exact_expected_hessian()` follows the arithmetic instead of the owning
  class, through the new generic `expected_hessian_exact()`. Reading the owner
  is not sufficient and two families prove it: `pseudohuber_distrib()`
  registers a method that calls the numerical `expected_derivative()` and then
  replaces the two components vanishing by symmetry, and
  `skewnormal2_distrib()` registers the chain onto `skewnormal1_distrib()`,
  whose expected information is the base class's quadrature -- costing 5220 ms
  where the parent it chains onto costs 2230. Both answered "written out"
  about a quadrature, and the consequences were live: `fit_distrib()` rejected
  a legitimate `fisher_scoring(approx = )` on them with a message stating the
  family computes its expected information in closed form, and its
  standard-error branch entered a multi-second quadrature believing it a
  formula.

# distributions7 0.23.0

* The centered skew normal rejects its derivatives at ZERO SKEWNESS, naming
  the reason and the parametrization that has none. Its map to the direct
  parameters runs through the cube root of `gamma1`, whose derivative is
  unbounded there: the first derivatives of the log-density survive the limit
  -- they approach a finite value from both sides, the map's factor
  cancelling -- and the second ones grow like `gamma1^(-2/3)`, which is a
  property of the CENTERED parametrization and not of the family. Until now
  the resulting `NA` reached a comparison several frames further on and the
  message named nothing. Patching the first order alone would have been
  worse: the generic that a marginal criterion reads is the second.
  ⚠️ The point matters because it is where a hyperparameter STARTS: the
  bounds on `gamma1` are symmetric, so the midpoint rule puts the starting
  value at exactly the one point the family has no derivatives at.

* `distrib_cross2_y()`, `distrib_grad_y_hess()` and `distrib_hess_y_hess()`
  are closed form for the multivariate Student t, which was the last family a
  marginal criterion could not estimate a correlated prior's matrix with.
  Unlike the gaussian's, this family's response Hessian DEPENDS on the
  observation, so the first and third return one matrix per row.

* All four orders are written on ONE set of pieces (`mvt_dpieces()`): the
  response derivatives are `-c w` and `-c Sigma^-1 + 2d ww'`, so everything
  follows from the first and second derivatives of `s = nu + q`, of `w`, of
  `Sigma^-1` and of the scalars `c` and `d` that follow from `s`. The second
  derivatives all vanish except four, and three of those share one middle
  matrix `A_k Sigma^-1 A_l + A_l Sigma^-1 A_k - A_kl`. `distrib_cross_y()` was
  rewritten to read the same pieces, so its existing checks validate them at
  first order -- the licence the toolkit uses for an order it cannot check
  directly.

* Verified against ONE difference of the analytic quantity below each, at
  p = 2 and 3: 1e-10 to 2e-10 throughout. And the gaussian limit is reached AT
  THE RATE 1/nu -- a factor of 1e4 in nu divides every gap by 1e4, in all
  three derivatives and at both dimensions.

# distributions7 0.22.0

* `distrib_cross_y()` is closed form for the multivariate Student t as well.
  The response gradient is `-c w`, so every component carries the derivative
  of the weight beside the gaussian term it multiplies, and the degrees of
  freedom contribute `-(q - p) w / (nu + q)^2`. Nothing here is obstructed:
  the log-density carries no distribution function, only `lgamma`, a
  logarithm and a quadratic form, each elementary in `nu`. Against numDeriv on
  the analytic response gradient, 1e-9 to 1e-10 at p = 2, 3, 4; and the whole
  block becomes the gaussian's AT THE RATE 1/nu -- 1.17e-03, 1.17e-05,
  1.17e-07 at nu of 1e4, 1e6, 1e8 -- which an arithmetic accident does not do.

* `distrib_cross2_y()`, `distrib_grad_y_hess()` and `distrib_hess_y_hess()`
  are closed form for the multivariate gaussian, which is what a marginal
  criterion reads to estimate the covariance of a correlated random effect.
  The response Hessian is `-Sigma^-1`, so it does not depend on the
  observation and does not depend on the mean at all: every component of the
  first two involving a mean is exactly zero and the rest are one matrix
  rather than one per row.
  ⚠️ Each is checked against ONE difference of the analytic quantity below it,
  never two in a row. A nested reference reported gaps of 0.3 on correct code.

* `mv_derived()` on a fixed-parameter wrapper delegates to the family instead
  of falling to the base method, which reports the distinct entries of the
  covariance: a centered prior was being read on a scale its own family does
  not use. The Jacobian keeps the columns of the free parameters alone.

# distributions7 0.21.0

* `fixed()` accepts a MULTIVARIATE family, in a third wrapper class beside
  the continuous and discrete ones. Holding the mean components at zero
  leaves the matrix parameter alone, which is what a centered prior on a
  random effect is. Every method splices and delegates as the other two do,
  and the generics a multivariate family rejects by design -- the
  distribution function, the quantile -- are inherited unregistered and go on
  rejecting. `has_mv_support()` and `has_mv_grad_y()` unwrap first, so a
  wrapper's delegation does not turn a family's refusal into a TRUE.

* `distrib_cross_y()` is closed form for the multivariate gaussian, in both
  parametrizations: `Sigma^-1 e_j` for the mean and `Sigma^-1 A_k w` for the
  matrix, one n-by-p matrix per parameter. It was refused for want of a
  consumer that fixed its shape, and a penalty whose prior is this family is
  that consumer. Agreement with numDeriv on the analytic response gradient is
  1e-9 to 1e-10 at p = 2, 3, 4.

* `distrib_expected_hessian.MvStudentTDistrib` has a help page again: its
  roxygen block had fused with the one above it, which had lost its function,
  so the page had never been written.

# distributions7 0.20.0

* The gaussian's derivatives are written in `z = (y - mu)/sigma` and
  `1/sigma`, never in a positive power of the scale.

  Written as `(res^2 - sigma^2)/sigma^3`, the score loses its denominator to
  overflow before the ratio itself becomes unrepresentable: at
  `sigma = 8e102` it returned exactly 0 where the value is `-1/sigma`, which
  on the link scale is -1, and at `1e200` it returned `NaN`. Zero is what a
  stopping rule reads as stationarity, so a run that had wandered out there
  reported `converged = TRUE` at a point that is not a maximum. The second
  derivative carried `sigma^4` and failed from `1.2e77`, the fourth
  `sigma^6` and from `4.4e61`. The variance chart is rewritten in `1/v` for
  the same reason; the precision chart already had no such form.

  Every component now agrees with the expression it replaces wherever that
  one held, and with the algebra beyond it. What remains out of reach above
  `1.3e154` is the LINK-SCALE second order, and there the cause is the chain
  rule rather than the kernel: it forms `(h')^2`, which overflows, against a
  component of order `1/sigma^2`, which underflows. The result is `NaN`
  rather than a plausible number, and a test pins that.

* The beta-binomial's mass stays a probability at any shapes.

  `lchoose(n, y) + lbeta(y + a, n - y + b) - lbeta(a, b)` evaluated left to
  right adds a term of order one to a beta function of order `1e23`, which
  annihilates it, and the subtraction then leaves exactly zero. Every mass
  came back as one, the support summed to 11 instead of 1, and the
  log-likelihood was 0 -- which beats any real fit. Fitted to binomial data,
  which asks the shapes to run to infinity at a fixed ratio, one start in
  six landed there and won.

  The shifts are integers, so each log-gamma difference is an exact sum of
  logarithms and the mass follows from three of them without forming any
  quantity larger than `n log(alpha + beta)`. The route is taken when the
  ordinary one is no longer accurate, at a threshold derived from the size
  of the terms it cancels rather than chosen. The mass now sums to one to
  1e-13 at shapes from 1 to 1e23 and reaches the binomial limit to 1e-13.

* `fit_distrib()` rejects a run of a discrete family whose log-likelihood is
  positive.

  A product of probabilities cannot exceed one, so such a run has left the
  region where the mass function is computable rather than found a better
  fit. It is discarded before the comparison, not ranked below the others:
  the objective is what breaks ties, and a number that is not a likelihood
  wins every tie it is allowed to enter. A log-likelihood of exactly zero is
  left alone, being what a legitimate fit against a support boundary
  reports.

* Six more families carry a method-of-moments estimate, taking the count to
  37 of 42: `skewnormal1` and `skewnormal2` from the first three moments,
  `weibull3` and `pig2` by carrying the estimate of the chart they
  reparametrize across its map, and `betabinom1` and `betabinom2` by reading
  the intra-class correlation off the variance. Each recovers its parameters
  to better than one per cent on 4e5 draws. The five without one --
  `skewt`, `pseudohuber`, `gengamma1`, `gengamma2`, `enet` -- have no closed
  inversion of their moments and keep the interpretation route.

  A family carrying a fixed constant spells it in its name, as in
  `"beta-binomial [size=10]"`, so the bracketed part is now dropped before
  the name is used as a lookup key: without that no beta-binomial of any
  size could match its own entry.

# distributions7 0.19.0

* A univariate family starts from the DATA, not from a draw over its
  parameters' domains.

  The base `distrib_start()` never looked at `y`. That is harmless while the
  response is of order one and fatal when it is not: on a response of mean
  919 and standard deviation 169 the draws are of order one, the first step
  is taken where the residuals are hundreds of standard deviations wide, and
  the scale runs to the largest representable double. Measured on a
  gaussian, `fit_distrib()` recovered N(5, 2) and N(50, 20) and FAILED on
  N(500, 200) -- a threshold in the scale of the data, not in the family.
  Downstream, `statmod(y ~ t, gaussian1_distrib(), Nile)` returned an
  intercept of 227.9 and a slope of exactly 0 where `lm` gives 1056.4 and
  -2.71.

  What makes one method serve every family is that they already declare
  `params_interpretation`. A parameter meaning a location starts at the
  sample median, one meaning a spread at the sample standard deviation or
  its square, one meaning degrees of freedom at the value the sample
  kurtosis implies, and anything else -- a shape, a dispersion, a
  probability -- keeps its draw, being of order one whatever the data. A
  family declaring nothing recognizable loses nothing. Values are clamped
  strictly inside their bounds, a sample median being able to land on a
  support boundary.

  Now recovered across five decades of scale, N(5,2) to N(50000,20000), and
  on gamma, lognormal, Student t, Weibull, Poisson and negative binomial
  responses centred in the hundreds. `statmod()` on the Nile agrees with
  `lm` to 4.3e-16.

* `moment_estimates()` gives a family its own method-of-moments estimate,
  which `distrib_start()` returns as the first starting value and refines by
  maximum likelihood from there. The interpretation route above is what a
  family without one falls back to. Each entry is checked against 2e5 draws
  from a known parameter, so the inversion is pinned to the family's own
  generator rather than to the algebra it was written from -- which is what
  catches a variance function transcribed from the wrong parametrization.

# distributions7 0.18.0

* `distrib_kernel(distrib, param)` returns the log-density, the score and
  the curvature in one parameter's unconstrained scale, with the family's
  methods and the link's resolved once and the chain rule applied to the
  single component wanted rather than to all of them.

  The generic route is right almost everywhere: it validates its arguments,
  aligns `theta` by name, dispatches, and assembles every component of the
  order asked for. A recursion that calls back once per observation can
  afford none of that, and a score-driven filter evaluates its score at a
  predictor it has just produced, so the call cannot be vectorized away.
  Measured on one observation, the kernel is 5.2x the generic for a
  gaussian score and 5.8x for its curvature; on a family whose own method
  is expensive the share is smaller (1.7x for the skew t).

  The caller takes on what the generic was doing: `theta` must already be
  in the family's order with values of an acceptable length, and nothing is
  checked. The inverse link is still clamped strictly inside its bounds,
  which is a correctness property and not an optimization -- an unclamped
  `exp(-800)` is zero, and a gaussian with a scale of exactly zero is not a
  distribution.

# distributions7 0.17.0

* The link-scale assembly no longer rebuilds its index layout on every call.
  `to_link_scale()` computed, for each component of the requested order, the
  multi-index list, a `unique` and a `tabulate`, and then a `sort` and a
  `paste` per term of the nested sum to spell the key of the parameter-scale
  component it needed. None of that depends on the values; it is a function
  of the parameter names and the order alone, and `link_scale_layout()` now
  computes it once and caches it.

  It was found by profiling a fitted score-driven model, where `paste`,
  `sort` and `unique` were among the leaders of the self time: a filter
  reaches this once per observation per iteration, so the names were being
  respelled millions of times in one fit. Measured end to end, a gas(1,1)
  fit at 2000 observations went from 109 to 85 seconds, and an ordinary
  gaussian fit at 100000 observations from 1.00 to 0.89. The components are
  bit-for-bit what they were, checked on four families at all four orders.

# distributions7 0.16.0

* `fit_distrib()` reports the observed information whenever the expected one
  has no closed form and the fit did not compute it.

  Which matrix the standard errors came from was decided by comparing
  `method` with the string `"newton"`. That argument has accepted an
  optimizers7 OBJECT since the delegation to that package, and an object is
  normalized internally to `"custom"`, so `newton()` failed the test and the
  expected information was assembled anyway -- by quadrature, for a family
  that does not write the expectation out. On a user-defined Gompertz, whose
  score grows like `exp(by)` while the density decays like
  `exp(-eta exp(by))`, that quadrature does not converge:
  `method = "newton"` fitted in 0.15 seconds and `method = newton()`, the
  same algorithm on the same data, had not returned after five minutes. The
  surrounding `tryCatch` does not help, a quadrature that fails to converge
  raising nothing.

  The condition now asks what is actually wanted: the expected information
  when the fit used it (Fisher scoring) or when the family writes it out and
  it costs one evaluation, and the observed one otherwise. Families that
  carry a closed-form expectation are unaffected, to the digit.

# distributions7 0.15.0

* `enet_distrib()` has a closed expected information, where it took one from
  the sampling fallback. Every piece was already in the family.

  In the two rates the density is an exponential family with sufficient
  statistics `-|z|` and `-z^2/2`, so `log Z` is its cumulant generating
  function and the information in `(a, c)` is exactly the Hessian of `log Z`
  -- which `.enet_logz_derivs()` computes for the observed Hessian already.
  The map to `(lambda, alpha)` is bilinear, so the information transforms by
  `J' I J` with no second-derivative term.

  The location is where the expected and the observed part company, and the
  reason is the kink. The observed second derivative in `mu` is `-c`, which
  misses the point mass `d sgn(z)/dz = 2 delta(z)` the density carries at its
  own location, exactly as the Laplace does; the information there is the
  variance of the score,

      I_mu_mu = a^2 + 2ac E|z| + c^2 E[z^2]
              = a^2 - 2ac dlogZ/da - 2c^2 dlogZ/dc,

  since `E|z| = -dlogZ/da` and `E[z^2] = -2 dlogZ/dc`. It reduces to
  `lambda^2` at `alpha = 1`, which is the Laplace's, and to `c` at
  `alpha = 0`, which is the Gaussian's, and both are asserted against the
  families that own them.

  Validated against the definition -- minus the outer product of the score,
  integrated against the density and SPLIT AT THE KINK, since a rule that
  straddles the corner measures the corner: 2.3e-14 to 1.1e-6 over five
  settings, the worse ones being the quadrature and not the formula. A Monte
  Carlo of 4e6 draws agrees within its own standard error at all of them.

* The multivariate validator's own gradient and second-derivative stencils
  take their nodes, weights and step from numericals7 rather than writing out
  a central difference with a step of their own. The rules are identical, so
  the checks report what they reported.

# distributions7 0.14.0

* The mixed grid reaches every shape that needs no new algebra. A family
  written out in its own parametrization with its own kernels --
  gaussian2, gaussian3, laplace2, invgauss2 -- is not a
  `reparametrize()` wrapper and does not inherit its methods, but its map
  onto a parent is already tabulated in `reparam_maps.R`, so the same
  chain rule applies and only the registration was missing.

* A family with no location keeps every scale formula: the derivation
  never used the location, only that sigma is a scale, so with
  `z = y/sigma` they hold as written. That covers the exponential, the
  Weibull's scale and the generalized Pareto's, and the exponential is
  closed outright, having nothing but a scale.

* The lognormal closes by the transformation: it is a gaussian at
  `t = log y`, and the transformation carries no parameter, so `t` does
  not move with theta and every theta-derivative is the gaussian's own
  at that point. What the response derivatives carry is the Jacobian.
  lognormal2 follows through its map, which is what closing a parent is
  for.

* `distrib_cross2_y()`, `distrib_grad_y_hess()` and
  `distrib_hess_y_hess()` reach 21 families each, up from 13. The eleven
  that remain are not location-scale and have no transformation onto
  one.

# distributions7 0.12.0

* A reparametrized family carries the whole mixed grid through its map:
  `distrib_cross2_y()` by the same first-order chain rule
  `distrib_cross_y()` takes, and the two second-order generics by the
  ordinary two-term expansion, which needs the map's second partials and
  the parent's first-order components. `reparam_tables()` already keys
  both, so nothing new is differentiated.

* The laplace was registered for `distrib_cross2_y()` and not for the two
  second-order generics. A census found it: a family half-registered
  answers at one order and falls back at the next without anything
  failing.

* The three generics now reach 13 families each, up from 9, 8 and 8.

# distributions7 0.11.0

* The location-scale identity closes the new mixed generics for the
  families it applies to. Where the response enters only through
  `z = (y - mu)/sigma`, every derivative of the response gradient and of
  the response curvature in the location and the scale is that family's
  own y-derivatives times a power of sigma, so nothing new is derived:
  `distrib_cross2_y()`, `distrib_grad_y_hess()` and
  `distrib_hess_y_hess()` go from 2, 1 and 1 closed families to 9, 8 and
  8. A family with a shape parameter beyond the two keeps its
  location-scale pairs closed and differences the rest, as the
  first-order block already does.

# distributions7 0.10.0

* `distrib_grad_y_hess()` and `distrib_hess_y_hess()` close the mixed
  grid: one or two derivatives in the response and TWO in the
  parameters. They are what the SECOND derivative of a marginal
  criterion needs of a penalty. Closed forms for the gaussian, one
  central difference of the analytic first-order component otherwise,
  and `fixed()` delegates them like every other derivative.

# distributions7 0.9.1

* `fixed()` delegates `distrib_cross2_y()`, subset to the free
  parameters as it does every other derivative. Without it every
  penalty built on `fixed(gaussian1_distrib(), mu = 0)` -- which is
  what a ridge and a random effect are -- reached the numerical
  fallback while the closed form sat one class away.

# distributions7 0.9.0

* `distrib_cross2_y()` completes the mixed surface: two derivatives in
  the response and one in each parameter, which is how the curvature of
  a log-density in the response moves with the parameters. A penalty is
  a negative log-density evaluated at the coefficients, so this is what
  a marginal likelihood needs to differentiate the determinant of a
  penalized information. Closed forms for the gaussian and the Student
  t, and one central difference of the analytic `distrib_hess_y()` for
  everything else.

# distributions7 0.8.0

## Derivatives of the distribution function

* `skewnormal1_distrib()` is closed in the shape as well as in the
  location and the scale, at every order, and `skewnormal2_distrib()`
  follows through the map. The distribution function is
  \eqn{\Phi(z) - 2T(z, \alpha)} and Owen's \eqn{T} has elementary
  partial derivatives,
  \eqn{\partial T/\partial h = -\varphi(h)(\Phi(\alpha h) - 1/2)} and
  \eqn{\partial T/\partial\alpha = \varphi(h)\varphi(\alpha h)/(1+\alpha^{2})},
  so the integral in its definition is differentiated away at the first
  order and never has to be differentiated again. Everything above is a
  product of normal densities, Hermite polynomials and a rational
  function of the shape.

* The check that the first identity is the right one is that it returns
  the density: \eqn{\partial F/\partial z = 2\varphi(z)\Phi(\alpha z)}. At
  \eqn{\alpha = 0} the location and scale components agree with the
  gaussian's to 1e-15, and the gaussian reaches them by another route.

* Every cdf order now leaves the same eight families on the stencil, and
  all eight are mathematical obstructions: the derivative of an
  incomplete gamma or beta in its shape is hypergeometric (gamma,
  chi-squared, beta, generalized gamma) and the von Mises distribution
  function is itself a quadrature. A test asserts that list, so a family
  added without a route joins it visibly.

# distributions7 0.7.0

## Derivatives of the distribution function

* `invgauss1_distrib()` and `enet_distrib()` have closed derivatives at
  all four orders, from one route: a distribution function of the form
  \eqn{c_0 + \sum_k s_k e^{w_k}\Phi(x_k)} is a Leibniz split between
  the weight and the tail, with a Faa di Bruno on each side. A family
  supplies, per term, the partial derivatives of the log weight and
  those of the argument.

  The inverse gaussian's \eqn{\Phi(a) + e^{c}\Phi(b)} has all three of
  \eqn{a}, \eqn{b} and \eqn{c} separable in the mean and the
  dispersion, so their mixed partials are products of one-variable ones.
  The elastic net's halves are truncated Gaussians, and its \eqn{s} and
  \eqn{x} are likewise separable in \eqn{\lambda} and \eqn{\alpha};
  its weight is written through the Mills ratio the family already
  carries, \eqn{w = -\log M(x) + x^{2}/2}, so its derivatives come from
  the same \eqn{G} the density uses.

* The weight is never formed on its own. \eqn{e^{2/(\phi\mu)}} is
  \code{Inf} at ordinary settings -- 2000 in the exponent at
  \eqn{\mu = 0.01}, \eqn{\phi = 0.1} -- exactly where \eqn{\Phi(b)}
  underflows, so the two are combined as
  \code{exp(w + pnorm(x, log.p = TRUE))} and the fourth derivative comes
  back finite.

* `invgauss2_distrib()` takes its Hessian through the mapped route too.
  Registering only the gradient there was right while the parent
  differenced its own second order and is not now that it does not.

* Every one of the four cdf surfaces now leaves the same nine families
  on the stencil, and all nine are obstructions or correct refusals: the
  derivative of an incomplete gamma or beta in its shape is
  hypergeometric (gamma, chi-squared, beta, generalized gamma), the von
  Mises distribution function is itself a quadrature, and
  `skewnormal2_distrib()` is refused by the gate while its parent's
  shape components are differenced.

# distributions7 0.6.0

## Derivatives of the distribution function

* `gpd_distrib()` has closed derivatives at all four orders, its
  Hessian included, from the exponential-survival route. The form is
  what makes it work: writing \eqn{u = \xi q/\sigma} and
  \eqn{\Lambda(u) = \log(1+u)/u}, the exponent is
  \eqn{L = -(q/\sigma)\Lambda(u)} and carries no division by the shape,
  so the exponential limit \eqn{\xi \to 0} is an ordinary point of the
  formula rather than a branch of the code.

* `Lambda` and its four derivatives come from the recursion
  \eqn{u\Lambda^{(r)} + r\Lambda^{(r-1)} = (-1)^{r-1}(r-1)!/(1+u)^{r}}
  above \eqn{\lvert u\rvert = 1/2} and from the Taylor series below it.
  The crossover is measured, not chosen: the recursion divides by
  \eqn{u} and subtracts nearly equal quantities, and its fourth
  derivative is wrong by a factor of \eqn{10^{39}} at
  \eqn{u = 10^{-14}}, by 1.7 at \eqn{10^{-4}} and by \eqn{3\times10^{-8}}
  at \eqn{10^{-2}}, while the two agree to \eqn{10^{-16}} at the switch.

* At \eqn{\xi = 0} the scale components equal the exponential family's
  exactly, at every order. That family reaches them through
  \eqn{L = -q/\mu} and shares no arithmetic with the series, so it is
  the reference near zero, where a stencil in the shape is not one: at
  \eqn{\xi = 10^{-6}} the step is a thousand times the value.

* A negative shape bounds the support above, at \eqn{\sigma/\lvert\xi\rvert},
  and the derivatives are zero past it. The support is declared by the
  family rather than read off the bounds, the endpoint depending on a
  parameter.

# distributions7 0.5.0

## Derivatives of the distribution function

* A family whose survival function is an exponential of something
  elementary now gets all four orders from one identity:
  \eqn{S = e^{L}} gives \eqn{\partial^{I}F = -S\,B_{I}(L)}, the complete
  Bell polynomial in the partial derivatives of \eqn{L}. A family states
  \eqn{L} and its partials and nothing else. `exponential_distrib()` and
  `weibull1_distrib()` are served by it, and `weibull3_distrib()` follows
  through the reparametrization wrapper; the route replaces the two
  written-out orders each carried before, which it reproduces.

* `laplace2_distrib()` reaches the two new orders as the Laplace at
  \eqn{\sigma = 1/\lambda}, through the mapped route.

* The upper tail of an exponential-survival family is exact wherever the
  logarithm is representable. `log S` is `L`, so its derivatives are
  `L`'s own and need no division by `1 - F`, which is exactly one in
  double precision past `q/mu = 37`: the first derivative of an
  exponential's log survival at `q = 700` was `Inf` and is now `700`.

* `weibull_cdf_deriv()` is gone with the two methods that were its only
  callers.

# distributions7 0.4.0

## Derivatives of the distribution function

* The third and fourth derivatives are closed for thirteen more families,
  taking the count from four to seventeen. Three routes did it, and none
  derived anything new.

  A family written as a map of another carries the parent's through one
  Faa di Bruno pass; orders one and two already did this and the new
  orders use `chain_assemble()`, the enumeration the reparametrized
  parameter derivatives run on, so no second copy of the partition sum
  exists. That closes `gaussian2_distrib()`, `gaussian3_distrib()` and
  every `reparametrize()` wrapper whose parent is exact.

  The mapped route now admits a transformation of the response as well
  as of the parameters. A lognormal is a gaussian at `log q` and the
  transformation carries no parameter, so the derivatives in the
  parameters are the gaussian's with the point substituted; that closes
  `lognormal1_distrib()`, and `lognormal2_distrib()` follows from it
  through the wrapper.

  `gumbel_distrib()` joins the location-scale families, and
  `student_t1_distrib()`, `pseudohuber_distrib()`,
  `skewnormal1_distrib()` and `skewt_distrib()` get the location and
  scale components from that construction with the shape components
  still differenced, as at the two orders below.

* The gate is the one orders one and two use: a chain rule is taken only
  when the parent is exact at every order up to the one asked for, so a
  differenced quantity is never reported as a closed form. Measured
  against the partial-expectation integral, which shares no code with
  any of the three routes, the fully closed families agree to 4e-15 and
  the partial ones to the stencil's own 3e-6.

# distributions7 0.3.0

## Derivatives

* The third and fourth response derivatives are closed for every
  continuous family. Eighteen were taking a finite-difference stencil,
  all of them families whose response is not a pure location; each
  already carried a closed second response derivative, and the third is
  the same elementary function differentiated once more. The log-density
  of each is a sum of terms in `log(y)`, `log(1 - y)`, a power of `y`,
  a logarithm of an affine function of `y`, or a cosine, so the terms
  are written once and each family is a sum of them.

* `reparametrize()` carries the third and fourth response derivatives to
  the parent, as it already carried the first and second: the map acts on
  the parameters and the derivative is taken in the response, so the two
  do not interact.

* The generalized Pareto's coefficient is written as `xi^k + xi^(k-1)`
  rather than `(1 + 1/xi) * xi^k`. It is the same number and stays finite
  as the shape goes to zero, where the family is exponential and every
  order above the first is exactly zero.

  Measured against one differentiation of the analytic second response
  derivative, the third order agrees to 1e-11 and the fourth to 4e-5,
  which is each reference's own accuracy. A test walks the namespace and
  fails if a continuous family is left on the stencil.

# distributions7 0.2.0

## Derivatives

* Every one of the 46 families is now analytic to fourth order in the
  parameters. The three that were still on the numerical fallback --
  `dirichlet_distrib()`, `multinomial_distrib()` and
  `mvstudent_t_distrib()` -- have closed third and fourth derivatives.
  The two simplex-valued log-densities are a sum of terms each depending
  on one coordinate, so the chain rule is one univariate partition sum
  per coordinate; the Student t splits into the mean-and-matrix part,
  which reuses the gaussian's expansion of the derivative of an inverse,
  and a part in the degrees of freedom, which is elementary. In each
  case the same assembly run at orders one and two reproduces the
  hand-written score and information, to 9e-16 and 4e-15.

* The mixed derivative `distrib_cross_y()` of `skewnormal1_distrib()`
  and `pseudohuber_distrib()` is closed in the shape parameter as well
  as in the location and the scale. The response reaches the skew
  normal's shape only through `alpha * z` and the pseudo-Huber's only
  through `D`. The skew normal assembles all three components from a
  single evaluation of the inverse Mills ratio rather than through
  `distrib_grad_y()` and `distrib_hess_y()`, which evaluate it twice
  more: 77.5 ms to 19.4 at n = 1e5.

## Internals

* `mvg_ptensors()` takes the pieces rather than the distribution, so the
  Student t reuses one copy of the expansion, and its accessor answers
  for the empty multiset.

# distributions7 0.1.0

## Plots

* `plot()` on a univariate distribution draws one curve per element of a
  parameter given as a vector, so
  `plot(gaussian1_distrib(), list(mu = 0, sigma = c(1, 2, 4)))` is three
  densities on one panel. The settings are separated by color and by line
  type together, which keeps them apart in a printed copy that has no color;
  the parameters that vary are named in a legend, placed on whichever side
  the mass leaves emptier, and those held fixed are stated once in the
  title. A discrete family is drawn as several sets of stems, shifted
  sideways so that equal masses at one support point stay countable, and
  separated by symbol rather than by line type: a dashed stem reads as a
  broken one, and at a support of any size the panel fills with fragments.

  Every component must have length one or the same `k`. A length that merely
  divides `k` is rejected rather than recycled, since a partial setting is
  far more likely to be a mistake than a request. The horizontal window
  covers every setting rather than the first.

  This meaning is available because a plot has no data to recycle against;
  the density and derivative generics read a vector component as one value
  per observation, which is a different question asked of the same object.
  A multivariate family, whose picture is already a matrix of panels with no
  axis left to overlay settings on, rejects a vector component instead.

## Families

* One name per parametrization. A family with several parametrizations
  carries a number on each -- `gaussian1`/`gaussian2`/`gaussian3`,
  `gamma1`/`gamma2`, `negbin1`/`negbin2` after Cameron and Trivedi,
  `weibull1`/`weibull3` after gamlss's WEI and WEI3 with `weibull2`
  deliberately empty, and ten further groups. The reference index lists
  what is present.

* Univariate families added: weibull, gumbel, skewnormal, skewt,
  exponential, geometric, chisq, betabinom, NB1, the generalized Pareto,
  the generalized gamma, von Mises, the Poisson-inverse Gaussian in both
  gamlss parametrizations, and `enet_distrib()`, the elastic-net prior.

* Multivariate families: the gaussian and the Student t, whose matrix
  parameter comes from `parameters7`, and the Dirichlet and the
  multinomial, whose simplex parameter does. The matrix parameter is
  flattened into scalars with identity links, so `align_theta()`,
  `deriv_names()`, the link scale and `fit_distrib()` need no special case;
  the constraint lives in the structure, where it belongs. The base class
  sits beside `continuous_distrib` and `discrete_distrib` rather than under
  either, the one-dimensional defaults registered there -- a cdf by
  quadrature, a quantile by root finding -- having no counterpart in
  several dimensions.

* `mv_summary()` reports the quantities a reader reads rather than the
  coordinates: standard deviations and correlations with delta-method
  standard errors, each interval built on the scale that keeps the quantity
  in its own set (log for a standard deviation, Fisher's z for a
  correlation) and mapped back. A precision parametrization adds the
  conditional standard deviations and the partial correlations.

* Wrappers: `fixed()`, holding parameters at known values and the only one
  that removes parameters; `folded()`, for the absolute value; and
  `truncated()`, `zero_inflated()`, `zero_adjusted()` and
  `transformation()` from earlier. Neither zero wrapper can be stacked on
  the other -- truncating at zero cancels `(1 - zeta)` between the
  numerator and the truncation constant, so `zeta` leaves the likelihood
  entirely -- and both are rejected by the constructor, along with a
  discrete parent carrying too few support points to identify the extra
  probability.

## Derivatives

* Every univariate family is analytic to fourth order, observed, except
  the skew t's components in `nu`, which cannot be: the density carries
  `T_{nu+1}` and the derivative of a Student t distribution function in its
  degrees of freedom has no elementary form. Those come from one five-point
  stencil applied to an analytic quantity, never from a difference of a
  difference.

* Orders three and four are closed form for every wrapper. Each wrapper's
  log-likelihood is the parent's log-density plus, or instead of, `log L`
  for some parameter-dependent `L`, so two partition sums cover all of
  them: the complete Bell polynomial, and the moment-to-cumulant relation,
  which needs only the ratios `d^B L / L`.

* `distrib_deriv3_cdf()` and `distrib_deriv4_cdf()` complete that
  sequence. The routes are the ones the orders below use: a discrete family
  sums the identity exactly and a continuous one applies one product
  stencil to its analytic distribution function. Nothing new was derived --
  the quantity summed is the complete Bell polynomial and the conversion to
  the log scale the moment-to-cumulant relation, both already in the
  package for the wrappers and now in `partition_sums.R`, where they belong.
  The general forms reproduce the written-out orders one and two exactly,
  on both tails and both scales.

* `distrib_deriv3_y()` and `distrib_deriv4_y()` do the same for the
  response. A family whose response enters only as `y - mu` needs no
  formula of its own: `d^k l / dy^k = (-1)^k d^k l / dmu^k`, so fourteen
  families inherit these orders from derivatives they already have, often
  compiled ones. The identity is checked at orders one and two, where both
  sides are written independently and it holds exactly. A family on a half
  line takes one stencil of the order asked for, with the step halved
  because the stencil reaches two steps either side of a support boundary.

* `distrib_grad_cdf()` and `distrib_hess_cdf()`, the derivatives of the
  distribution function, which is what a censored likelihood and a quantile
  residual's standard error need. Closed form for twelve of the original
  fourteen families; gamma and beta have none, the shape direction being
  hypergeometric.

* `distrib_cross_y()`, the mixed response-parameter derivatives, closed
  form for every continuous family. Where the response enters only through
  `z = (y - mu)/sigma` the identity `d2l/dy dmu = -l_yy` and
  `d2l/dy dsigma = -z l_yy - l_y/sigma` closes nine families at once.

* `reparametrize()`, building a family from another through a map, with
  Faa di Bruno over partitions and the map's partials as hand-written keyed
  tables (`map_derivs`); without tables, one stencil per partial. Measured
  against a family written out in full it costs 1.5x at the gradient and
  6.6x at order four, so it is the user-facing route and new families are
  written out.

* Jets are removed from every production path. Generic jet composition
  measured at 2x to 36x the hand-written closed forms on the PIG kernels;
  the mechanical transcriptions survive in the tests as independent
  references.

## Fitting

* `fit_distrib()` delegates its optimization to `optimizers7`. Fisher
  scoring is `newton()` with the expected information passed as `he`, and
  `method` accepts an optimizer object, so nothing in this package
  implements a descent loop.

* `fisher_scoring(approx =, nsim =, criterion =, maxit =)` replaces the
  loose `approx` and `nsim` arguments: how the expected information is
  approximated is a property of Fisher scoring and had no business sitting
  beside optimizers that never look at it. A strategy chosen where it would
  be ignored is rejected.

* `maxit` and `tol` leave the signature. With `method = <an optimizer>`
  they were silently discarded, so a call setting both got no complaint and
  no effect from the second; the budget and the stopping rule now live on
  the method.

* The objective is `-l(eta)/n`. The maximizer and every Newton step are
  unchanged, the factor canceling in `H^-1 g`; what changes is what a
  threshold means, an absolute gradient tolerance on a summed score asking
  of a sample of ten million an accuracy per observation ten million times
  finer than of a sample of ten. `ll_hat` and the information are
  recomputed unscaled at the optimum, so `logLik()`, AIC, BIC and every
  standard error are untouched.

* `distrib_start()` computes a starting value from the data. The default is
  the old random draw from the parameter domains; the multivariate gaussian
  returns the sample mean and covariance, its own maximum likelihood
  estimate for an unstructured matrix.

* The restart loop keeps the best result rather than the last: a converged
  run beats a non-converged one, and among runs of equal status the lower
  objective wins. A run that reached a point is no longer discarded as a
  failure.

* The default tolerance is `1e-6`, for the reason recorded in optimizers7's
  own notes: the attainable gradient is bounded below by
  `sqrt(2*lambda*eps*|f*|)`, and a log-likelihood is of order one at its
  maximum.

* `confint()`, with `scale = c("parameter", "link")`, recomputing at any
  level from the stored estimate and standard error. The link-scale
  interval is the one computed and the parameter-scale table is its image
  under the inverse link, so `print()` shows both.

* A fit records its elapsed time, accumulated over every start and every
  fallback, and the score per observation at the point it stopped.

## Numerical layer

* The enumerations, stencils, batched quadrature and series, and special
  functions come from `numericals7`. `expectation()`'s integrand contract
  is elementwise in the response and the parameters jointly, so every
  parameter combination shares one batched call.

* Random numbers come from generalized ratio-of-uniforms, recentered at the
  mode and normalized. A density diverging at an edge is transformed away,
  and the exponent is measured by the same probe that detects the
  divergence rather than searched for: Gamma shape 0.4 went from 27 ms a
  draw to 0.8 microseconds. Discrete families invert the cumulative table,
  which is exact.

## Validation and documentation

* `check_distrib()` runs thirteen numerical checks on a continuous family,
  twelve on a discrete one and a nine-check battery on a multivariate one.
  It is aware of atoms, so a mixed distribution such as
  `zero_adjusted(gamma1_distrib())` is not reported as four failures on
  correct code; it allows for a kink where `params_smooth` declares one;
  and a gradient made five per cent wrong is still caught.

* `mv_reference_draw()` supplies the proposal the normalization check
  integrates against. The gaussian proposal does not fail loudly on the
  Dirichlet -- `chol()` accepts the singular covariance -- and returns
  2.0e-08 for an integral that is 1.

* Every family's constructor page displays its density, taken from the
  book's catalog, whose transcription is checked against `distrib_pdf()` at
  every render.

* Three vignettes -- defining a distribution, fitting a model, derivatives
  and the link scale -- a README with badges, a pkgdown site and continuous
  integration on five platforms.
