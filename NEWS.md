# distributions7 0.57.0

* **`fixed()` delegates `distrib_deriv3_y()` and `distrib_deriv4_y()` to its
  parent**, as it already did for the first two response derivatives. Without
  the two methods a fixed continuous family took the base class's stencil on
  the log-density, so `fixed(student_t1_distrib(), mu = 0)` -- the
  heavy-tailed prior of `penalties7` -- differenced where its parent carries
  the closed form. The route is pinned by an identity against the parent,
  which a stencil cannot satisfy. Continuous families only: a discrete family
  has no response derivative to delegate to.

# distributions7 0.56.0

* **`check_distrib()`'s cdf row leaves out a grid point in the last places of
  a non-zero bound, as its response row already does.** The response row has
  excluded a draw within `128 |b| eps` of such a bound since 0.55.0, where the
  spacing of doubles is absolute and the reference compares nothing; the cdf
  row, whose grid is the deciles, did not, on the reading that those deciles
  sit nowhere near a bound. Measured, that is true of every family but one:
  beta2 is 4.9e+04 units away at worst and the two von Mises 3.2e+14, while
  **beta1 reaches zero units — the grid point IS the bound.** It carries
  \eqn{\mathrm{shape}_2 = (1-\mu)\phi}, which inside `generate_random_theta()`'s
  own box reaches \eqn{0.1 \times 0.1 = 0.01}, and its upper decile is then
  \eqn{1 - O(e^{-c/\mathrm{shape}_2})}, which stops being representable below
  about 0.06. Over 20000 draws of that generator the row **failed 59 times
  (0.30%) on analytic code, with statistics up to 8.006**, every one between
  0.5 and 85 units of the bound, and a further 185 (0.93%) had the grid
  collapse onto the bound and the row was not emitted at all.

* Measured on the four worst parameter values of that census, before and
  after: 8.006, 3.451, 1.86e-01 and 4.43e-04, all FAIL, become 7.6e-09,
  3.9e-07, 4.2e-08 and 1.1e-09, all OK, one grid point of fifteen left out and
  counted in each. Over 400 draws of `generate_random_theta()` the row now
  fails **0 times against 1** and is **emitted every time against 4 silent
  drops**, its worst statistic 8.5e-06 against 1.96e-02. ⚠️ Where the grid is
  not near a bound nothing moves: beta1 at `mu = 0.5, phi = 0.5` reads
  1.664e-10 and at `mu = 0.3, phi = 3` 7.959e-11, identical before and after.
  A discrete family compares an exact difference of the distribution function
  against the mass, with no step at all, so nothing is left out there. The
  non-regression net that captures `check_distrib()`'s table and printed text
  for six families over five seeds sees no leaf move from this at all: the
  parameter values it carries put no decile near a bound.

* **A numerical derivative in the PARAMETERS chooses its step near a bound,
  as the response direction already does.** `fd_steps()` cut its step to 49%
  of the distance to the nearest finite bound, and where the log-density is
  singular in that parameter there the error of such a difference stops
  falling once the cut binds. It now offers a second candidate,
  `to_bound = "scale"`, which scales the step ON that distance, and the new
  `fd_stable_step()` evaluates the quotient at both, each also at half its own
  step, and keeps the one that agrees better with itself. Measured over a
  census of 324 cells at distances from 1 to 1e-8 from a bound, against the
  families' own analytic derivatives: the gradient goes from **199 to 307**
  components within 1e-6 and the diagonal Hessian from **106 to 238**; on a
  gamma at a dispersion of 1e-6 the gradient's relative error goes from
  3.5e-01 to 1.1e-06.

* ⚠️ **The choice is made for the whole vector, not observation by
  observation, which is the opposite of `fd_stable_quotient()` in the response
  direction.** There each observation carries its own \eqn{y} and therefore
  its own step, so a per-observation choice is a choice between two quantities
  that genuinely differ; here the parameter is usually one number for the
  whole sample and the two candidates estimate the same thing, so the
  variation of the reading across observations is rounding. Measured, the
  whole-vector choice puts 307 gradients and 238 Hessian components within
  1e-6 against 285 and 229, and on the 46 cells where the two differ it is the
  better on 44 (on the gpd the per-observation choice returns `Inf` where the
  whole-vector one reads 8.3e-08).

* ⚠️ **Where the two steps coincide nothing moves, by construction rather than
  by tolerance.** That is every parameter with no finite bound and every
  positive parameter at or above one: the clamped quotient is returned at once
  and is `identical()` to what the previous release produced, 61 cells of 324
  taking that path. Elsewhere the choice costs four quotients where one would
  do.

* **What the choice costs in time**, at 500 observations, before and after.
  Where the two steps coincide, a gaussian pays the extra call and nothing
  else — gradient 0.188 to 0.230 ms, Hessian 0.375 to 0.426, fourth
  derivative 1.172 to 1.148 — and a gamma at a dispersion of 2 reads 0.499 to
  0.549, 1.178 to 1.230 and 7.120 to 7.047. Where they differ it is the
  measured 4 quotients against 1: a gamma at a dispersion of 0.3 goes 0.522 to
  1.318 ms on the gradient, 1.154 to 1.957 on the Hessian and 6.551 to 10.174
  at fourth order, and a beta, whose two parameters both carry two finite
  bounds, 0.637 to 1.632, 1.464 to 2.295 and 1.382 to 1.867. The higher orders
  cost proportionally less because the choice is paid once per parameter.
  `fd_stable_step(need_value = FALSE)` is what keeps the coinciding case free
  there: those callers read only the step, so nothing is evaluated at all —
  without it a gamma at a dispersion of 2 paid 9.087 ms for the two Hessians
  of a choice that was never in doubt.

* **All seven callers take the chosen step**: `numerical_gradient()`,
  `numerical_hessian()`, `numerical_deriv3()`, `numerical_deriv4()`,
  `numerical_cross_y()`, `numerical_cross2_y()`, `numerical_theta2_y()` and
  `numerical_dexpected_hessian()`. Each chooses with its own `h_rel` and its
  own stencil; a single choice made at first order and reused would cost 2
  cells of 324 on the Hessian, the two tests agreeing on 244 of the 263 cells
  where the steps differ at all. For the Hessian and the two higher orders the
  step is chosen once per parameter and then used by every component that
  parameter enters, the mixed ones included, so the choice is paid once per
  parameter rather than once per component.

* ⚠️ **`numerical_dexpected_hessian()` does not choose its step under
  `approx = "mc"`**, and the link scale, which carries no finite bound, has
  nothing to choose. A Monte Carlo expected information is not the same number
  twice, so the reading would be of that noise: measured on pig1, the quotient
  moves by 2.75 relative between two evaluations at one step, against 17.3
  between the step and its half.

* ⚠️ **THE CHOICE MAKES THE NUMERICAL INFORMATION A DISCONTINUOUS FUNCTION OF
  THE POINT, and it costs a fit nothing.** Two nearby points can take
  different candidates, so the Hessian jumps where a fixed step would move
  smoothly: measured on a density-only Gompertz whose two parameters both lie
  below one, perturbing a parameter by \eqn{10^{-10}} moves the numerical
  Hessian by 5.7e-05 relative where the clamped step moves it by 9.2e-07. A
  margin in the choice -- keeping the scaled step only where its reading is
  several times better -- is free in accuracy (307 and 238 unchanged at
  margins of 2, 5 and 10) and removes the jump at the smallest perturbations
  but not at \eqn{10^{-8}}, so it is not taken and the discontinuity is
  declared instead. What licenses that is a measurement on fits rather than on
  points: over nine density-only families whose parameters fall below one, all
  nine converge under both rules, **in the same 140 iterations**, with vcov()
  against the analytic family at 5.49e-07 and 5.59e-07 and its reproducibility
  from a displaced start at 4.00e-07 and 4.04e-07. An optimizer does not
  sample the Hessian at points \eqn{10^{-10}} apart.

* ⚠️ **What this costs in the ORDINARY regime, which the census near the bound
  could not see.** Over 2539 components at the parameter values
  `generate_random_theta()` produces, 174 improve and 44 worsen, the median
  ratio among those that move being 0.53. The worst case is `enet`'s
  `alpha_alpha_alpha` at third order, from 1.16e-08 to 7.34e-08, and no
  component crosses 1e-4 or 1e-3 in either direction; at 1e-6, one leaves and
  six enter. `check_distrib()` inherits the change, since it takes these
  functions as its reference in the parameter direction, and 20 of its rows
  move on `beta1` and `enet`: `deriv4` on beta1 from 2.38e-07 to 1.91e-06 and
  on enet from 7.67e-07 to 3.12e-06, every verdict unchanged.

# distributions7 0.55.0

* **`check_distrib()`'s response-derivative row no longer fails correct
  families near a support bound.** Its reference was a central difference
  whose step was cut to 0.49 of the distance to the bound, and on a
  log-density carrying \eqn{(a-1)\log y} the relative error of that difference
  is \eqn{(h/d)^2/3} once the cut binds: a plateau of 9.4e-02 on the first
  derivative and 1.4e-01 on the second, which halving the step cannot expose
  because both steps are cut alike. Measured over 32 continuous families,
  three parameter values and five seeds, the row failed 35 times in 480, all
  on families singular at a bound (gamma1 14, gamma2 10, chisq 4, beta1 3,
  beta2 3, lognormal1 1), whose derivatives agree with Richardson to 1e-10.
  The reference is now `fd_stable_quotient()`, internal: it evaluates the
  difference at the clamped step and at a step scaled on the distance to the
  bound, each also at half its step, and keeps per observation the one that
  agrees better with itself, the clamped step on a tie. Over the same census
  the row fails nowhere, its worst at 5.0e-05.

* ⚠️ Neither step is right everywhere, which is why both are read. On a sweep
  to within 1e-8 of a bound the new reference fails none of 264 points of the
  families singular there, where the clamped step fails 165 on the first
  derivative and 198 on the second; on the families smooth at the bound it
  fails 1 and 47 of 144, where the clamped step fails 0 and 43. It costs four
  differences where one did, except on a family with no finite bound, where
  the two steps coincide and the clamped difference is returned at once.
  `fd_steps_y()` takes `to_bound = c("clamp", "scale")` for the second step,
  its default unchanged.

* **Near a bound that is not zero the reference divides by the steps it
  actually takes, and a draw in the last places of that bound is left out.**
  Near zero the spacing of doubles is relative; near any other bound it is
  absolute, and within about 1e-10 of it a step that is not a whole number of
  units is not the step its evaluation points lie at, while the scaled step
  falls below one unit and rounds to zero. The quotients now divide by
  \eqn{(y+h)-y} and \eqn{y-(y-h)}, through the internal `fd_first_taken()` and
  `fd_second_taken()`, and no candidate step is shorter than
  \eqn{2\lvert y\rvert\epsilon}. Their relative error still grows as about
  \eqn{5.4/d^2} in the distance \eqn{d} counted in units of the spacing, so a
  draw closer than \eqn{128\lvert b\rvert\epsilon} to a bound \eqn{b \ne 0},
  256 units, is left out of the row and counted in its detail, as a Monte
  Carlo draw exactly on a bound is. Measured on 2000 samples of 100 draws from
  `beta1_distrib()`, singular at 1: at `mu = 0.5, phi = 0.5` the row failed in
  486 samples, fails in 10 on the steps taken and in none with the draws left
  out (30 of 200000); at `mu = 0.9, phi = 1`, where 2.5 per cent of the draws
  land exactly on 1, it failed in 1999, fails in 1285 and then in none, with
  4.4 per cent of the draws left out. Over the census of 32 families no
  verdict moves and no draw is left out, the worst statistic going from
  3.0e-05 to 5.0e-05, and on the distribution function's grid no statistic
  moves by more than 1.4e-10.

* **The row comparing the distribution function with the density uses the
  same reference.** It differenced at 1e-5 max(1, |x|), a step not kept inside
  the support, and failed on the grid points nearest a bound: 3 rows over the
  same families and parameter values (gamma1 twice and gamma2 once, the worst
  at 1.1e-02), and 66 of 144 and 85 of 264 points of a sweep toward the bounds
  of the families smooth and singular there. It fails none of either now, its
  worst at 9.2e-05. ⚠️ The relative step moves from 1e-5 to
  \eqn{\epsilon^{1/3}} = 6.06e-06 on every family, so the row's statistic moves
  wherever it is computed: on `laplace_distrib()`, `laplace2_distrib()` and
  `enet_distrib()`, whose density has a kink the grid crosses, it falls by the
  ratio of the two steps, 0.6055, and on `gaussian1_distrib()` from 7.3e-12 to
  6.1e-12.

* **A check whose statistic has no value is reported as not run.** A row
  whose statistic came out `NaN` or `NA` without its computation raising
  compares nothing, and reporting it as a failure made a correct family fail:
  the generalized Pareto at \eqn{\xi < -1/2}, where the expected information
  does not exist, failed in 5 runs of 5. Such a row is left out of the table
  and recorded with its reason in the attribute `"skipped"`, and `print()`
  lists it as not run. An infinite statistic stays in the table and fails,
  being a component that overflows where its reference does not; so do a
  density or a distribution function that is not finite, whose checks are
  defined on those values, and a check whose computation raised.

* **A Monte Carlo draw exactly on a finite bound is left out of the
  expected-information row and counted in its detail.** A log-density
  singular at the bound is not finite there, nor is the score, and one such
  draw made the whole Monte Carlo mean `NaN`: on `beta1_distrib()` with
  `nsim = 2e5`, 5 to 11 draws per seed land exactly on 1 in double precision,
  and the row came back `NaN` in 3 seeds of 3 where without them its statistic
  reads 0.51, 0.75 and 1.9. A point mass the family declares through
  `distrib_atoms()` stays in, and a score that is not finite at an interior
  point still leaves the row without a statistic.

* ⚠️ **Nothing a fit reads moves.** The battery of twelve fitted models (168
  comparisons) and the net of the smoothers' move to numericals7 (529 leaves)
  are `identical()`. In the net of this batch every check table moves, 88
  leaves: the distribution-function row on all 44 entries; the response row on
  gamma1 in five seeds of five (four from a failure to a pass), on beta1 in
  five of five (all five from a failure to a pass, three of them only once
  the quotients divide by the steps taken), and by rounding on the families
  with no finite bound, every verdict unchanged (gaussian1 in two seeds,
  laplace, laplace2, enet and the injected families built on laplace); and the
  expected-information row on beta1 in two seeds, from a failure at `NaN` to a
  pass with one and two draws left out. The 1.4e-10 above is the census's:
  near a non-zero bound the steps taken move the distribution-function row
  more, beta1 at `mu = 0.5, phi = 0.5` reading 1.7e-10 where the nominal
  quotient read 5.3e-09 and the old step 1.6e-05. The fits of that net are
  identical against statmodels7 0.124.0; what moves there with 0.125.0 is
  statmodels7's own.

# distributions7 0.54.0

* **A family declares WHERE its log-density is not smooth, and the order of
  differentiability follows from the declaration rather than being a second
  thing to get wrong.** `kink_decomposition()` is a new generic returning the
  one non-smooth piece as `c(theta) * phi(v(y, theta))`: `phi` names the
  composition, `v` its argument and the derivatives of that argument in the
  parameters, `coef` the multiplier in front. The base method returns `NULL`,
  so every family is smooth until it says otherwise, and three declare —
  `laplace`, `laplace2` and `enet`, all for the same reason, an absolute value
  of the residual.

* `params_order()` reads the order off it: the order of the composition for
  every parameter the kink's argument moves with, and `Inf` for the rest. It is
  DEDUCED and not declared, which is why `v` carries its own derivatives — the
  set of non-smooth parameters is \eqn{\{p : \partial v/\partial\theta_p \ne
  0\}}, measured at probe values rather than assumed. A Huber likelihood puts
  its kink at \eqn{\lvert y-\mu\rvert - k\sigma}, which moves with the
  location, the scale and the cut-off together, and that test names all three
  without the family listing them.

* The order the five compositions leave, each measured as the first derivative
  that jumps across the origin and kept only where the jump does not shrink
  with the step: 0 for \eqn{\lvert v\rvert} and \eqn{(v)_+}, 1 for
  \eqn{(v)_+^2}, 2 for \eqn{(v)_+^3}, and -1 for \eqn{1\{v>0\}}, which is not
  continuous at all and is reported apart — there a smoothing would repair a
  discontinuity of the log-density rather than a missing derivative, and the
  unsmoothed problem has no well-defined maximum.

* `check_kink()` holds a declaration to the family that declares it, which is
  what makes it worth trusting. The load-bearing check is the JUMP: the score
  of \eqn{c\,\phi(v)} jumps by \eqn{c\,\Delta\phi'\,\partial v/\partial\theta_p}
  as `y` crosses the kink, and the family's own score is written independently
  of the declaration, so the two can be compared. Measured, `laplace` and
  `laplace2` agree at 0.00e+00 and `enet` at 1.0e-06. Two further checks
  isolate `dv` against a difference of `v`, and assert that the parameters the
  declaration calls smooth really are — without the third a declaration naming
  too few parameters would pass.

* ⚠️ **Nothing consumes any of this and no behaviour moves.** Against a battery
  of twelve fitted models spanning `linpar`, `s`, `te`, `random`, `ridge`,
  `lasso`, a distributional model, a Poisson, a Gamma, `gas`, `seg` and
  `regime`, the log-likelihood, the coefficients, the hyperparameters, the
  effective degrees of freedom, `vcov`, the convergence flag and the
  certificate are `identical()`: 168 comparisons and no difference.

* ⚠️ **A WRAPPER REPORTS `NA` AND NOT `Inf`.** No wrapper propagates the
  declaration yet, so `truncated(laplace_distrib())` returns `NULL` from
  `kink_decomposition()` while `params_smooth` still records the kink.
  `params_order()` reports `NA` there, which says the order has not been
  established; `Inf` would say the parameter is smooth, which it is not. A
  consumer must read `NA` as unusable. `fixed(laplace_distrib(), mu = 0)` is
  the case that genuinely becomes smooth, the wrapper removing the only
  non-smooth parameter, and it reports `Inf`.

* ⚠️ **The elastic net's kink is real at every `alpha` and its SIZE is
  \eqn{2\lambda\alpha}**, which is why a detector reading the second Bartlett
  identity finds the family in one sweep and not in another: at a small `alpha`
  the family is nearly Gaussian and the missing curvature disappears into the
  scale of the rest of the matrix. The declaration says so where a measurement
  of the curvature cannot.

* **`log(1+w) - w` is computed from the ratio where forming `w` loses it, and
  three call sites carried the same defect.** `psi_Ew(z - 1)` represents
  \eqn{1+w} only to an absolute `eps`, so its error is about \eqn{eps/(2z)}
  and grows without bound as `z` goes to zero; below \eqn{z \approx 1.1e-16}
  the subtraction has lost the argument outright, `z - 1` is exactly `-1`, and
  `log1p(-1)` is `-Inf`. `psi_Ew2(opw, w)` takes the ratio as well, which both
  call sites already compute from their own inputs -- `y/mu` for the gamma,
  \eqn{(y+\theta)/(\theta+\mu)} for the negative binomial -- so no
  approximation is involved and the crossover is DERIVED: setting the two
  spellings' errors equal gives \eqn{1/(2\,opw) = \lvert\log opw\rvert}, i.e.
  \eqn{opw \approx 0.35}.

* ⚠️ **It was reachable from an ordinary fit and the degradation was gradual
  before it was catastrophic**: 2.9e-10 at `z = 1e-8`, 2.6e-05 at 1e-14, then
  `-Inf`. A gamma of shape 0.236 puts 6.3 observations of 20000 below 1e-15,
  and one such observation takes the whole summed score with it -- measured,
  13 of 40 seeds of `generate_random_theta()` left `gamma1`'s `phi_phi` at
  `-Inf`, and 0 of 40 do now. `negbin2`'s dispersion score read `-Inf` at a
  zero count with \eqn{\theta \le 1e-16}, where the value is
  \eqn{\log\theta + 1}; \eqn{\theta} rides a log link whose floor is 1.9e-77,
  and a count model with no overdispersion drives it there.

* The repaired values agree with a central difference of the analytic
  log-density at 1.0e-08 -- the reference's own accuracy -- from `z = 1e-6`
  down to 1e-40, and `negbin2` agrees with the closed form at **0.00e+00** for
  every \eqn{\theta} at or below 1e-16.

* ⚠️ **The cost was measured rather than assumed, and it is nothing where
  nothing was wrong.** On the summed score of a gamma sample the switch is
  bit-identical at shapes 1 and 4, moves 2.0e-12 at shape 0.5, and repairs
  3.2e-03 at shape 0.236. Against a battery of twelve fitted models one
  quantity moved: the certificate's mode error on a Gamma fit, 2.290433e-14 to
  2.290433e-14, a relative 3.8e-10 on a number eleven orders below the 1e-3 it
  is compared against. The log-likelihood, the coefficients, `vcov`, the
  effective degrees of freedom, the convergence flag and the certificate state
  are `identical()`.

* ⚠️ **The third site is the C-callable twin**, `d7_ccallable.cpp`, which
  mirrors `gamma1_parts` expression by expression and is held to it by
  `identical()`; repairing one and not the other would have turned that twin
  test red. The rule this file records for a shape of mistake -- grep for the
  shape, then measure each occurrence -- is what turned one repair into three,
  and the `negbin` site was found that way rather than by anything failing.

# distributions7 0.53.0

* The centered skew normal's map to the direct parametrization is **finite at
  every skewness the link can produce**, where it returned `Inf` from a
  predictor of about 36 upward. `bounded_link()` keeps `gamma1` strictly
  inside `(-0.9952717, 0.9952717)` and hands back the last double below the
  bound; `1 - delta^2`, formed there as `1 - (mu_z/b)^2`, rounds to exactly
  zero, so the shape came out `1/0` where its value is `-1.36e8` and
  perfectly representable. `sn_one_minus_delta2()` evaluates the same quantity
  as `1 - (|gamma1|/gamma_max)^(2/3)` through `log` and `expm1`, an exact
  identity with no cancellation left in it.

* **The derivative tables carried the same defect one step later.**
  `md_skewnormal2()` formed `D = b^2 + (b^2-1) r^2`, which is that same
  `b^2 (1 - delta^2)` written so that it cancels: at the last representable
  skewness it reaches `-1.11e-16` where the quantity is `9.42e-17`, so
  `sqrt(D)` was `NaN` and every derivative order with it. Both places read one
  shared helper now, so the value and its derivatives follow from one
  algebraic form, where the map and the tables had been written from two.

* A skewness the map cannot carry is **reported in the centered family's own
  terms**. Every probability function of `skewnormal2_distrib()` evaluates the
  parent at the mapped parameters, and the parent validates what it is handed
  against its own domains, so a caller who wrote `gamma1` read an error naming
  `alpha` and `"skew normal1"`, neither of which the call mentions.
  `sn2_theta()` is the one place the delegation happens and now checks what it
  is about to pass on, through `sn2_reject_unmappable()`, which names
  `"skew normal2"`, the centered parameter responsible and its value.

* ⚠️ **No ordinary fit moves**, and the controls are transcriptions of the
  expressions replaced rather than the package read back to itself. The map
  agrees with the old one to `1.15e-14` over 27 settings of
  `(sigma, gamma1)`; the four derivative orders, summed as a fit reads them,
  agree **exactly** at `gamma1` of `0.4`, `-0.7` and `1e-4`, the two forms
  coinciding bit for bit there, and to `7.85e-11` at `0.9`. Against an
  independent route, the round trip through `skewnormal1_distrib()`'s own
  `skewness()`, the map reproduces the skewness it was given with a gap of
  `0` at the last representable value. `check_distrib()` passes all thirteen
  checks and the guard costs 0.06 per cent of one gradient.

* ⚠️ The check states a property of the delegation rather than guarding a
  value the public surface reaches today: `gamma1` is bounded on the
  constructor and every generic validates it before dispatch, so a skewness
  outside its domain is already reported correctly. What stays reachable is a
  `sigma` and a `gamma1` each inside its own domain whose implied scale or
  location leaves the doubles, from `sigma = 1.1e308` upward.

# distributions7 0.52.0

* The skew t's derivatives in `(mu, sigma, alpha)` are **closed form** at
  orders three and four: ten of the twenty and fifteen of the thirty-five. The
  location and the scale reach the log-density only through
  `z = (y - mu)/sigma` and the shape only through `w = alpha u(z)`, so two
  observations close the block. Differentiating in the shape never leaves it,
  `d^c/dalpha^c l = u^c Lam^(c)(w)`; and for any function of `z` alone,
  `d^a_mu d^b_sigma F = (-1)^(a+b) sigma^-(a+b) P_{a,b}(z)` with
  `P_{a,b+1} = (a+b) P_{a,b} + z P'_{a,b}`. Every piece is elementary:
  `u` and `log t_nu` are rational up to one square root, and
  `Lam^(k)` comes from the Riccati recursion for `Q = t_{nu+1}/T_{nu+1}`,
  exactly as a skew normal's inverse Mills ratio does. `skewt_msa_tower()`,
  `skewt_msa_component()` and `skewt_msa_derivs()`.

* The ten fourth-order components carrying **exactly one** `nu` are one
  five-point difference along `nu` of the closed-form third derivative beside
  them (`skewt_msa_nu1()`), which is the rule `distrib_hessian()`'s mixed
  components already follow, read one order up. The generic construction took a
  mixed **second** difference of the Hessian instead, and a second difference
  amplifies rounding by `h^-2`.

* `numerical_deriv3()` and `numerical_deriv4()` take `skip`, a character vector
  of components to leave `NULL` for a caller that supplies them in closed form.
  The names and their order do not move. With `skip = NULL`, the default, both
  are `identical()` to what they returned before.

* The licence for the orders that cannot be checked against a hand-written form
  is that the **same assembly at orders one and two reproduces the family's
  own score and Hessian**, derived separately and already under
  `check_distrib()`: measured at 7.7e-17 to 1.7e-16. Against Richardson applied
  to the hand-written Hessian, over fifteen settings of `(nu, alpha)` with
  `nu` from 3 to 50, the closed forms read 3.6e-08 at order three and 2.4e-07
  at order four, which is the reference's own floor.

* What it buys, measured over that grid against the route it replaces: order
  four went from **2.4e-05 to 2.4e-07**, and the ten components carrying one
  `nu` are between **20 and 203 times** closer to Richardson. Order three shows
  no gain there and is not claimed to: a first difference of an analytic
  Hessian already sits at Richardson's own floor, so the reference cannot see
  the difference. What order three buys is the two orders above it: the ten
  fourth-order components carrying one `nu`, and the fifth order below.

* ⚠️ It also closes a **nested difference** at the fifth order, which nothing
  had noticed. `numerical_deriv5()` differences `distrib_deriv4()`, and for a
  component free of `nu` that fourth derivative was itself a SECOND difference
  of the Hessian in the same variables -- a difference of a difference, which
  this package forbids everywhere. Measured with two accuracies of the same
  first difference, sharing only their centre node: over the 21 order-5
  components free of `nu` the two agree to **2.5e-10** where the route replaced
  gave **4.8e-03, 3.6e-03 and 1.2e-02** at `nu` of 5, 8 and 20. Seven orders of
  magnitude, and it is what order three buys, since order three shows no gain
  of its own.

* And it is **cheaper**, the fifteen closed components no longer being
  differenced at all: at `n = 20000` `distrib_deriv4()` costs 5.05 s where the
  generic construction alone costs 17.79 s, which is where the page's earlier
  "about sixteen seconds" came from.

* `has_exact_deriv4()` answers **`TRUE`** for the skew t, which is a decision
  recorded here rather than a consequence of the closed forms above: with part
  of the family exact, what the predicate should say about it is a choice. It
  owns its fourth-order method, so the default owner reading gives that answer
  and no override is registered for it. The predicate is one logical for the
  whole family and says the family **supplies** the derivative, not that every
  component of it is a closed form: fifteen of the thirty-five are, the twenty
  carrying `nu` are single stencils, and a model needing a fourth derivative of
  this family gets one.

* ⚠️ What that costs is a measurement rather than a promise, and it depends on
  `nu`. `check_distrib(orders = 1:5)` now emits an order-5 row for the skew t;
  `check_distrib()`'s own `rel()` floors its denominator at 1, so for a small
  component it reads an **absolute** error, and the stencil noise in `nu` falls
  as `nu` grows. Swept over `nu` in 3, 5, 8, 20, 50 and `alpha` in -2, 0.5, 3,
  the row reads **3.7e-03 to 4.6e-03 at `nu` = 3**, where it fails against the
  default tolerance of 1e-3, is marginal at 5, and reads 4.4e-05 to 4.3e-04 at
  8 and 5.1e-06 to 3.7e-05 at 20, where it passes with orders of margin. Both
  ends are asserted in the tests. `orders` still defaults to `1:4`, so no
  existing caller meets that row.

* ⚠️ Read **per component against its own scale** the fifth order of this
  family is untrustworthy at every `nu`, and that is the statement to carry
  rather than the row's verdict: the two accuracies of `numerical_deriv5()`,
  which share only their centre node, disagree by 5e-02 to 1.1 over that grid,
  while the 21 components free of `nu` agree to **2.5e-10** -- the range every
  family whose fourth order is genuinely analytic sits in, 1.3e-10 to 1.2e-07.

# distributions7 0.51.0

* `to_link_scale()` reaches the **fifth** order, so the `stop()` that capped it
  at four is gone. Three edits and no new machinery: `link_scale_layout()` and
  `deriv_index_list()` were already generic in the order, and what capped the
  surface was `bell_partial()`'s table, `inverse_link_derivs()`'s switch and
  `link_scale_lower_orders()`'s ladder. The fifth derivative of an inverse link
  the first of those consumes is what linkfunctions7 0.4.0 delivered.

  The order-5 partial Bell row was checked before it was written, against a
  construction of the same polynomials as the coefficient of \eqn{t^5/5!} in
  \eqn{(\sum_m x_m t^m/m!)^j/j!}, which shares no arithmetic with a table:
  the two agree **exactly**, and the coefficients sum to the Bell number
  \eqn{B_5 = 52}. Both are tests, along with a negative control that puts one
  coefficient 5 per cent out.

* **The link-scale fifth now has two routes, and they agree.** `distrib_deriv5()`
  differentiates in \eqn{\eta} the order-4 component already on the link scale,
  which needs neither \eqn{h_5} nor the Bell row; `to_link_scale()` can instead
  carry the parameter-scale fifth over by Faa di Bruno at order 5, which needs
  both. They share no arithmetic beyond the parameter-scale fourth, so their
  agreement is a check rather than one expression twice: measured over six
  families, between 2.9e-09 and 9.6e-08.

  ⚠️ **What ships is unchanged, and the measurement is why.** Against
  Richardson on the analytic fourth the differencing route reads 1.1e-10 to
  2.7e-08 and the chain-rule route 3.0e-09 to 8.1e-08, so the one already in
  place is between 3 and 33 times the closer -- the chain mixes a numerical
  fifth with four analytic lower orders and accumulates more rounding than one
  difference does. It is also the more robust of the two, needing no \eqn{h_5},
  so it does not degrade for a user-defined link whose fifth derivative is a
  numerical fallback. The chain route is 1.4x to 1.75x faster, on a quantity
  costing two to nine milliseconds, which does not buy back the digits.

* An identity link leaves every order alone, which is the control that costs
  nothing and would catch a chain applied where it should not be: both routes
  reproduce the parameter scale exactly (0.000e+00).

# distributions7 0.50.0

* `distrib_deriv5()`, the **fifth**-order derivatives of the log-likelihood,
  on the parameter scale and on the link scale. The order exists because each
  order of differentiating a score-driven filter's predictor through its own
  recursion draws in one more order of the family -- the curvature reaches the
  third, the directional third derivative the fourth, and the outer Hessian of
  a model carrying such a term the fifth. Writing it inside the term that
  needs it would be the private helper this package exists to abolish.

  No family computes it in closed form. Every one reaches `numerical_deriv5()`,
  which applies **one** central difference to the analytic fourth. That is the
  rule the rest of the derivative surface obeys and it is worth saying that
  `numerical_deriv4()` deliberately does not: it differences the analytic
  Hessian twice, which it can afford because a second difference of an exact
  quantity is still only two orders removed from one. A fifth built the same
  way would not be.

  The nodes, the weights and the step are all \pkg{numericals7}'s.
  `fd_derivative()` is not called directly for the reason
  `numDeriv_grad()` already records -- its `f` maps a vector of points to the
  values at those points, while this one reads a whole named list at each node
  -- so the offsets and weights are taken from the library instead and
  `accuracy` is an argument rather than new arithmetic. A node of zero weight
  is skipped, so the default rule costs two evaluations of the fourth order per
  parameter, and grouping the components by the index being differentiated
  makes the whole order cost `2p` rather than two per component.

  Measured against Richardson extrapolation on the analytic fourth, over seven
  families on both scales: between **7.6e-11 and 1.9e-08**. So the numerical
  fifth is not a placeholder for the architecture -- it is accurate enough to
  build the outer Hessian on.

* The link scale needs neither the fifth derivative of a link nor an order-5
  entry in `bell_partial()`. Differentiating in \eqn{\eta} the order-4
  component **already on the link scale** -- itself analytic, being a Faa di
  Bruno chain of analytic pieces -- gives the fifth on that scale directly. The
  `stop()` at `link_scale.R:150` therefore stays where it is until analytic
  fifth derivatives arrive.

* `has_exact_deriv4()` says whether the quantity being differenced is itself
  analytic, which is what the fifth order is worth. The default reads the class
  that owns the registered `distrib_deriv4()` method; a wrapper overrides it and
  asks its parent, a wrapper's fourth derivative being a partition sum over the
  parent's first four. `check_distrib()` gains an order-5 row, emitted only
  where the predicate is `TRUE` and reached by `orders = 1:5`, so no existing
  caller's row count moves.

  ⚠️ **A family may own its fourth-order method and still not be exact in every
  component**, which is the shape `expected_hessian_exact()` already records
  for the pseudo-Huber and skewnormal2. `skewt_distrib()` is that case here,
  and the predicate answers `TRUE` for it by the owner reading: it owns the
  method, while its density carries \eqn{T_{\nu+1}}, whose derivative in the
  degrees of freedom has no elementary form, so its higher orders are built
  from single stencils on analytic quantities. Differencing them amplifies
  whatever noise they carry by \eqn{1/h}. Measured at
  \eqn{\mu = 0, \sigma = 1, \alpha = 0.5}, the order-5 components from the
  three-point rule and from the five-point one -- which share only their centre
  node -- disagree by 1.8e-02 at \eqn{\nu = 5}, **1.7e-01** at 8 and 2.2e-02
  at 20, where every family whose fourth order is genuinely analytic agrees
  between 1.3e-10 and 1.2e-07.

  ⚠️ A census over the other 38 univariate families found one further
  disagreement, at the generalized Pareto, and **it was the probe and not the
  family**: `generate_random_theta()` had put the shape near the crossover of
  its own series branch. At fixed values all three routes agree between 1.6e-10
  and 6.8e-08 and its order 4 is right to 2.0e-08.

* A family defined with nothing but `distrib_pdf()` still answers at the fifth
  order -- the promise that a distribution needs only a density is not
  withdrawn -- but there the fourth is itself a difference, so the fifth is a
  difference of a difference. `has_exact_deriv4()` is `FALSE` for it and for
  every wrapper of it, and `check_distrib()` emits no verdict rather than one
  nobody earned. Measured on a density-only gaussian the fifth is 3e+04
  relative to the analytic family's, which is not a derivative of anything.

# distributions7 0.49.0

* THE PURE-DISPERSION DERIVATIVES AT ORDERS THREE AND FOUR no longer cancel,
  in BOTH negative binomials and both observed and expected. Each family
  tends to a Poisson from its own side -- negbin2 as `theta` runs away,
  negbin1 as `theta` goes to zero and the size `r = mu/theta` runs away -- so
  every derivative in the dispersion vanishes there and is written as a
  difference of terms that agree to leading order, with the chain dividing by
  the highest power at these orders.

* negbin2, OBSERVED. The shift is `y`, a count, and the term in
  `(y - mu)/s^k` that accompanies each polygamma difference is itself a sum
  over the same range, so the two merge and factorize: with `t = th + j` and
  `s = th + mu`,

      l_ttt  = sum_{j<y} 2(mu-j)(s^2 + s t + t^2)/(t^3 s^3)
               - mu^2 (3th + mu)/(th^2 s^3)
      l_tttt = -sum_{j<y} 6(mu-j)(t+s)(t^2+s^2)/(t^4 s^4)
               + 2 mu^2 (6th^2 + 4 th mu + mu^2)/(th^3 s^4)

  Measured, the spelling replaced is 4.15e-02 out at `theta = 1e7` at order 3
  and 3.09e-02 at order 4, where these are within 1.80e-13 over 105 cells
  spanning `y` to 500, `mu` to 100 and `theta` from 0.1 to 1e9. The sum costs
  `y` terms, so it is taken only where the direct form loses digits, which is
  where `theta` is large against `y` and the sum is therefore SHORT: the
  threshold is 100, measured -- the direct form is within 3.9e-10 at
  `theta/y <= 100` and reaches 5e-02 at 1e6.

* negbin2, EXPECTED. The closed term is `mu` times a constant and therefore
  `E[sum_{j<Y} c]`, so it merges into the summand rather than being added
  afterwards, exactly as `nb_E_ltt` does at order two. The composition it
  replaces LOSES ITS SIGN at `theta = 1e7`, reading -3.23e-32 where the value
  is +4.80e-34 at order 3 and +9.55e-39 where it is -2.88e-40 at order 4;
  after, `value * theta^p` settles on its asymptote to eight figures
  (47.99999, 48.00000 at `mu = 4`; 29999.99, 30000.00 at `mu = 100`).

* negbin1. ⚠️ THE SIZE NEED NOT APPEAR AT ALL, which is what the recursion in
  the powers of `r` was fighting. Since `lgamma(y+r) - lgamma(r)` is
  `sum_{i<y} log(r + i)` and `r = mu/theta`, the `y log(theta)` each such term
  carries cancels EXACTLY against the `C(theta)` part, leaving

      l = sum_{i<y} log(mu + i theta) - log(y!)
          - (mu/theta) log1p(theta) - y log1p(theta),

  whose variable does not run away. Verified against `dnbinom` to between 0
  and 2.2e-11, the last being R's own loss of the Poisson limit. Every
  derivative is then elementary, and one pass over `i` serves every component
  of an order, the denominator power being the order itself; the pass is
  compiled (`negbin1_psums_cpp`). Measured at `mu = 4`, `y = 3`, the route
  replaced reads 2.770e-01 at `theta = 1e-3` where the value is 2.797e-01,
  -2.37 at 1e-4, -1.97e+09 at 1e-6 and -1.46e+17 at 1e-8, and at order four
  8.97 at 1e-3 where the value is -1.59; the new one converges on 9/32 and
  -1.5984375 with the gap falling like `theta` (1.6e-2, 1.6e-3, 1.6e-4,
  1.6e-5).

* `nb1_M_derivs()` is the one composite piece, `M(theta) = log1p(theta)/theta`
  and its derivatives, which has a removable singularity at zero: the
  recursion divides by `theta` and the series converges only below one. ⚠️ The
  crossover is MEASURED -- the two agree to between 1.9e-16 and 3.4e-13 over
  `theta` from 0.01 to 0.5 and each fails on its own side, the recursion by
  2.3e+08 at order four at `theta = 1e-6` and the series by 1.1e-02 at 0.8.
  It is the shape the generalized Pareto's `Lambda` already carries.

* ⚠️ The negbin1 route is taken only below `nb1_exact_cut()`, which is 1: the
  sum costs `y` terms an observation against the O(1) of the polygamma route,
  and above the crossover that route is within 4e-11. On a real design of
  57600 counts summing to 6.24e6 the whole pass is a few milliseconds
  compiled, against 23 s for the same pass written in R -- which is why it is
  compiled and not left in the assembly.

* `test-negbin-higher-cancellation.R` pins all three, each with the form it
  replaced as a negative control.

# distributions7 0.48.0

* The expected polygamma quantities of BOTH negative binomials are summed
  as DIFFERENCES through an exact recurrence, where before they were
  expectations with `psi^(n)(theta)` subtracted at the call site. The two
  agree to leading order wherever the family tends to a Poisson -- negbin2
  as `theta` runs away, negbin1 as `theta` goes to zero and the size
  `r = mu/theta` runs away -- so the subtraction lost the answer. The shift
  is `Y`, which is a count, and for an integer shift

      psi^(n)(x + k) - psi^(n)(x) = (-1)^n n! sum_{j<k} 1/(x + j)^(n+1),

  terms of ONE SIGN, accumulated beside the mass the loop already carries.
  It is the device `psi_diff.h` states for the observed derivatives, applied
  to the expected ones.

* It is also the whole cost. One `trigamma` or `psigamma` call per term
  becomes one division, and where `psi_T_rest` is below its crossover it
  becomes one division instead of TWO `trigamma` calls, the second always at
  the same argument. Measured on the 57600 fitted means of a real design,
  one evaluation of each helper:

  |                            | before   | after   |        |
  |----------------------------|----------|---------|--------|
  | `nb_E_trigamma` (negbin2)  | 3.570 s  | 0.070 s | 51x    |
  | `nb_E_psigamma` order 3    | 6.500 s  | 0.070 s | 93x    |
  | `nb_E_psigamma` order 4    | 6.600 s  | 0.080 s | 83x    |
  | `nb1_E_Pr`, theta = 0.2    | 0.360 s  | 0.060 s | 6x     |
  | `nb1_E_Pr`, theta = 100    | 40.14 s  | 0.250 s | 161x   |

  The negbin1 spread is `psi_T_rest`: with the repetition loop sized by
  elapsed time, a term costs 1451 ns at `r = 0.5` and 472 ns at `r = 6`
  against 17.4 ns and 4.5 ns here, and 8.2 ns against 2.2 ns above
  `r = 100`, where that helper is already its asymptotic series.

* ⚠️ THE MASS SEED CARRIED A CANCELLATION OF ITS OWN, in `nb_E_trigamma`
  and `nb_E_psigamma`. `log P(Y = 0)` was written
  `theta * (log(theta) - log(theta + mu))`, and at `theta = 1.585e5` with
  `mu = 0.1` those two logarithms are 11.9736 apiece while their difference
  is 6.31e-07: seven digits lost on the seed, which the multiplicative
  recurrence then carries to every term. It is `-theta * log1p(mu/theta)`.
  Measured with nothing else changed, that one line reads 3.52e-07 out
  where `log1p` reads 5.74e-13. negbin1's seed was already correct.

* What the two together buy, over 277 cells spanning `mu` from 0.1 to 1e5
  and `theta` from 0.05 to 1e6, against the same recurrence summed in R
  with `dnbinom`: the route replaced is 2.08e-03 out at `mu = 0.1`,
  `theta = 5.012e5` at order 2, 1.04e-03 at order 3 and 6.95e-04 at
  order 4, where the kernels are now 3.07e-12 or better at every cell and
  every order. Away from the Poisson limit the two routes agree to 1e-10
  or better, which is what says the difference is the cancellation and not
  a change of formula.

* `test-negbin-polygamma-recurrence.R` pins it, with the route replaced
  transcribed in R as the negative control: that transcription reproduces
  the compiled measurement to the digit (2.08e-03 and 3.52e-07 at the two
  cells) and must FAIL where the kernels pass, so the tolerance cannot be
  met by a route that has gone back to subtracting at the end.

* ⚠️ TWO THINGS ARE NOT REPAIRED HERE. The OBSERVED third and fourth
  derivatives still write `R::psigamma(y + th, n) - R::psigamma(th, n)`
  directly, which is the same cancellation at the same shift -- the gap
  0.36.0 records.

* AND negbin2's EXPECTED INFORMATION no longer loses its sign, which closes
  what 0.32.0 records as unrepaired. `E[l_theta_theta]` used to be assembled
  from the polygamma difference plus `mu/(theta(theta+mu))`, each of order
  `mu/theta^2` while their sum is of order `mu^2/(2 theta^4)` -- thirteen
  digits at `theta = 1e6`, which no accuracy in the summand reaches. The
  second term is itself an expectation over the same support, since
  `E[sum_{j<Y} c] = c mu`, so the two merge into ONE accumulation whose term
  is `(theta(2j - mu) + j^2) / (theta (theta+mu) (theta+j)^2)`, a quotient of
  exact polynomials.

* The asymptote it is checked against is DERIVED and not fitted: expanding
  both pieces in `1/theta`, the `theta^-3` term contributes `mu^2/theta^4`
  and the `theta^-4` term `-3mu^2/(2 theta^4)`, so the information tends to
  `mu^2/(2 theta^4)`. Over `mu` in {1, 4, 100} and `theta` from 10 to 1e8 the
  new form is positive at all 24 cells and its gap to that asymptote falls
  monotonically -- at `mu = 1`, 5.22e-06, 3.02e-07, 1.79e-08 at `theta` 1e6,
  1e7, 1e8 -- where the composition reads -1.63e-25 at `mu = 4`,
  `theta = 1e7` and -5.63e-28 at `mu = 100`, `theta = 1e8`. Over a grid of
  56 cells the worst relative error is 4.15e-04 against the composition's
  2.47e+07.

* It costs nothing, measured BACK TO BACK: two whole fits of the same model
  run one after the other give 263.6 s for the composition and 264.8 s for
  this form, with the outer trajectory identical evaluation for evaluation.
  ⚠️ Runs taken at different moments of the same session read 270.3 and
  271.3 against 292.8 and 292.2 -- eight per cent, which is the machine and
  not the change, and only running the two arms consecutively separates them.
  Per term the two forms are within 1.4% of each other at every mean from 1
  to 5e4, with the merged one never the slower.

* ⚠️ THE TRADE, stated. Merging moves the cancellation from order `theta`
  to order `theta/mu`, so it costs accuracy in the one corner where `mu` is
  far larger than `theta`: at `theta = 100` and `mu` from 4.0e4 to 1.7e5 the
  composition reads 3.45e-12 to 2.28e-11 against the merged sum's 6.98e-10 to
  3.96e-09, a factor of 58 to 265 between two quantities that are both
  negligible for an information matrix. `test-subnormal-seed.R` carries that
  tolerance at 1e-7 where it was 1e-9, with the numbers beside it. A dispatch
  by regime would recover both and is not taken: what it would buy is 1e-11
  against 1e-9.

* ⚠️ Two references that do NOT settle this, recorded so they are not
  reached for again. A sum taken in R DIVERGES from the kernel past
  `theta = 1e5` -- its own gap to the asymptote grows from 2.02e-05 to
  1.07e-01 while the kernel's falls -- so there it is the reference that
  fails. And the Poisson limit is not a reference at all: under a Poisson
  mass the answer is exactly `3 mu^2/(2 theta^4)`, a factor of three, because
  `E[Y(Y-1)]` is `mu^2` there against `mu^2(1 + 1/theta)` here and the
  difference is precisely the term the derivation turns on. That the factor
  comes out as exactly 3 is what confirms the derivation.

# distributions7 0.47.0

* The negative binomial with a linear variance computes each quantity and
  no more. `nb1_parts()` takes the order asked for and returns after the
  first-order block, where before it wrote every component whatever the
  caller read: `negbin1_gradient_cpp` needs `P`, `Q` and the two
  derivatives of the size, and was paying for `Pr` as well, which is the
  only place `psi'` enters and below the crossover of `psi_T_rest` is two
  `trigamma` calls an observation. Measured at 28800 observations the
  gradient goes from 0.0117 s to 0.0050 s, and the value, the score and
  the observed hessian are `identical()` to the previous release over
  eight regimes of the parameters.

* `nb1_E_Pr()`, the series the expected information rests on, sums the
  mass through its own recurrence and stops on the accumulated mass,
  where before it located a far-tail quantile with `R::qnbinom` and then
  summed at least a hundred terms with one `R::dnbinom` each. That
  quantile call is the one `nb_E_trigamma()` removed one family over for
  the reason `d7_par.h` states: its search reaches `pbeta`, whose warning
  path calls into the R API and kills the process from a worker thread.
  It was also the cost, and with it gone `negbin1_expected_hessian_cpp`
  runs through `d7::par_for` like every other kernel here, bit-identical
  at any thread count. On a 28800-cell design with an offset one
  evaluation goes from 1.2767 s to 0.3356 s, and to 0.0467 s over eight
  threads. Against a reference summed in R the worst relative error over
  98 points spanning `mu` from 0.01 to 20000 and `theta` from 1e-4 to 40
  is 1.07e-09, against the old kernel's 1.04e-09.

* ⚠️ The log scale is left on the threshold the loop switches at and not
  on `exp(lpk)` merely being nonzero. The window where that matters is
  narrow and ordinary: at `mu = 918.832`, `theta = 0.5` the size is
  1837.7 and `log P(Y = 0)` is -745.1, so the mass at zero is a subnormal
  with almost no significand, and seeding the multiplicative recurrence
  there put `mu_theta` at -2.73e-02 where it is 1.2089e-04. Across that
  whole window the kernel now agrees with the reference to 3e-13.

* ⚠️ The same line was wrong in negbin2, in BOTH helpers that carry the
  recurrence -- `nb_E_trigamma()`, which the expected information rests
  on, and `nb_E_psigamma()`, which the expected third and fourth
  derivatives rest on. It was found by grepping for the shape after
  repairing negbin1 rather than by anything failing. Measured against a
  sum taken in R over the family's own mass, at `mu = 1.702e5`,
  `theta = 100`, where `log P(Y = 0)` is -744.0:
  `E[trigamma(Y + theta)]` was 3.7e-02 out, and the error runs smoothly
  up to that edge (2.3e-06 at -734, 1.4e-03 at -739). The window is
  reached by an ordinary count model with a large mean and a mild
  overdispersion, and nothing about the failure is loud.

* ⚠️ A memo went with the parallelization, and what it was worth is
  stated rather than buried. The loop this replaced carried the
  expectation across CONSECUTIVE equal parameters, which fires wherever a
  design repeats a mean in adjacent rows and nowhere at all under an
  offset -- measured on those same 28800 cells, 28800 distinct means of
  28800, so it never hit once and was what kept the loop sequential.
  Where it does fire the change costs: 100 distinct means over 28800 rows
  go from 0.0048 s to 0.4175 s sequentially, and to 0.0590 s over eight
  threads.

* ⚠️ Two things this does NOT repair, both measured and both older than
  the change. The expected information's `theta` block is three terms of
  order 1e13 summing to order 1e2 at `theta = 1e-6`, so a change to the
  series at the 1e-10 level moves it entirely and no accuracy in the
  summand reaches it. And `mu_theta` at a large mean is a difference of
  two numbers near 4/3, carrying four to five digits of cancellation of
  its own.

# distributions7 0.46.0

* The two Poisson-inverse gaussian families get ONE COMPILED KERNEL PER
  GENERIC, as every other compiled family has: `pig1_pdf_cpp`,
  `pig1_gradient_cpp`, `pig1_hessian_cpp`, `pig1_deriv3_cpp`,
  `pig1_deriv4_cpp` and their `pig2` twins, in place of the single
  fifteen-column block every method read. The shape was inherited rather
  than chosen: the first version propagated a bivariate fourth-order jet,
  which computes every partial together by construction, and when the jet
  was replaced by explicit closed forms the block's signature stayed. So
  the value cost four orders.

  Each kernel now recycles a scalar parameter and tests the support itself,
  which is what `gaussian.cpp` and `betabinom.cpp` already do, and returns a
  named list; the ten methods are one line each, with nothing between them
  and their kernel. **`pig_hd_block()` is gone**, and with it `pig_memo` and
  `pig_hd_forget()`: pig was the only compiled family carrying an R helper
  on that path. The mathematics that page documented -- the closed form of
  \eqn{S_y(\alpha)}, which kernel runs, the jet twin -- moved onto
  `pig1_distrib()`, where a reader looks for it.

  `pig1_hd_cpp` and `pig2_hd_cpp` survive, returning all fifteen columns:
  they are the reference the per-generic kernels are held to and the side
  the jet twin compares against.

  **No formula was retyped and the arithmetic is unchanged.** The body was
  reordered, each statement placed at the greater of its own derivative
  order and the order of everything it reads, and each stage returns once it
  has written its own block. Against the previous release, over both
  families and all five generics plus a scalar parameter and rows off the
  support, **twelve comparisons of twelve are `identical()`**.

  ⚠️ **What the split buys is smaller than the count of quantities
  suggests.** At n = 20000 the kernels run at 22.8, 23.0, 23.3, 23.6 and
  24.8 ms against the block's 27.9, so the value is 1.22x. The cost is not
  the derivative tables: it is \eqn{S_y}, which is \eqn{O(y)} per
  observation with three `lgammafn` calls a term and which every order pays
  in full. Measured at a response of exactly zero, where that sum is empty,
  the algebra alone runs 4.16 ms at order zero against 6.27 at order four;
  at a mean of 3 the same pair is 24.87 against 27.49, and at a mean of 30
  it is 116.03 against 117.43.

  What the split does remove is the R layer: the helper's mask, its three
  subsets, its NaN matrix and its scatter were 24 per cent of a gradient
  call at that size and are now 2 per cent.

  ⚠️ **And it costs the memo, which is stated rather than hidden.** The
  helper kept a depth-one memo of the point it last evaluated, worth 1.19x
  to 1.27x on a real regression, because one point is read several times an
  iteration. Without it the same fits cost 1.032 s and 0.848 s against
  0.926 s and 0.775 s. The duplication is not this family's and was closed
  where it happens, in `statmodels7::iwls_pieces()`: with that, the same
  fits cost 0.944 s and 0.753 s, the second being better than the memo ever
  gave. Counted, a fit went from three evaluations of the score an
  iteration to two.

  ⚠️ **The families that reach the outer-product route are FOUR** --
  `pig1`, `pig2`, `skewnormal1` and `skewt` -- and not every family, which
  is what a first draft of this entry claimed. Every other family has a
  closed-form expected information registered on its own class, so `approx`
  is never read and no score is recomputed; `pseudohuber` and
  `skewnormal2` register their own quadrature, which is not the outer
  product either.

# distributions7 0.45.0

* The log-mass of `pig1_distrib()` and `pig2_distrib()` no longer cancels in
  the Poisson limit. It carries `1/sigma + psi(alpha)`, and with
  `psi(alpha) = -alpha + log S_y(alpha)` the pair `1/sigma - alpha` is a
  difference of two quantities of size `alpha` whose value is of size `mu`.
  Both are now written in closed form: `-2 mu/(1 + sqrt(1 + 2 sigma mu))` for
  pig1, and `-mu [alpha/(mu + r)] [1 + mu/(alpha + r)]` with
  `r = sqrt(mu^2 + alpha^2)` for pig2, which follows from
  `alpha - r = -mu^2/(alpha + r)`. `psi_derivs()` reports `log S` beside
  `psi`, so recovering it by adding `alpha` back -- the same cancellation --
  is not needed.

  Measured at `mu = 0.3`: the direct form is 1.6e-4 out at `alpha = 1e12`,
  reads exactly zero at 1e16 -- a probability of one -- and +128 at 1e18,
  which is one ulp of 1e18, so the mass over the support summed to `Inf`.
  The closed form reaches `dpois()` to 5.7e-8 at `alpha = 1e6`, to 2.8e-17 at
  1e16, and holds there to `alpha = 1e300`, with the mass summing to 1
  throughout. pig1 behaves the same as `sigma` goes to 1e-300.

  It was reachable from a fit and not only from a probe. On a count model
  with a province in the dispersion equation, one province carried no
  overdispersion, `alpha` ran to the flat region, and past 1e13 the reward
  was the cancellation rather than the data: the fit parked at
  `alpha = 1.3e103` and reported a log-likelihood of -4985.04 where the same
  point evaluated correctly is -5239.635. The 254.60 is the honest
  contribution of that province's 1344 cells, -282.14, read as -27.54.

* `pig2_hd_cpp()` is finite at every `alpha` a `log_link()` can produce. Its
  derivatives divided by `r^3`, `r^5` and `r^7` and multiplied by `alpha^3`
  and `alpha^4`, and `alpha^3` passes `double.xmax` at `alpha = 5.6e102`,
  where `r_mma` read `Inf - Inf` and the whole Hessian came back `NaN`. The
  partials of `r` are written in `r_m` and `r_a`, which lie in `[0, 1]`,
  against powers of `1/r`; `r` itself comes from `std::hypot`, finite where
  `mu^2 + alpha^2` is not. The Faa di Bruno expansion is homogeneous -- a
  term of order k carries `F_k` and exactly k factors `S` -- so it is written
  in `Q_k = F_k / b^k`, each linear in `b = 1/sigma`, against `T_x = b S_x`;
  `F_4 = 6 y b^4 + 24 b^5` alone left the doubles at `alpha = 4.5e61`. All
  four orders are finite at `alpha = 1e300`.

  Measured on 2e5 points, the two changes together cost nothing: pig2 221.2
  ms against 222.4, pig1 214.9 against 226.9.

* The score of `pig2_distrib()` in `alpha` is closed for the same reason the
  value is. Since `r^2 - m^2 = alpha^2` the cancelling pair is
  `1/sigma - alpha = (r - m) - alpha`, and `r - alpha = m^2/(r + alpha)`, so
  its derivative is `-m^2/(r (r + alpha))`, a product of two bounded factors;
  the other piece, `-y asinh(m/alpha)`, gives `y m/(alpha r)`. Both are of
  size `alpha^-2`, which is what the score is. Assembled as
  `-y/alpha + psi'(alpha) + F_1 S_a` the three terms are each of size one:
  measured at `mu = 0.08`, the score was noise past `alpha = 1e8` and read
  exactly -1 past 1e162, where `1/alpha^2` leaves the doubles. `psi_derivs()`
  reports `A1 = psi'(alpha) + 1` beside `psi'`, that being a third
  cancellation. The closed form agrees with `numDeriv` to 7.8e-08 and tracks
  `-m^2/(2 alpha^2)` at a ratio of 1.000000 from `alpha = 1e4` to 1e100.

  It changes what a fit sees. On the model above the aggregate score over
  the province's 1344 cells changes sign between `alpha = 3.16` (+0.31) and
  `alpha = 10` (-0.73), so the likelihood asks for a finite dispersion; the
  direct form reported no slope at all there.

* `test-boundary-cancellation.R` gains the Poisson-inverse gaussian, the only
  case there whose cancellation is in the value rather than in a derivative.
  As the others it carries the negative control: the form that was replaced
  is more than 90 per cent out at every `alpha` from 1e16 up, so reverting to
  it fails the test rather than passing it silently.

# distributions7 0.44.0

* `approx = "opg"` computes the outer product of the observed scores, and is
  the default where a family has no closed-form expected information. It used
  to be an alias of `"bartlett"`: the two are readings of the same identity,
  `E[l_ij] = -E[l_i l_j]`, and differ by orders of magnitude in cost, because
  `"bartlett"` evaluates the expectation -- a sum over the support for a
  discrete family, a quadrature for a continuous one -- while `"opg"` reads
  the integrand at the observation.

  Measured on `pig1_distrib()` with `mu` varying by observation, one call at
  n = 1000: the identity evaluated 1.67 million rows of the family's
  derivative kernel and took 18.4 s, against n rows and 0.04 s. At n = 500 a
  `statmod()` regression goes from 89.06 s to 0.64 s, same log-likelihood,
  coefficients agreeing to 1.2e-06. Six families reach this route: pig1,
  pig2, pseudohuber, skewnormal1, skewnormal2 and skewt.

* `expected_by_opg()` is the new function; order 2 only, a higher order
  routing to `expected_by_bartlett()`, the outer product of scores being the
  second-order identity with no counterpart above it.

* `expected_hessian_exact()` is EXPORTED. It says whether a family writes its
  expected information out, which is what a consumer reads to decide which
  matrix to report; `statmodels7`'s `vcov()` now reads it.

* `fit_distrib()` reads the standard errors off the expected information only
  where the family writes it out, and off the observed Hessian otherwise --
  the condition is about the family and no longer about the optimizer that
  ran. A scoring step may be driven by any positive definite matrix, the
  score being exact, so the step takes the cheap one; a standard error is a
  different question. Measured on a PIG regression at n = 500, the outer
  product reports standard errors 5.7 per cent from those of the exact
  expectation where the observed information reports 0.6 per cent.

* `vcov()` on a `distrib_fit` takes `information` (`"fit"`, `"observed"` or
  `"expected"`) and `approx`, and recomputes at the optimum for the last two.

# distributions7 0.43.0

* `mvstudent_t_distrib()` is split into two numbered families, as
  `mvgaussian_distrib()` was at 0.42.0: `mvstudent_t1_distrib()` parametrizes
  the **scale matrix** and `mvstudent_t2_distrib()` its **inverse**. There is
  no alias. `MvStudentTDistrib` gains an `inverted` property and two concrete
  subclasses, and keeps every method.

  Unlike the gaussian this is a capability the family did not have: the
  inverse parametrization is new, not renamed.

* The constructor's page says what the two matrices ARE, because neither is a
  moment. The scalings are written out,

      Var(Y)    = nu / (nu - 2)  Sigma
      Var(Y)^-1 = (nu - 2) / nu  Sigma^-1

  so `mvstudent_t2_distrib()`'s matrix is the precision of the response only
  up to (nu - 2) / nu, which is 2/3 at nu = 6, and a coordinate named
  `omega_log_L1` builds Sigma^-1 and not the precision. The prefixes stay the
  gaussian's -- they say which of the two matrices a free value builds -- and
  the page says what those matrices mean.

  Two readings carry over to the response with no factor, both measured
  against the response's own covariance in the tests: the **correlations** of
  Sigma are the response's, and the **partial correlations** read off
  Sigma^-1 are the response's, a positive multiple cancelling out of a ratio.
  One does not: for a gaussian 1/Omega_jj is a conditional variance, here it
  is the Schur complement of the SCALE matrix. `mv_summary()` prints it as
  `cscale_vj`, by the rule that already makes the diagonal quantities
  `scale_sd_vj`.

* `mv_matrix_pieces()` is the one place a parametrization carrying the other
  side is transported: the matrix, its inverse, the log-determinant and the
  derivative arrays, all of Sigma. Both `mvg_pieces()` and `mvt_pieces()` read
  it, so the two families share the arithmetic instead of a copy of it.

  ⚠️ The log-determinant's derivatives are transported there too, and that is
  what the change was for. The sign flip used to be written by hand at each
  call site, and the Student t's new gradient was **14 relative** out while
  its Hessian was right, because the gradient read
  `parameters7::param_dlogdet()` directly. `check_distrib()` caught it -- the
  gradient and the score-mean-zero checks failed while the Hessian passed,
  which is the signature of exactly that term. Measured before assuming: the
  gaussian's own orders three and four, which also read the log-determinant
  directly, are CORRECT for the inverted form (2.66e-10 against a control of
  1.03e-10), so nothing there needed repair.

* `mv_marginal()` and `distrib_start()` keep the parametrization they are
  given, for the t as for the gaussian.

# distributions7 0.42.0

* `mvgaussian_distrib()` is split into two numbered families:
  `mvgaussian1_distrib()` parametrizes the covariance and
  `mvgaussian2_distrib()` the precision. The two arguments `sigma` and `omega`
  are gone with it, each constructor taking only its own. There is no alias:
  the old name is removed, as `covstructs7` was when it became `parameters7`.

  It is the convention the package already applies to `gaussian1`, `gaussian2`
  and `gaussian3` -- one name per parametrization, because a parametrization
  decides what a linear predictor acts on -- and the book had been arguing for
  weeks that the two are different families while the package gave one name.

  The split is not cosmetic, and the measurement is the reason. Where the
  matrix parametrization is closed under inversion the two describe the same
  laws in different coordinates: on 400 four-dimensional observations with
  `log_cholesky(4)` both reach -1867.271395 and differ by 5e-13. Where it is
  not, they are different models: the same measurement with `ar1(4)` gives
  -1868.65 against -1982.61, a gap of **113.95**, the inverse of an AR(1)
  covariance being tridiagonal and not an AR(1) at any parameters. A test
  pins both readings.

  `MvGaussian1Distrib` and `MvGaussian2Distrib` are the two classes, both
  subclasses of `MvGaussianDistrib`, which keeps every method and the
  `inverted` property: the two differ in what their free values describe, not
  in how anything is computed.

* Neither family comments on the parametrization it is given.
  `mvgaussian2_distrib(p, parameters7::ar1(p))` is the model whose precision
  has the entries rho^|i-j|; it is a legitimate law and is built in silence.
  To write the AR(1) process on the precision side, pass the family whose
  value is that inverse, `parameters7::ar1_inv(p)`.

* `mv_marginal()` keeps the parametrization it came from: the marginal of a
  precision-parametrized gaussian reports a precision, its free values read
  off the inverse of the covariance submatrix.

* Nothing reads `parameters7`'s `role`, which that package removes at 0.19.0.
  A census taken before the change found zero uses of it here; the side has
  always been recorded by `inverted` and reported through
  `params_interpretation`.

# distributions7 0.41.0

* Four moment methods read `theta` by position and now read it by name. Every
  generic here takes `theta` as a **named** list, and `align_theta()` reorders
  it, strips stray names off the values and validates it against the open
  `params_bounds`. No moment generic aligns before it dispatches, so each
  method has to, and 165 of 168 did.

* The consequence was a silently different number:

  ```r
  mean(enet_distrib(), list(mu = 0.3, lambda = 2, alpha = 0.7))
  #> 0.3
  mean(enet_distrib(), list(alpha = 0.7, lambda = 2, mu = 0.3))
  #> 0.7          the same three values, named in another order
  ```

  and on the von Mises in its resultant length the mean read 0.2594 one way
  and 0.3862 the other. Nothing warned, and every other generic in the package
  accepts any order, so there was no reason for a caller to write them in the
  family's own.

* The four are `mean`, `variance` and `skewness` of the elastic net, which
  reach `theta` through `.enet_parts()`, and `mean` of `vonmises2`, which reads
  `theta[[1]]` and calls `vm2_parts()`. `std_dev` of the elastic net is
  `sqrt(variance(...))` and is fixed by `variance`: it read 0.68749872 in the
  family's order and 0.69481641 reversed.

* The methods that delegate to `moment()` were already order-independent and
  are untouched. `moment()` does not align either, but it reaches its value
  through `expectation()` and `distrib_pdf()`, whose generic aligns before it
  dispatches.

* `.enet_parts()` and `vm2_parts()` are not the place for the fix: they receive
  `theta` and not the distribution, and `align_theta()` needs both.

* Aligning also **validates**, which is the behavior the other 165 have: a
  component outside its open domain, or a missing one, is now an error where
  these four returned a number. At a `theta` already in the family's own order
  nothing moves, aligning an ordered list being the identity.

* Found by three instruments that each miss what the others see, which is why
  all three are kept: a probe over the four moment generics, a probe over every
  generic, and a static scan of the method table for a body that indexes
  `theta` without aligning it. The first reported three, the second four, the
  third one.

# distributions7 0.40.0

* Every method of every generic accepts `...` where its generic declares one.
  0.29.0 did this for `distrib_pdf`, for a measured reason: the von Mises
  family needed a thread count carried down to `log_bessel_i`, and 43 of 45
  methods signalled `unused argument` on an argument they do not read. Eight
  generics were still split, and 197 of their methods refused.

* The argument for closing it is uniformity rather than extensibility. 27 of
  the package's 35 generics were already uniform, every method accepting; the
  eight that were not accepted on some families and refused on others, so the
  same call succeeded on one family and failed on another with no way for a
  reader to predict which. `distrib_rng` accepted on 5 methods of 56 and
  `distrib_quantile` on 11 of 44, and the ones that accepted were mostly the
  wrappers and the reparametrized families, which were forced to by having to
  forward to a parent.

* 171 inline method definitions and nine shared helper bodies are edited.
  `loc_scale_grad_cdf` serves five families and `trunc_cdf` two, so nine edits
  close 26 registrations; the helpers are named functions, so roxygen
  generates a `\usage` for them and the `@param ...` is required there rather
  than conventional.

* The cost, stated: an argument a method does not read is now absorbed rather
  than reported, so `distrib_cdf(d, q, theta, lower.tial = FALSE)` is silently
  ignored where it used to stop. That is the cost `distrib_pdf` has paid since
  0.29.0 and every one of the other 26 uniform generics pays; what changes is
  that it is now paid consistently.

* `mv_location` and `mv_sigma` declare no `...` at all, so their methods
  cannot be given one without changing the generic. They are left as they are.

* No behavior beyond the absorption changes, and the suite is unmoved at 6430
  passing.

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

  the two Bernoulli expansions canceling term by term at \eqn{\nu^{-2}}.
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

  `f1` splits into two canceling pairs rather than one, the second being
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
  profiled at **80.8 per cent** of a fit whose concentration is modeled, the
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
  information is a written-out formula -- the license the skew t already has
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
  canceling -- and the second ones grow like `gamma1^(-2/3)`, which is a
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
  first order -- the license the toolkit uses for an order it cannot check
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
  responses centered in the hundreds. `statmod()` on the Nile agrees with
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
