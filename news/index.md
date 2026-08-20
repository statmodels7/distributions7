# Changelog

## distributions7 0.27.4

- The cross-count twins of the parallel kernels compare at a tolerance
  of 1e-13 instead of
  [`identical()`](https://rdrr.io/r/base/identical.html). The Windows CI
  runner’s R runtime returns per-thread last bits from its own polygamma
  path: one ulp, deterministic, the same value at the same index across
  three independent binaries – with the sequential branch routed through
  the worker’s own function, and again with that function noinline –
  which is the opposite of a race signature and not something package
  code can bind. The decomposition guarantee (no reduction is ever
  split) stands, and the tolerance still fails a split reduction or a
  data race by ten orders.

## distributions7 0.27.3

- The worker’s loop in `d7::par_for()` is marked noinline, so the
  sequential branch and the parallel one execute the single compiled
  copy they both call. Routing the sequential branch through the same
  source function (0.27.2) had not been enough: the compiler inlined it
  at each call site and optimized the two copies apart, and the Windows
  CI runner went on producing one-ulp differences between one and two
  threads in the negbin kernels. Bit-identity across counts has to be a
  property of the binary, not of the source.

## distributions7 0.27.2

- `d7::par_for()`’s sequential branch runs through the worker over the
  whole range instead of writing a loop of its own, so both branches
  execute the same compiled function: the bit-identity across thread
  counts becomes a property of the code rather than of the optimizer.
  The Windows CI runner had reported last-bit differences between one
  and two threads in the two negbin kernels whose per-element arithmetic
  wraps `R::psigamma` calls, which no reordering of the decomposition
  can explain and which does not reproduce on this machine’s compiler.

## distributions7 0.27.1

- The negative binomial’s expected series helpers (`nb_E_trigamma`,
  `nb_E_psigamma`) no longer call `qnbinom` or `dnbinom` inside the
  parallel bodies: the quantile’s search reaches `pbeta`, whose warning
  path calls into the R API, and a warning raised from a worker thread
  trips R’s C-stack check and killed the test process on four of the
  five CI platforms. The series now stops on its own accumulated mass at
  the point the quantile located, with the underflowing head carried in
  log scale; the switch back to the multiplicative recurrence waits
  until the mass is comfortably normal, since seeding it at a subnormal
  was measured to carry a 2.5x error to the mode at `mu = theta = 1e4`.
  Values agree with the previous sizing to 1e-10 or better across nine
  regimes, and the rule in `d7_par.h` now names Rmath’s p/q functions as
  off limits inside a worker.

## distributions7 0.27.0

- The per-observation derivative kernels of the transcendental compiled
  families – poisson, negbin2, negbin1, beta1, student_t1, invgauss1,
  betabinom1 and gamma2, forty-five kernels over fourteen files – run
  their loops through the same d7 driver gaussian1 and gamma1 already
  use, at the transcendental threshold (128), with `threads` on their
  derivative methods. The loop bodies are untouched, so the results are
  bit-identical at any count, which the suite asserts kernel by kernel
  with [`identical()`](https://rdrr.io/r/base/identical.html); every
  file was read before conversion for the shared-buffer shape that would
  have made a mechanical pass a data race (none carried one).

## distributions7 0.26.0

- Scalar C entry points for the fast route of a score-driven filter
  (piano_parallel.txt, section 2a), registered with
  `R_RegisterCCallable`: `d7_scalar_id` keyed by the family’s S7 class
  name (gaussian1 and gamma1; an unknown name answers -1 and the
  consumer keeps its R callbacks) and `d7_score_curv`, the score and the
  (k, k) second derivative of the log-density in one parameter at one
  observation on the parameter scale, mirroring the family’s own vector
  kernels expression by expression. A twin test holds them against
  [`distrib_gradient()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.md)
  and
  [`distrib_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.md)
  with [`identical()`](https://rdrr.io/r/base/identical.html). The
  remaining compiled families take the same few lines each when a
  measurement names them; the sixteen in vectorized R have no C body to
  point at and stay on the callbacks.

## distributions7 0.25.0

- `fit_distrib(threads = numericals7::n_threads())` accepts the
  toolkit’s thread policy. The count travels as an argument down to the
  family’s compiled per-observation kernels; the process-level
  RcppParallel setting is sized at the fit’s entry and restored on exit.
  At the default, `n_threads(1)`, the code takes exactly the sequential
  path.
- The per-observation derivative kernels of
  [`gaussian1_distrib()`](https://statmodels7.github.io/distributions7/reference/gaussian1_distrib.md)
  (all four orders) and
  [`gamma1_distrib()`](https://statmodels7.github.io/distributions7/reference/gamma1_distrib.md)
  (all four orders, observed and expected) run their loops through one
  RcppParallel driver, decomposed over the elements of the output:
  observation i’s derivatives are computed and written in full by one
  thread, so no reduction is split and the result is bit-identical at
  any thread count, which a test asserts with
  [`identical()`](https://rdrr.io/r/base/identical.html). Their
  derivative methods take a `threads` argument (default 1) through the
  generics’ dots. Below a measured internal threshold a kernel stays
  sequential whatever the count says; the remaining compiled families
  take the same one-line conversion when a measurement names them.

## distributions7 0.24.0

- [`distrib_dexpected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_dexpected_hessian.md)
  is the derivative of the expected information in the parameters, one
  component per pair `(a, b)` and differentiating parameter `c`. It
  exists for a marginal criterion whose penalized matrix carries the
  EXPECTED information: that matrix enters through its determinant, so
  its gradient asks for `dK/dbeta`, which is `-l'''` with the observed
  information and `-dE[l'']/deta` with the expected one – and the two
  are different objects, because differentiating an expectation moves
  the measure as well as the integrand. The missing piece,
  `E[l_ab l_c]`, is a mixed moment no Bartlett identity isolates: the
  third ties the SYMMETRIZED sum, not the single term. The components
  are symmetric in `(a, b)` and NOT in `c`, so they are keyed by
  [`dexpected_names()`](https://statmodels7.github.io/distributions7/reference/dexpected_names.md)
  rather than by the sorted triples
  [`deriv_names()`](https://statmodels7.github.io/distributions7/reference/deriv_names.md)
  uses at order three.

- The default method is ONE central difference of the family’s own
  expected information, which is a single stencil on an analytic
  quantity wherever that information is a written-out formula – the
  licence the skew t already has for its degrees of freedom, and not the
  nested differencing the package forbids. Validated against the
  gaussian’s hand-written components (7.2e-11) and against the identity
  `E[l_abc] + E[l_ab l_c]` computed by quadrature on a beta, every one
  of whose six components is non-zero (8.5e-10).

- It REFUSES where the expected information is itself approximated, and
  the reason is cost rather than accuracy: measured at 100 observations
  the six families that approximate it cost 1880 to 147300 ms against a
  median of 0.183 ms for the thirty-four that do not, so 2p of those
  calls per evaluation is not a slower route but an unusable one.

- [`has_exact_expected_hessian()`](https://statmodels7.github.io/distributions7/reference/has_exact_expected_hessian.md)
  follows the arithmetic instead of the owning class, through the new
  generic
  [`expected_hessian_exact()`](https://statmodels7.github.io/distributions7/reference/expected_hessian_exact.md).
  Reading the owner is not sufficient and two families prove it:
  [`pseudohuber_distrib()`](https://statmodels7.github.io/distributions7/reference/pseudohuber_distrib.md)
  registers a method that calls the numerical
  [`expected_derivative()`](https://statmodels7.github.io/distributions7/reference/expected_derivative.md)
  and then replaces the two components vanishing by symmetry, and
  [`skewnormal2_distrib()`](https://statmodels7.github.io/distributions7/reference/skewnormal2_distrib.md)
  registers the chain onto
  [`skewnormal1_distrib()`](https://statmodels7.github.io/distributions7/reference/skewnormal1_distrib.md),
  whose expected information is the base class’s quadrature – costing
  5220 ms where the parent it chains onto costs 2230. Both answered
  “written out” about a quadrature, and the consequences were live:
  [`fit_distrib()`](https://statmodels7.github.io/distributions7/reference/fit_distrib.md)
  rejected a legitimate `fisher_scoring(approx = )` on them with a
  message stating the family computes its expected information in closed
  form, and its standard-error branch entered a multi-second quadrature
  believing it a formula.

## distributions7 0.23.0

- The centered skew normal rejects its derivatives at ZERO SKEWNESS,
  naming the reason and the parametrization that has none. Its map to
  the direct parameters runs through the cube root of `gamma1`, whose
  derivative is unbounded there: the first derivatives of the
  log-density survive the limit – they approach a finite value from both
  sides, the map’s factor cancelling – and the second ones grow like
  `gamma1^(-2/3)`, which is a property of the CENTERED parametrization
  and not of the family. Until now the resulting `NA` reached a
  comparison several frames further on and the message named nothing.
  Patching the first order alone would have been worse: the generic that
  a marginal criterion reads is the second. ⚠️ The point matters because
  it is where a hyperparameter STARTS: the bounds on `gamma1` are
  symmetric, so the midpoint rule puts the starting value at exactly the
  one point the family has no derivatives at.

- [`distrib_cross2_y()`](https://statmodels7.github.io/distributions7/reference/distrib_cross2_y.md),
  [`distrib_grad_y_hess()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y_hess.md)
  and
  [`distrib_hess_y_hess()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y_hess.md)
  are closed form for the multivariate Student t, which was the last
  family a marginal criterion could not estimate a correlated prior’s
  matrix with. Unlike the gaussian’s, this family’s response Hessian
  DEPENDS on the observation, so the first and third return one matrix
  per row.

- All four orders are written on ONE set of pieces
  ([`mvt_dpieces()`](https://statmodels7.github.io/distributions7/reference/mvt_dpieces.md)):
  the response derivatives are `-c w` and `-c Sigma^-1 + 2d ww'`, so
  everything follows from the first and second derivatives of
  `s = nu + q`, of `w`, of `Sigma^-1` and of the scalars `c` and `d`
  that follow from `s`. The second derivatives all vanish except four,
  and three of those share one middle matrix
  `A_k Sigma^-1 A_l + A_l Sigma^-1 A_k - A_kl`.
  [`distrib_cross_y()`](https://statmodels7.github.io/distributions7/reference/distrib_cross_y.md)
  was rewritten to read the same pieces, so its existing checks validate
  them at first order – the licence the toolkit uses for an order it
  cannot check directly.

- Verified against ONE difference of the analytic quantity below each,
  at p = 2 and 3: 1e-10 to 2e-10 throughout. And the gaussian limit is
  reached AT THE RATE 1/nu – a factor of 1e4 in nu divides every gap by
  1e4, in all three derivatives and at both dimensions.

## distributions7 0.22.0

- [`distrib_cross_y()`](https://statmodels7.github.io/distributions7/reference/distrib_cross_y.md)
  is closed form for the multivariate Student t as well. The response
  gradient is `-c w`, so every component carries the derivative of the
  weight beside the gaussian term it multiplies, and the degrees of
  freedom contribute `-(q - p) w / (nu + q)^2`. Nothing here is
  obstructed: the log-density carries no distribution function, only
  `lgamma`, a logarithm and a quadratic form, each elementary in `nu`.
  Against numDeriv on the analytic response gradient, 1e-9 to 1e-10 at p
  = 2, 3, 4; and the whole block becomes the gaussian’s AT THE RATE 1/nu
  – 1.17e-03, 1.17e-05, 1.17e-07 at nu of 1e4, 1e6, 1e8 – which an
  arithmetic accident does not do.

- [`distrib_cross2_y()`](https://statmodels7.github.io/distributions7/reference/distrib_cross2_y.md),
  [`distrib_grad_y_hess()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y_hess.md)
  and
  [`distrib_hess_y_hess()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y_hess.md)
  are closed form for the multivariate gaussian, which is what a
  marginal criterion reads to estimate the covariance of a correlated
  random effect. The response Hessian is `-Sigma^-1`, so it does not
  depend on the observation and does not depend on the mean at all:
  every component of the first two involving a mean is exactly zero and
  the rest are one matrix rather than one per row. ⚠️ Each is checked
  against ONE difference of the analytic quantity below it, never two in
  a row. A nested reference reported gaps of 0.3 on correct code.

- [`mv_derived()`](https://statmodels7.github.io/distributions7/reference/mv_derived.md)
  on a fixed-parameter wrapper delegates to the family instead of
  falling to the base method, which reports the distinct entries of the
  covariance: a centered prior was being read on a scale its own family
  does not use. The Jacobian keeps the columns of the free parameters
  alone.

## distributions7 0.21.0

- [`fixed()`](https://statmodels7.github.io/distributions7/reference/fixed.md)
  accepts a MULTIVARIATE family, in a third wrapper class beside the
  continuous and discrete ones. Holding the mean components at zero
  leaves the matrix parameter alone, which is what a centered prior on a
  random effect is. Every method splices and delegates as the other two
  do, and the generics a multivariate family rejects by design – the
  distribution function, the quantile – are inherited unregistered and
  go on rejecting.
  [`has_mv_support()`](https://statmodels7.github.io/distributions7/reference/has_mv_support.md)
  and
  [`has_mv_grad_y()`](https://statmodels7.github.io/distributions7/reference/has_mv_grad_y.md)
  unwrap first, so a wrapper’s delegation does not turn a family’s
  refusal into a TRUE.

- [`distrib_cross_y()`](https://statmodels7.github.io/distributions7/reference/distrib_cross_y.md)
  is closed form for the multivariate gaussian, in both
  parametrizations: `Sigma^-1 e_j` for the mean and `Sigma^-1 A_k w` for
  the matrix, one n-by-p matrix per parameter. It was refused for want
  of a consumer that fixed its shape, and a penalty whose prior is this
  family is that consumer. Agreement with numDeriv on the analytic
  response gradient is 1e-9 to 1e-10 at p = 2, 3, 4.

- `distrib_expected_hessian.MvStudentTDistrib` has a help page again:
  its roxygen block had fused with the one above it, which had lost its
  function, so the page had never been written.

## distributions7 0.20.0

- The gaussian’s derivatives are written in `z = (y - mu)/sigma` and
  `1/sigma`, never in a positive power of the scale.

  Written as `(res^2 - sigma^2)/sigma^3`, the score loses its
  denominator to overflow before the ratio itself becomes
  unrepresentable: at `sigma = 8e102` it returned exactly 0 where the
  value is `-1/sigma`, which on the link scale is -1, and at `1e200` it
  returned `NaN`. Zero is what a stopping rule reads as stationarity, so
  a run that had wandered out there reported `converged = TRUE` at a
  point that is not a maximum. The second derivative carried `sigma^4`
  and failed from `1.2e77`, the fourth `sigma^6` and from `4.4e61`. The
  variance chart is rewritten in `1/v` for the same reason; the
  precision chart already had no such form.

  Every component now agrees with the expression it replaces wherever
  that one held, and with the algebra beyond it. What remains out of
  reach above `1.3e154` is the LINK-SCALE second order, and there the
  cause is the chain rule rather than the kernel: it forms `(h')^2`,
  which overflows, against a component of order `1/sigma^2`, which
  underflows. The result is `NaN` rather than a plausible number, and a
  test pins that.

- The beta-binomial’s mass stays a probability at any shapes.

  `lchoose(n, y) + lbeta(y + a, n - y + b) - lbeta(a, b)` evaluated left
  to right adds a term of order one to a beta function of order `1e23`,
  which annihilates it, and the subtraction then leaves exactly zero.
  Every mass came back as one, the support summed to 11 instead of 1,
  and the log-likelihood was 0 – which beats any real fit. Fitted to
  binomial data, which asks the shapes to run to infinity at a fixed
  ratio, one start in six landed there and won.

  The shifts are integers, so each log-gamma difference is an exact sum
  of logarithms and the mass follows from three of them without forming
  any quantity larger than `n log(alpha + beta)`. The route is taken
  when the ordinary one is no longer accurate, at a threshold derived
  from the size of the terms it cancels rather than chosen. The mass now
  sums to one to 1e-13 at shapes from 1 to 1e23 and reaches the binomial
  limit to 1e-13.

- [`fit_distrib()`](https://statmodels7.github.io/distributions7/reference/fit_distrib.md)
  rejects a run of a discrete family whose log-likelihood is positive.

  A product of probabilities cannot exceed one, so such a run has left
  the region where the mass function is computable rather than found a
  better fit. It is discarded before the comparison, not ranked below
  the others: the objective is what breaks ties, and a number that is
  not a likelihood wins every tie it is allowed to enter. A
  log-likelihood of exactly zero is left alone, being what a legitimate
  fit against a support boundary reports.

- Six more families carry a method-of-moments estimate, taking the count
  to 37 of 42: `skewnormal1` and `skewnormal2` from the first three
  moments, `weibull3` and `pig2` by carrying the estimate of the chart
  they reparametrize across its map, and `betabinom1` and `betabinom2`
  by reading the intra-class correlation off the variance. Each recovers
  its parameters to better than one per cent on 4e5 draws. The five
  without one – `skewt`, `pseudohuber`, `gengamma1`, `gengamma2`, `enet`
  – have no closed inversion of their moments and keep the
  interpretation route.

  A family carrying a fixed constant spells it in its name, as in
  `"beta-binomial [size=10]"`, so the bracketed part is now dropped
  before the name is used as a lookup key: without that no beta-binomial
  of any size could match its own entry.

## distributions7 0.19.0

- A univariate family starts from the DATA, not from a draw over its
  parameters’ domains.

  The base
  [`distrib_start()`](https://statmodels7.github.io/distributions7/reference/distrib_start.md)
  never looked at `y`. That is harmless while the response is of order
  one and fatal when it is not: on a response of mean 919 and standard
  deviation 169 the draws are of order one, the first step is taken
  where the residuals are hundreds of standard deviations wide, and the
  scale runs to the largest representable double. Measured on a
  gaussian,
  [`fit_distrib()`](https://statmodels7.github.io/distributions7/reference/fit_distrib.md)
  recovered N(5, 2) and N(50, 20) and FAILED on N(500, 200) – a
  threshold in the scale of the data, not in the family. Downstream,
  `statmod(y ~ t, gaussian1_distrib(), Nile)` returned an intercept of
  227.9 and a slope of exactly 0 where `lm` gives 1056.4 and -2.71.

  What makes one method serve every family is that they already declare
  `params_interpretation`. A parameter meaning a location starts at the
  sample median, one meaning a spread at the sample standard deviation
  or its square, one meaning degrees of freedom at the value the sample
  kurtosis implies, and anything else – a shape, a dispersion, a
  probability – keeps its draw, being of order one whatever the data. A
  family declaring nothing recognizable loses nothing. Values are
  clamped strictly inside their bounds, a sample median being able to
  land on a support boundary.

  Now recovered across five decades of scale, N(5,2) to N(50000,20000),
  and on gamma, lognormal, Student t, Weibull, Poisson and negative
  binomial responses centred in the hundreds. `statmod()` on the Nile
  agrees with `lm` to 4.3e-16.

- [`moment_estimates()`](https://statmodels7.github.io/distributions7/reference/moment_estimates.md)
  gives a family its own method-of-moments estimate, which
  [`distrib_start()`](https://statmodels7.github.io/distributions7/reference/distrib_start.md)
  returns as the first starting value and refines by maximum likelihood
  from there. The interpretation route above is what a family without
  one falls back to. Each entry is checked against 2e5 draws from a
  known parameter, so the inversion is pinned to the family’s own
  generator rather than to the algebra it was written from – which is
  what catches a variance function transcribed from the wrong
  parametrization.

## distributions7 0.18.0

- `distrib_kernel(distrib, param)` returns the log-density, the score
  and the curvature in one parameter’s unconstrained scale, with the
  family’s methods and the link’s resolved once and the chain rule
  applied to the single component wanted rather than to all of them.

  The generic route is right almost everywhere: it validates its
  arguments, aligns `theta` by name, dispatches, and assembles every
  component of the order asked for. A recursion that calls back once per
  observation can afford none of that, and a score-driven filter
  evaluates its score at a predictor it has just produced, so the call
  cannot be vectorized away. Measured on one observation, the kernel is
  5.2x the generic for a gaussian score and 5.8x for its curvature; on a
  family whose own method is expensive the share is smaller (1.7x for
  the skew t).

  The caller takes on what the generic was doing: `theta` must already
  be in the family’s order with values of an acceptable length, and
  nothing is checked. The inverse link is still clamped strictly inside
  its bounds, which is a correctness property and not an optimization –
  an unclamped `exp(-800)` is zero, and a gaussian with a scale of
  exactly zero is not a distribution.

## distributions7 0.17.0

- The link-scale assembly no longer rebuilds its index layout on every
  call.
  [`to_link_scale()`](https://statmodels7.github.io/distributions7/reference/to_link_scale.md)
  computed, for each component of the requested order, the multi-index
  list, a `unique` and a `tabulate`, and then a `sort` and a `paste` per
  term of the nested sum to spell the key of the parameter-scale
  component it needed. None of that depends on the values; it is a
  function of the parameter names and the order alone, and
  [`link_scale_layout()`](https://statmodels7.github.io/distributions7/reference/link_scale_layout.md)
  now computes it once and caches it.

  It was found by profiling a fitted score-driven model, where `paste`,
  `sort` and `unique` were among the leaders of the self time: a filter
  reaches this once per observation per iteration, so the names were
  being respelled millions of times in one fit. Measured end to end, a
  gas(1,1) fit at 2000 observations went from 109 to 85 seconds, and an
  ordinary gaussian fit at 100000 observations from 1.00 to 0.89. The
  components are bit-for-bit what they were, checked on four families at
  all four orders.

## distributions7 0.16.0

- [`fit_distrib()`](https://statmodels7.github.io/distributions7/reference/fit_distrib.md)
  reports the observed information whenever the expected one has no
  closed form and the fit did not compute it.

  Which matrix the standard errors came from was decided by comparing
  `method` with the string `"newton"`. That argument has accepted an
  optimizers7 OBJECT since the delegation to that package, and an object
  is normalized internally to `"custom"`, so
  [`newton()`](https://statmodels7.github.io/optimizers7/reference/newton.html)
  failed the test and the expected information was assembled anyway – by
  quadrature, for a family that does not write the expectation out. On a
  user-defined Gompertz, whose score grows like `exp(by)` while the
  density decays like `exp(-eta exp(by))`, that quadrature does not
  converge: `method = "newton"` fitted in 0.15 seconds and
  `method = newton()`, the same algorithm on the same data, had not
  returned after five minutes. The surrounding `tryCatch` does not help,
  a quadrature that fails to converge raising nothing.

  The condition now asks what is actually wanted: the expected
  information when the fit used it (Fisher scoring) or when the family
  writes it out and it costs one evaluation, and the observed one
  otherwise. Families that carry a closed-form expectation are
  unaffected, to the digit.

## distributions7 0.15.0

- [`enet_distrib()`](https://statmodels7.github.io/distributions7/reference/enet_distrib.md)
  has a closed expected information, where it took one from the sampling
  fallback. Every piece was already in the family.

  In the two rates the density is an exponential family with sufficient
  statistics `-|z|` and `-z^2/2`, so `log Z` is its cumulant generating
  function and the information in `(a, c)` is exactly the Hessian of
  `log Z` – which `.enet_logz_derivs()` computes for the observed
  Hessian already. The map to `(lambda, alpha)` is bilinear, so the
  information transforms by `J' I J` with no second-derivative term.

  The location is where the expected and the observed part company, and
  the reason is the kink. The observed second derivative in `mu` is
  `-c`, which misses the point mass `d sgn(z)/dz = 2 delta(z)` the
  density carries at its own location, exactly as the Laplace does; the
  information there is the variance of the score,

  ``` R
  I_mu_mu = a^2 + 2ac E|z| + c^2 E[z^2]
          = a^2 - 2ac dlogZ/da - 2c^2 dlogZ/dc,
  ```

  since `E|z| = -dlogZ/da` and `E[z^2] = -2 dlogZ/dc`. It reduces to
  `lambda^2` at `alpha = 1`, which is the Laplace’s, and to `c` at
  `alpha = 0`, which is the Gaussian’s, and both are asserted against
  the families that own them.

  Validated against the definition – minus the outer product of the
  score, integrated against the density and SPLIT AT THE KINK, since a
  rule that straddles the corner measures the corner: 2.3e-14 to 1.1e-6
  over five settings, the worse ones being the quadrature and not the
  formula. A Monte Carlo of 4e6 draws agrees within its own standard
  error at all of them.

- The multivariate validator’s own gradient and second-derivative
  stencils take their nodes, weights and step from numericals7 rather
  than writing out a central difference with a step of their own. The
  rules are identical, so the checks report what they reported.

## distributions7 0.14.0

- The mixed grid reaches every shape that needs no new algebra. A family
  written out in its own parametrization with its own kernels –
  gaussian2, gaussian3, laplace2, invgauss2 – is not a
  [`reparametrize()`](https://statmodels7.github.io/distributions7/reference/reparametrize.md)
  wrapper and does not inherit its methods, but its map onto a parent is
  already tabulated in `reparam_maps.R`, so the same chain rule applies
  and only the registration was missing.

- A family with no location keeps every scale formula: the derivation
  never used the location, only that sigma is a scale, so with
  `z = y/sigma` they hold as written. That covers the exponential, the
  Weibull’s scale and the generalized Pareto’s, and the exponential is
  closed outright, having nothing but a scale.

- The lognormal closes by the transformation: it is a gaussian at
  `t = log y`, and the transformation carries no parameter, so `t` does
  not move with theta and every theta-derivative is the gaussian’s own
  at that point. What the response derivatives carry is the Jacobian.
  lognormal2 follows through its map, which is what closing a parent is
  for.

- [`distrib_cross2_y()`](https://statmodels7.github.io/distributions7/reference/distrib_cross2_y.md),
  [`distrib_grad_y_hess()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y_hess.md)
  and
  [`distrib_hess_y_hess()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y_hess.md)
  reach 21 families each, up from 13. The eleven that remain are not
  location-scale and have no transformation onto one.

## distributions7 0.12.0

- A reparametrized family carries the whole mixed grid through its map:
  [`distrib_cross2_y()`](https://statmodels7.github.io/distributions7/reference/distrib_cross2_y.md)
  by the same first-order chain rule
  [`distrib_cross_y()`](https://statmodels7.github.io/distributions7/reference/distrib_cross_y.md)
  takes, and the two second-order generics by the ordinary two-term
  expansion, which needs the map’s second partials and the parent’s
  first-order components.
  [`reparam_tables()`](https://statmodels7.github.io/distributions7/reference/reparam_tables.md)
  already keys both, so nothing new is differentiated.

- The laplace was registered for
  [`distrib_cross2_y()`](https://statmodels7.github.io/distributions7/reference/distrib_cross2_y.md)
  and not for the two second-order generics. A census found it: a family
  half-registered answers at one order and falls back at the next
  without anything failing.

- The three generics now reach 13 families each, up from 9, 8 and 8.

## distributions7 0.11.0

- The location-scale identity closes the new mixed generics for the
  families it applies to. Where the response enters only through
  `z = (y - mu)/sigma`, every derivative of the response gradient and of
  the response curvature in the location and the scale is that family’s
  own y-derivatives times a power of sigma, so nothing new is derived:
  [`distrib_cross2_y()`](https://statmodels7.github.io/distributions7/reference/distrib_cross2_y.md),
  [`distrib_grad_y_hess()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y_hess.md)
  and
  [`distrib_hess_y_hess()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y_hess.md)
  go from 2, 1 and 1 closed families to 9, 8 and
  8.  A family with a shape parameter beyond the two keeps its
      location-scale pairs closed and differences the rest, as the
      first-order block already does.

## distributions7 0.10.0

- [`distrib_grad_y_hess()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y_hess.md)
  and
  [`distrib_hess_y_hess()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y_hess.md)
  close the mixed grid: one or two derivatives in the response and TWO
  in the parameters. They are what the SECOND derivative of a marginal
  criterion needs of a penalty. Closed forms for the gaussian, one
  central difference of the analytic first-order component otherwise,
  and
  [`fixed()`](https://statmodels7.github.io/distributions7/reference/fixed.md)
  delegates them like every other derivative.

## distributions7 0.9.1

- [`fixed()`](https://statmodels7.github.io/distributions7/reference/fixed.md)
  delegates
  [`distrib_cross2_y()`](https://statmodels7.github.io/distributions7/reference/distrib_cross2_y.md),
  subset to the free parameters as it does every other derivative.
  Without it every penalty built on `fixed(gaussian1_distrib(), mu = 0)`
  – which is what a ridge and a random effect are – reached the
  numerical fallback while the closed form sat one class away.

## distributions7 0.9.0

- [`distrib_cross2_y()`](https://statmodels7.github.io/distributions7/reference/distrib_cross2_y.md)
  completes the mixed surface: two derivatives in the response and one
  in each parameter, which is how the curvature of a log-density in the
  response moves with the parameters. A penalty is a negative
  log-density evaluated at the coefficients, so this is what a marginal
  likelihood needs to differentiate the determinant of a penalized
  information. Closed forms for the gaussian and the Student t, and one
  central difference of the analytic
  [`distrib_hess_y()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y.md)
  for everything else.

## distributions7 0.8.0

### Derivatives of the distribution function

- [`skewnormal1_distrib()`](https://statmodels7.github.io/distributions7/reference/skewnormal1_distrib.md)
  is closed in the shape as well as in the location and the scale, at
  every order, and
  [`skewnormal2_distrib()`](https://statmodels7.github.io/distributions7/reference/skewnormal2_distrib.md)
  follows through the map. The distribution function is and Owen’s has
  elementary partial derivatives, and , so the integral in its
  definition is differentiated away at the first order and never has to
  be differentiated again. Everything above is a product of normal
  densities, Hermite polynomials and a rational function of the shape.

- The check that the first identity is the right one is that it returns
  the density: . At the location and scale components agree with the
  gaussian’s to 1e-15, and the gaussian reaches them by another route.

- Every cdf order now leaves the same eight families on the stencil, and
  all eight are mathematical obstructions: the derivative of an
  incomplete gamma or beta in its shape is hypergeometric (gamma,
  chi-squared, beta, generalized gamma) and the von Mises distribution
  function is itself a quadrature. A test asserts that list, so a family
  added without a route joins it visibly.

## distributions7 0.7.0

### Derivatives of the distribution function

- [`invgauss1_distrib()`](https://statmodels7.github.io/distributions7/reference/invgauss1_distrib.md)
  and
  [`enet_distrib()`](https://statmodels7.github.io/distributions7/reference/enet_distrib.md)
  have closed derivatives at all four orders, from one route: a
  distribution function of the form is a Leibniz split between the
  weight and the tail, with a Faa di Bruno on each side. A family
  supplies, per term, the partial derivatives of the log weight and
  those of the argument.

  The inverse gaussian’s has all three of , and separable in the mean
  and the dispersion, so their mixed partials are products of
  one-variable ones. The elastic net’s halves are truncated Gaussians,
  and its and are likewise separable in and ; its weight is written
  through the Mills ratio the family already carries, , so its
  derivatives come from the same the density uses.

- The weight is never formed on its own. is at ordinary settings – 2000
  in the exponent at , – exactly where underflows, so the two are
  combined as and the fourth derivative comes back finite.

- [`invgauss2_distrib()`](https://statmodels7.github.io/distributions7/reference/invgauss2_distrib.md)
  takes its Hessian through the mapped route too. Registering only the
  gradient there was right while the parent differenced its own second
  order and is not now that it does not.

- Every one of the four cdf surfaces now leaves the same nine families
  on the stencil, and all nine are obstructions or correct refusals: the
  derivative of an incomplete gamma or beta in its shape is
  hypergeometric (gamma, chi-squared, beta, generalized gamma), the von
  Mises distribution function is itself a quadrature, and
  [`skewnormal2_distrib()`](https://statmodels7.github.io/distributions7/reference/skewnormal2_distrib.md)
  is refused by the gate while its parent’s shape components are
  differenced.

## distributions7 0.6.0

### Derivatives of the distribution function

- [`gpd_distrib()`](https://statmodels7.github.io/distributions7/reference/gpd_distrib.md)
  has closed derivatives at all four orders, its Hessian included, from
  the exponential-survival route. The form is what makes it work:
  writing and , the exponent is and carries no division by the shape, so
  the exponential limit is an ordinary point of the formula rather than
  a branch of the code.

- `Lambda` and its four derivatives come from the recursion above and
  from the Taylor series below it. The crossover is measured, not
  chosen: the recursion divides by and subtracts nearly equal
  quantities, and its fourth derivative is wrong by a factor of at , by
  1.7 at and by at , while the two agree to at the switch.

- At the scale components equal the exponential family’s exactly, at
  every order. That family reaches them through and shares no arithmetic
  with the series, so it is the reference near zero, where a stencil in
  the shape is not one: at the step is a thousand times the value.

- A negative shape bounds the support above, at , and the derivatives
  are zero past it. The support is declared by the family rather than
  read off the bounds, the endpoint depending on a parameter.

## distributions7 0.5.0

### Derivatives of the distribution function

- A family whose survival function is an exponential of something
  elementary now gets all four orders from one identity: gives , the
  complete Bell polynomial in the partial derivatives of . A family
  states and its partials and nothing else.
  [`exponential_distrib()`](https://statmodels7.github.io/distributions7/reference/exponential_distrib.md)
  and
  [`weibull1_distrib()`](https://statmodels7.github.io/distributions7/reference/weibull1_distrib.md)
  are served by it, and
  [`weibull3_distrib()`](https://statmodels7.github.io/distributions7/reference/weibull3_distrib.md)
  follows through the reparametrization wrapper; the route replaces the
  two written-out orders each carried before, which it reproduces.

- [`laplace2_distrib()`](https://statmodels7.github.io/distributions7/reference/laplace2_distrib.md)
  reaches the two new orders as the Laplace at , through the mapped
  route.

- The upper tail of an exponential-survival family is exact wherever the
  logarithm is representable. `log S` is `L`, so its derivatives are
  `L`’s own and need no division by `1 - F`, which is exactly one in
  double precision past `q/mu = 37`: the first derivative of an
  exponential’s log survival at `q = 700` was `Inf` and is now `700`.

- `weibull_cdf_deriv()` is gone with the two methods that were its only
  callers.

## distributions7 0.4.0

### Derivatives of the distribution function

- The third and fourth derivatives are closed for thirteen more
  families, taking the count from four to seventeen. Three routes did
  it, and none derived anything new.

  A family written as a map of another carries the parent’s through one
  Faa di Bruno pass; orders one and two already did this and the new
  orders use
  [`chain_assemble()`](https://statmodels7.github.io/distributions7/reference/chain_assemble.md),
  the enumeration the reparametrized parameter derivatives run on, so no
  second copy of the partition sum exists. That closes
  [`gaussian2_distrib()`](https://statmodels7.github.io/distributions7/reference/gaussian2_distrib.md),
  [`gaussian3_distrib()`](https://statmodels7.github.io/distributions7/reference/gaussian3_distrib.md)
  and every
  [`reparametrize()`](https://statmodels7.github.io/distributions7/reference/reparametrize.md)
  wrapper whose parent is exact.

  The mapped route now admits a transformation of the response as well
  as of the parameters. A lognormal is a gaussian at `log q` and the
  transformation carries no parameter, so the derivatives in the
  parameters are the gaussian’s with the point substituted; that closes
  [`lognormal1_distrib()`](https://statmodels7.github.io/distributions7/reference/lognormal1_distrib.md),
  and
  [`lognormal2_distrib()`](https://statmodels7.github.io/distributions7/reference/lognormal2_distrib.md)
  follows from it through the wrapper.

  [`gumbel_distrib()`](https://statmodels7.github.io/distributions7/reference/gumbel_distrib.md)
  joins the location-scale families, and
  [`student_t1_distrib()`](https://statmodels7.github.io/distributions7/reference/student_t1_distrib.md),
  [`pseudohuber_distrib()`](https://statmodels7.github.io/distributions7/reference/pseudohuber_distrib.md),
  [`skewnormal1_distrib()`](https://statmodels7.github.io/distributions7/reference/skewnormal1_distrib.md)
  and
  [`skewt_distrib()`](https://statmodels7.github.io/distributions7/reference/skewt_distrib.md)
  get the location and scale components from that construction with the
  shape components still differenced, as at the two orders below.

- The gate is the one orders one and two use: a chain rule is taken only
  when the parent is exact at every order up to the one asked for, so a
  differenced quantity is never reported as a closed form. Measured
  against the partial-expectation integral, which shares no code with
  any of the three routes, the fully closed families agree to 4e-15 and
  the partial ones to the stencil’s own 3e-6.

## distributions7 0.3.0

### Derivatives

- The third and fourth response derivatives are closed for every
  continuous family. Eighteen were taking a finite-difference stencil,
  all of them families whose response is not a pure location; each
  already carried a closed second response derivative, and the third is
  the same elementary function differentiated once more. The log-density
  of each is a sum of terms in `log(y)`, `log(1 - y)`, a power of `y`, a
  logarithm of an affine function of `y`, or a cosine, so the terms are
  written once and each family is a sum of them.

- [`reparametrize()`](https://statmodels7.github.io/distributions7/reference/reparametrize.md)
  carries the third and fourth response derivatives to the parent, as it
  already carried the first and second: the map acts on the parameters
  and the derivative is taken in the response, so the two do not
  interact.

- The generalized Pareto’s coefficient is written as `xi^k + xi^(k-1)`
  rather than `(1 + 1/xi) * xi^k`. It is the same number and stays
  finite as the shape goes to zero, where the family is exponential and
  every order above the first is exactly zero.

  Measured against one differentiation of the analytic second response
  derivative, the third order agrees to 1e-11 and the fourth to 4e-5,
  which is each reference’s own accuracy. A test walks the namespace and
  fails if a continuous family is left on the stencil.

## distributions7 0.2.0

### Derivatives

- Every one of the 46 families is now analytic to fourth order in the
  parameters. The three that were still on the numerical fallback –
  [`dirichlet_distrib()`](https://statmodels7.github.io/distributions7/reference/dirichlet_distrib.md),
  [`multinomial_distrib()`](https://statmodels7.github.io/distributions7/reference/multinomial_distrib.md)
  and
  [`mvstudent_t_distrib()`](https://statmodels7.github.io/distributions7/reference/mvstudent_t_distrib.md)
  – have closed third and fourth derivatives. The two simplex-valued
  log-densities are a sum of terms each depending on one coordinate, so
  the chain rule is one univariate partition sum per coordinate; the
  Student t splits into the mean-and-matrix part, which reuses the
  gaussian’s expansion of the derivative of an inverse, and a part in
  the degrees of freedom, which is elementary. In each case the same
  assembly run at orders one and two reproduces the hand-written score
  and information, to 9e-16 and 4e-15.

- The mixed derivative
  [`distrib_cross_y()`](https://statmodels7.github.io/distributions7/reference/distrib_cross_y.md)
  of
  [`skewnormal1_distrib()`](https://statmodels7.github.io/distributions7/reference/skewnormal1_distrib.md)
  and
  [`pseudohuber_distrib()`](https://statmodels7.github.io/distributions7/reference/pseudohuber_distrib.md)
  is closed in the shape parameter as well as in the location and the
  scale. The response reaches the skew normal’s shape only through
  `alpha * z` and the pseudo-Huber’s only through `D`. The skew normal
  assembles all three components from a single evaluation of the inverse
  Mills ratio rather than through
  [`distrib_grad_y()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.md)
  and
  [`distrib_hess_y()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y.md),
  which evaluate it twice more: 77.5 ms to 19.4 at n = 1e5.

### Internals

- [`mvg_ptensors()`](https://statmodels7.github.io/distributions7/reference/mvg_ptensors.md)
  takes the pieces rather than the distribution, so the Student t reuses
  one copy of the expansion, and its accessor answers for the empty
  multiset.

## distributions7 0.1.0

### Plots

- [`plot()`](https://rdrr.io/r/graphics/plot.default.html) on a
  univariate distribution draws one curve per element of a parameter
  given as a vector, so
  `plot(gaussian1_distrib(), list(mu = 0, sigma = c(1, 2, 4)))` is three
  densities on one panel. The settings are separated by color and by
  line type together, which keeps them apart in a printed copy that has
  no color; the parameters that vary are named in a legend, placed on
  whichever side the mass leaves emptier, and those held fixed are
  stated once in the title. A discrete family is drawn as several sets
  of stems, shifted sideways so that equal masses at one support point
  stay countable, and separated by symbol rather than by line type: a
  dashed stem reads as a broken one, and at a support of any size the
  panel fills with fragments.

  Every component must have length one or the same `k`. A length that
  merely divides `k` is rejected rather than recycled, since a partial
  setting is far more likely to be a mistake than a request. The
  horizontal window covers every setting rather than the first.

  This meaning is available because a plot has no data to recycle
  against; the density and derivative generics read a vector component
  as one value per observation, which is a different question asked of
  the same object. A multivariate family, whose picture is already a
  matrix of panels with no axis left to overlay settings on, rejects a
  vector component instead.

### Families

- One name per parametrization. A family with several parametrizations
  carries a number on each – `gaussian1`/`gaussian2`/`gaussian3`,
  `gamma1`/`gamma2`, `negbin1`/`negbin2` after Cameron and Trivedi,
  `weibull1`/`weibull3` after gamlss’s WEI and WEI3 with `weibull2`
  deliberately empty, and ten further groups. The reference index lists
  what is present.

- Univariate families added: weibull, gumbel, skewnormal, skewt,
  exponential, geometric, chisq, betabinom, NB1, the generalized Pareto,
  the generalized gamma, von Mises, the Poisson-inverse Gaussian in both
  gamlss parametrizations, and
  [`enet_distrib()`](https://statmodels7.github.io/distributions7/reference/enet_distrib.md),
  the elastic-net prior.

- Multivariate families: the gaussian and the Student t, whose matrix
  parameter comes from `parameters7`, and the Dirichlet and the
  multinomial, whose simplex parameter does. The matrix parameter is
  flattened into scalars with identity links, so
  [`align_theta()`](https://statmodels7.github.io/distributions7/reference/align_theta.md),
  [`deriv_names()`](https://statmodels7.github.io/distributions7/reference/deriv_names.md),
  the link scale and
  [`fit_distrib()`](https://statmodels7.github.io/distributions7/reference/fit_distrib.md)
  need no special case; the constraint lives in the structure, where it
  belongs. The base class sits beside `continuous_distrib` and
  `discrete_distrib` rather than under either, the one-dimensional
  defaults registered there – a cdf by quadrature, a quantile by root
  finding – having no counterpart in several dimensions.

- [`mv_summary()`](https://statmodels7.github.io/distributions7/reference/mv_summary.md)
  reports the quantities a reader reads rather than the coordinates:
  standard deviations and correlations with delta-method standard
  errors, each interval built on the scale that keeps the quantity in
  its own set (log for a standard deviation, Fisher’s z for a
  correlation) and mapped back. A precision parametrization adds the
  conditional standard deviations and the partial correlations.

- Wrappers:
  [`fixed()`](https://statmodels7.github.io/distributions7/reference/fixed.md),
  holding parameters at known values and the only one that removes
  parameters;
  [`folded()`](https://statmodels7.github.io/distributions7/reference/folded.md),
  for the absolute value; and
  [`truncated()`](https://statmodels7.github.io/distributions7/reference/truncated.md),
  [`zero_inflated()`](https://statmodels7.github.io/distributions7/reference/zero_inflated.md),
  [`zero_adjusted()`](https://statmodels7.github.io/distributions7/reference/zero_adjusted.md)
  and
  [`transformation()`](https://statmodels7.github.io/distributions7/reference/transformation.md)
  from earlier. Neither zero wrapper can be stacked on the other –
  truncating at zero cancels `(1 - zeta)` between the numerator and the
  truncation constant, so `zeta` leaves the likelihood entirely – and
  both are rejected by the constructor, along with a discrete parent
  carrying too few support points to identify the extra probability.

### Derivatives

- Every univariate family is analytic to fourth order, observed, except
  the skew t’s components in `nu`, which cannot be: the density carries
  `T_{nu+1}` and the derivative of a Student t distribution function in
  its degrees of freedom has no elementary form. Those come from one
  five-point stencil applied to an analytic quantity, never from a
  difference of a difference.

- Orders three and four are closed form for every wrapper. Each
  wrapper’s log-likelihood is the parent’s log-density plus, or instead
  of, `log L` for some parameter-dependent `L`, so two partition sums
  cover all of them: the complete Bell polynomial, and the
  moment-to-cumulant relation, which needs only the ratios `d^B L / L`.

- [`distrib_deriv3_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3_cdf.md)
  and
  [`distrib_deriv4_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3_cdf.md)
  complete that sequence. The routes are the ones the orders below use:
  a discrete family sums the identity exactly and a continuous one
  applies one product stencil to its analytic distribution function.
  Nothing new was derived – the quantity summed is the complete Bell
  polynomial and the conversion to the log scale the moment-to-cumulant
  relation, both already in the package for the wrappers and now in
  `partition_sums.R`, where they belong. The general forms reproduce the
  written-out orders one and two exactly, on both tails and both scales.

- [`distrib_deriv3_y()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3_y.md)
  and
  [`distrib_deriv4_y()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3_y.md)
  do the same for the response. A family whose response enters only as
  `y - mu` needs no formula of its own:
  `d^k l / dy^k = (-1)^k d^k l / dmu^k`, so fourteen families inherit
  these orders from derivatives they already have, often compiled ones.
  The identity is checked at orders one and two, where both sides are
  written independently and it holds exactly. A family on a half line
  takes one stencil of the order asked for, with the step halved because
  the stencil reaches two steps either side of a support boundary.

- [`distrib_grad_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_cdf.md)
  and
  [`distrib_hess_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_cdf.md),
  the derivatives of the distribution function, which is what a censored
  likelihood and a quantile residual’s standard error need. Closed form
  for twelve of the original fourteen families; gamma and beta have
  none, the shape direction being hypergeometric.

- [`distrib_cross_y()`](https://statmodels7.github.io/distributions7/reference/distrib_cross_y.md),
  the mixed response-parameter derivatives, closed form for every
  continuous family. Where the response enters only through
  `z = (y - mu)/sigma` the identity `d2l/dy dmu = -l_yy` and
  `d2l/dy dsigma = -z l_yy - l_y/sigma` closes nine families at once.

- [`reparametrize()`](https://statmodels7.github.io/distributions7/reference/reparametrize.md),
  building a family from another through a map, with Faa di Bruno over
  partitions and the map’s partials as hand-written keyed tables
  (`map_derivs`); without tables, one stencil per partial. Measured
  against a family written out in full it costs 1.5x at the gradient and
  6.6x at order four, so it is the user-facing route and new families
  are written out.

- Jets are removed from every production path. Generic jet composition
  measured at 2x to 36x the hand-written closed forms on the PIG
  kernels; the mechanical transcriptions survive in the tests as
  independent references.

### Fitting

- [`fit_distrib()`](https://statmodels7.github.io/distributions7/reference/fit_distrib.md)
  delegates its optimization to `optimizers7`. Fisher scoring is
  [`newton()`](https://statmodels7.github.io/optimizers7/reference/newton.html)
  with the expected information passed as `he`, and `method` accepts an
  optimizer object, so nothing in this package implements a descent
  loop.

- `fisher_scoring(approx =, nsim =, criterion =, maxit =)` replaces the
  loose `approx` and `nsim` arguments: how the expected information is
  approximated is a property of Fisher scoring and had no business
  sitting beside optimizers that never look at it. A strategy chosen
  where it would be ignored is rejected.

- `maxit` and `tol` leave the signature. With `method = <an optimizer>`
  they were silently discarded, so a call setting both got no complaint
  and no effect from the second; the budget and the stopping rule now
  live on the method.

- The objective is `-l(eta)/n`. The maximizer and every Newton step are
  unchanged, the factor canceling in `H^-1 g`; what changes is what a
  threshold means, an absolute gradient tolerance on a summed score
  asking of a sample of ten million an accuracy per observation ten
  million times finer than of a sample of ten. `ll_hat` and the
  information are recomputed unscaled at the optimum, so
  [`logLik()`](https://rdrr.io/r/stats/logLik.html), AIC, BIC and every
  standard error are untouched.

- [`distrib_start()`](https://statmodels7.github.io/distributions7/reference/distrib_start.md)
  computes a starting value from the data. The default is the old random
  draw from the parameter domains; the multivariate gaussian returns the
  sample mean and covariance, its own maximum likelihood estimate for an
  unstructured matrix.

- The restart loop keeps the best result rather than the last: a
  converged run beats a non-converged one, and among runs of equal
  status the lower objective wins. A run that reached a point is no
  longer discarded as a failure.

- The default tolerance is `1e-6`, for the reason recorded in
  optimizers7’s own notes: the attainable gradient is bounded below by
  `sqrt(2*lambda*eps*|f*|)`, and a log-likelihood is of order one at its
  maximum.

- [`confint()`](https://rdrr.io/r/stats/confint.html), with
  `scale = c("parameter", "link")`, recomputing at any level from the
  stored estimate and standard error. The link-scale interval is the one
  computed and the parameter-scale table is its image under the inverse
  link, so [`print()`](https://rdrr.io/r/base/print.html) shows both.

- A fit records its elapsed time, accumulated over every start and every
  fallback, and the score per observation at the point it stopped.

### Numerical layer

- The enumerations, stencils, batched quadrature and series, and special
  functions come from `numericals7`.
  [`expectation()`](https://statmodels7.github.io/distributions7/reference/expectation.md)’s
  integrand contract is elementwise in the response and the parameters
  jointly, so every parameter combination shares one batched call.

- Random numbers come from generalized ratio-of-uniforms, recentered at
  the mode and normalized. A density diverging at an edge is transformed
  away, and the exponent is measured by the same probe that detects the
  divergence rather than searched for: Gamma shape 0.4 went from 27 ms a
  draw to 0.8 microseconds. Discrete families invert the cumulative
  table, which is exact.

### Validation and documentation

- [`check_distrib()`](https://statmodels7.github.io/distributions7/reference/check_distrib.md)
  runs thirteen numerical checks on a continuous family, twelve on a
  discrete one and a nine-check battery on a multivariate one. It is
  aware of atoms, so a mixed distribution such as
  `zero_adjusted(gamma1_distrib())` is not reported as four failures on
  correct code; it allows for a kink where `params_smooth` declares one;
  and a gradient made five per cent wrong is still caught.

- [`mv_reference_draw()`](https://statmodels7.github.io/distributions7/reference/mv_reference_draw.md)
  supplies the proposal the normalization check integrates against. The
  gaussian proposal does not fail loudly on the Dirichlet –
  [`chol()`](https://rdrr.io/r/base/chol.html) accepts the singular
  covariance – and returns 2.0e-08 for an integral that is 1.

- Every family’s constructor page displays its density, taken from the
  book’s catalog, whose transcription is checked against
  [`distrib_pdf()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.md)
  at every render.

- Three vignettes – defining a distribution, fitting a model,
  derivatives and the link scale – a README with badges, a pkgdown site
  and continuous integration on five platforms.
