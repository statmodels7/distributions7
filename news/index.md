# Changelog

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
