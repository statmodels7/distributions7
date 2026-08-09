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
