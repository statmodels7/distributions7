# S7 Class for Maximum-Likelihood Fits

The class of the object
[`fit_distrib()`](https://statmodels7.github.io/distributions7/reference/fit_distrib.md)
returns. It holds the estimates on both the parameter scale and the link
scale, their standard errors and confidence limits, the maximized
log-likelihood with the two information criteria built on it, and a
record of what the optimizer did.

Both scales are kept because the fit is computed on one and read on the
other. The variance matrix is the inverse information at \\\hat\eta\\;
`vcov` is its image under the delta method and `ci` is `ci_eta` mapped
through \\g^{-1}\\, so a limit on the parameter scale respects the
parameter's domain by construction.

This page documents the raw S7 constructor. It validates nothing and is
not the way to build a fit; call
[`fit_distrib()`](https://statmodels7.github.io/distributions7/reference/fit_distrib.md).

## Usage

``` r
distrib_fit(
  distrib = distrib(),
  y = integer(0),
  n = integer(0),
  coefficients = integer(0),
  se = integer(0),
  ci = NULL,
  eta = integer(0),
  se_eta = integer(0),
  ci_eta = NULL,
  vcov = NULL,
  vcov_eta = NULL,
  loglik = integer(0),
  aic = integer(0),
  bic = integer(0),
  iterations = integer(0),
  converged = logical(0),
  method = character(0),
  criterion = character(0),
  note = character(0),
  counts = NULL,
  score = integer(0),
  elapsed = integer(0),
  level = integer(0)
)
```

## Arguments

- distrib:

  The fitted `distrib` object, carrying the parametrization and the
  links the estimates are expressed in.

- y:

  The observations the fit was computed from, kept so that the fitted
  distribution can be compared with the data by
  [`plot.distrib_fit()`](https://statmodels7.github.io/distributions7/reference/plot.distrib_fit.md)
  and resampled by
  [`simulate.distrib_fit()`](https://statmodels7.github.io/distributions7/reference/simulate.distrib_fit.md).

- n:

  The number of observations: the row count for a multivariate response
  and the length otherwise. `bic` and the printed header both read it.

- coefficients:

  Named numeric estimates on the parameter scale, one per parameter, in
  the order `distrib@params` gives.

- se:

  Named standard errors on the parameter scale, from the delta method.
  `NaN` where the corresponding variance came back negative or missing.

- ci:

  A two-column numeric matrix of confidence limits on the parameter
  scale, one row per parameter, columns `lower` and `upper`.

- eta:

  Named estimates on the link scale, \\\hat\eta = g(\hat\theta)\\.

- se_eta:

  Named standard errors on the link scale, the square roots of the
  diagonal of `vcov_eta`. These are the ones the fit computes; `se` is
  derived from them.

- ci_eta:

  A two-column numeric matrix of confidence limits on the link scale,
  symmetric about `eta`.

- vcov:

  The variance-covariance matrix on the parameter scale, with both
  dimnames set to the parameter names.

- vcov_eta:

  The variance-covariance matrix on the link scale, the inverse of the
  information at \\\hat\eta\\. Every entry is `NA` when the information
  could not be evaluated or inverted there; the estimates stand in that
  case and only the uncertainty is missing.

- loglik:

  The maximized log-likelihood, summed over observations and **not**
  divided by \\n\\.

- bic, aic:

  Information criteria built on the unscaled log-likelihood, \\-2\ell +
  2p\\ and \\-2\ell + p\log n\\ with \\p\\ the number of estimated
  parameters.

- iterations:

  How many iterations the run that was kept took.

- converged:

  Logical of length 1: whether a stopping rule confirmed convergence. A
  run that exhausted its iteration budget is not converged whatever
  point it reached.

- method:

  The optimization method actually used, which is not always the one
  asked for: `"Fisher scoring"` and `"Newton-Raphson"` fall back to
  `"BFGS"` when they fail, and this records which one produced the
  estimates.

- criterion:

  Which stopping rule ended the run, in the words optimizers7 reports
  it, such as `"gradient (max-norm) < 1e-06"`. Empty when none fired.

- note:

  Any remark the optimizer attached to the run, such as
  `"the line search found no acceptable step"`. Empty when there is
  none.

- counts:

  A named list of evaluation counts with components `f`, `g` and `h`:
  how many times the objective, its gradient and its Hessian were
  evaluated in the run that was kept.

- score:

  The max-norm of the score **per observation** at the reported optimum.
  This is the quantity the stopping rule tested, so it says how close to
  stationary a run ended and is the one number a non-converged fit is
  worth reading for.

- elapsed:

  Seconds spent optimizing, summed over every starting value and every
  fallback attempted, not just the run that was kept.

- level:

  The confidence level `ci` and `ci_eta` were built at.

## Value

An S7 object of class `distrib_fit`, with the properties above.

## Methods

Registered on this class:
[`coef()`](https://statmodels7.github.io/distributions7/reference/coef.distrib_fit.md),
[`confint()`](https://statmodels7.github.io/distributions7/reference/confint.distrib_fit.md),
[`logLik()`](https://statmodels7.github.io/distributions7/reference/logLik.distrib_fit.md),
[`plot()`](https://statmodels7.github.io/distributions7/reference/plot.distrib_fit.md),
[`print()`](https://statmodels7.github.io/distributions7/reference/print.distrib_fit.md),
[`simulate()`](https://statmodels7.github.io/distributions7/reference/simulate.distrib_fit.md),
[`vcov()`](https://statmodels7.github.io/distributions7/reference/vcov.distrib_fit.md)

[`coef()`](https://rdrr.io/r/stats/coef.html),
[`vcov()`](https://rdrr.io/r/stats/vcov.html) and
[`confint()`](https://rdrr.io/r/stats/confint.html) each take a `scale`
argument and report either scale from the stored components.

## See also

[`fit_distrib()`](https://statmodels7.github.io/distributions7/reference/fit_distrib.md),
which builds one;
[`print.distrib_fit()`](https://statmodels7.github.io/distributions7/reference/print.distrib_fit.md)
for the printed layout;
[`confint.distrib_fit()`](https://statmodels7.github.io/distributions7/reference/confint.distrib_fit.md)
to recompute an interval at another level.

## Examples

``` r
set.seed(1)
d <- gaussian1_distrib()
y <- distrib_rng(d, 200, list(mu = 1, sigma = 2))
fit <- fit_distrib(d, y)
S7::S7_inherits(fit, distrib_fit)
#> [1] TRUE

# The two scales are stored together, and eta is the link of the estimate.
rbind(parameter = coef(fit), link = coef(fit, scale = "link"))
#>                 mu     sigma
#> parameter 1.071079 1.8535432
#> link      1.071079 0.6170991
all.equal(fit@eta[["sigma"]], log(fit@coefficients[["sigma"]]))
#> [1] TRUE

# The parameter-scale interval is the link-scale one mapped through g^-1,
# so the lower limit of a scale stays positive whatever the sample.
fit@ci
#>          lower    upper
#> mu    0.814196 1.327963
#> sigma 1.680516 2.044386
all.equal(fit@ci[["sigma", "lower"]], exp(fit@ci_eta[["sigma", "lower"]]))
#> [1] TRUE

# What the optimizer did. 'score' is the max-norm of the score per
# observation at the point reported, which is what the rule tested.
c(iterations = fit@iterations, converged = fit@converged, score = fit@score)
#>   iterations    converged        score 
#> 2.000000e+00 1.000000e+00 7.864824e-11 
fit@criterion
#> [1] "gradient (max-norm) < 1e-06"
```
