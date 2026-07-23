# S7 Class for Maximum-Likelihood Fits

Object returned by
[`fit_distrib`](https://statmodels7.github.io/distributions7/reference/fit_distrib.md),
holding the estimates on both the parameter and the link scale together
with their uncertainty.

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
  level = integer(0)
)
```

## Arguments

- distrib:

  The fitted `distrib` object.

- y:

  The observations the fit was computed from, kept so that the fitted
  distribution can be compared with the data (see
  [`plot.distrib_fit`](https://statmodels7.github.io/distributions7/reference/plot.distrib_fit.md)).

- n:

  Number of observations.

- coefficients:

  Named estimates on the parameter scale.

- se:

  Standard errors on the parameter scale (delta method).

- ci:

  Matrix of confidence limits on the parameter scale.

- eta:

  Estimates on the link scale.

- se_eta:

  Standard errors on the link scale.

- ci_eta:

  Matrix of confidence limits on the link scale.

- vcov:

  Variance-covariance matrix on the parameter scale.

- vcov_eta:

  Variance-covariance matrix on the link scale.

- loglik:

  Maximised log-likelihood.

- aic, bic:

  Information criteria.

- iterations:

  Number of iterations used.

- converged:

  Logical convergence flag.

- method:

  Optimisation method actually used.

- level:

  Confidence level.

## Methods

Methods implemented for this class:
[`coef()`](https://statmodels7.github.io/distributions7/reference/coef.distrib_fit.md),
[`logLik()`](https://statmodels7.github.io/distributions7/reference/logLik.distrib_fit.md),
[`plot()`](https://statmodels7.github.io/distributions7/reference/plot.distrib_fit.md),
[`print()`](https://statmodels7.github.io/distributions7/reference/print.distrib_fit.md),
[`simulate()`](https://statmodels7.github.io/distributions7/reference/simulate.distrib_fit.md),
[`vcov()`](https://statmodels7.github.io/distributions7/reference/vcov.distrib_fit.md)
