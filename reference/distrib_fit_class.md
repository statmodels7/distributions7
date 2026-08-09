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

  Maximized log-likelihood.

- aic, bic:

  Information criteria.

- iterations:

  Number of iterations used.

- converged:

  Logical convergence flag.

- method:

  Optimization method actually used.

- criterion:

  Which stopping rule ended the run, as optimizers7 reports it.

- note:

  Any remark the optimizer attached to the run.

- counts:

  How many times the objective and its gradient were evaluated.

- score:

  The max-norm of the score **per observation** at the reported optimum,
  which is the quantity the stopping rule tested.

- elapsed:

  Seconds spent optimizing, summed over every starting value and every
  fallback attempted.

- level:

  Confidence level.

## Value

An object of class `distrib_fit`.

## Methods

Methods implemented for this class:
[`coef()`](https://statmodels7.github.io/distributions7/reference/coef.distrib_fit.md),
[`logLik()`](https://statmodels7.github.io/distributions7/reference/logLik.distrib_fit.md),
[`plot()`](https://statmodels7.github.io/distributions7/reference/plot.distrib_fit.md),
[`print()`](https://statmodels7.github.io/distributions7/reference/print.distrib_fit.md),
[`simulate()`](https://statmodels7.github.io/distributions7/reference/simulate.distrib_fit.md),
[`vcov()`](https://statmodels7.github.io/distributions7/reference/vcov.distrib_fit.md)

## See also

[`fit_distrib`](https://statmodels7.github.io/distributions7/reference/fit_distrib.md)

## Examples

``` r
set.seed(1)
y <- distrib_rng(gaussian1_distrib(), 200, list(mu = 1, sigma = 2))
fit <- fit_distrib(gaussian1_distrib(), y)
S7::S7_inherits(fit, distrib_fit)
#> [1] TRUE
coef(fit)
#>       mu    sigma 
#> 1.071079 1.853543 
logLik(fit)
#> 'log Lik.' -407.2075 (df=2)
```
