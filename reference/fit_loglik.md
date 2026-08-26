# Total Log-Likelihood

Sums the log-density over the observations, \\\sum_i \log f(y_i;
\theta)\\, and returns the one number the fit maximizes. Every
observation is read at the same \\\theta\\, which is the i.i.d.
assumption the fitting layer makes.

The value is **not** divided by \\n\\.
[`fit_distrib()`](https://statmodels7.github.io/distributions7/reference/fit_distrib.md)
scales it when it hands the objective to an optimizer, and recomputes it
unscaled here for the `loglik`, `aic` and `bic` a fit reports.

## Usage

``` r
fit_loglik(distrib, y, theta)
```

## Arguments

- distrib:

  An object inheriting from `distrib`.

- y:

  A numeric vector of observations, or the response matrix of a
  multivariate family.

- theta:

  A named list of parameters on the parameter scale, aligned to
  `distrib@params`.

## Value

A single number. It is `-Inf` when any observation has zero density, and
`NaN` when the density itself is not computable at `theta`.

## See also

[`distrib_pdf()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.md),
which supplies the terms;
[`logLik.distrib_fit()`](https://statmodels7.github.io/distributions7/reference/logLik.distrib_fit.md)
for the value a fit carries.
