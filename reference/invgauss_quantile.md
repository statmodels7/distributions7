# Inverse Gaussian Quantile With a Far-Tail Refinement

Calls
[`statmod::qinvgauss()`](https://rdrr.io/pkg/statmod/man/invgauss.html)
and, where the requested tail probability is below \\e^{-5}\\, corrects
its value with
[`tail_quantile_refine()`](https://statmodels7.github.io/distributions7/reference/tail_quantile_refine.md):
there `qinvgauss()` can stop on its iteration limit with a value off by
tens of per cent, and its warning is not passed on for those elements.

## Usage

``` r
invgauss_quantile(distrib, p, theta, mean, dispersion, lower.tail, log.p)
```

## Arguments

- distrib:

  An `InvGauss1Distrib` or `InvGauss2Distrib` object.

- p, lower.tail, log.p:

  As in
  [`distrib_quantile()`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.md).

- theta:

  The family's named list of parameters.

- mean, dispersion:

  The parameters of `qinvgauss()`.

## Value

A numeric vector of quantiles.

## Examples

``` r
d <- invgauss1_distrib()
th <- list(mu = 1.7, phi = 2.3)
x <- distributions7:::invgauss_quantile(d, -800, th, 1.7, 2.3, FALSE, TRUE)
distrib_cdf(d, x, th, lower.tail = FALSE, log.p = TRUE)
#> [1] -800
```
