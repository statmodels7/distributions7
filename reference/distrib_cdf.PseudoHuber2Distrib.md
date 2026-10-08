# Pseudo-Huber Distribution Function, Standard-Deviation Parametrization

Computes \\P(Y \le q)\\ by numerical integration of the density, taken
by the compiled rule of
[`compiled_cdf()`](https://statmodels7.github.io/distributions7/reference/compiled_cdf.md)
over the tail on the side of \\q\\ away from \\\mu\\.

## Arguments

- distrib:

  A `PseudoHuber2Distrib` object.

- q:

  A numeric vector of quantiles.

- theta:

  A named list with components `mu`, `sigma` and `nu`.

- lower.tail:

  Logical; if `FALSE`, \\P(Y \> q)\\ is returned.

- log.p:

  Logical; if `TRUE`, the logarithm is returned.

- ...:

  Unused.

## Value

A numeric vector of probabilities.

## See also

[`distrib_quantile.PseudoHuber2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.PseudoHuber2Distrib.md),
which inverts it.

## Examples

``` r
d <- pseudohuber2_distrib()
distrib_cdf(d, c(-1, 0.4, 2), list(mu = 0.4, sigma = 1.5, nu = 2))
#> [1] 0.1498899 0.5000000 0.8777881
```
