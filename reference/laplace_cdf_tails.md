# The Laplace Tails on the Log Scale

These are internal helpers of the two Laplace families.
`laplace_cdf_tails()` returns the probability that `lower.tail` and
`log.p` select at the standardized residual `r`, and computes each tail
from its own arm: \\\log(1/2) - \|r\|\\ for the far tail and
`log1p(-exp(-|r|)/2)` for the near one. `laplace_std_quantile()` returns
the standard Laplace quantile from the logarithm of the smaller tail,
\\\pm(\log P + \log 2)\\.

## Usage

``` r
laplace_cdf_tails(r, lower.tail, log.p)

laplace_std_quantile(p, lower.tail, log.p)
```

## Arguments

- r:

  A numeric vector of standardized residuals.

- lower.tail, log.p:

  Logical of length 1.

- p:

  A numeric vector of probabilities, or their logarithms.

## Value

A numeric vector.

## Examples

``` r
distributions7:::laplace_cdf_tails(c(-50, 0, 50), TRUE, TRUE)
#> [1] -5.069315e+01 -6.931472e-01 -9.643749e-23
distributions7:::laplace_std_quantile(1e-20, TRUE, FALSE)
#> [1] -45.35855
```
