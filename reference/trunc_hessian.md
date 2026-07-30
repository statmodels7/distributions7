# Observed Hessian of a Truncated Distribution

\\d\_{ij}\ell_T = H\_{ij}(y) - M\_{ij} + m_i m_j\\.

## Usage

``` r
trunc_hessian(distrib, y, theta, scale = c("parameter", "link"), ...)
```

## Arguments

- distrib:

  A truncated distribution object.

- y:

  A numeric vector of observations.

- theta:

  A named list of parameters.

## Value

A named list, one component per Hessian entry.

## See also

[`trunc_M`](https://statmodels7.github.io/distributions7/reference/trunc_M.md),
[`trunc_score_mean`](https://statmodels7.github.io/distributions7/reference/trunc_score_mean.md)
