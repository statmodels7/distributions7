# Second-Order Truncated Moment of the Parent's Derivatives

\\M\_{ij} = \mathbb{E}\_T\[H\_{ij} + s_i s_j\]\\, the quantity entering
the truncated Hessian as \\d\_{ij}\ell_T = H\_{ij}(y) - M\_{ij} + m_i
m_j\\.

## Usage

``` r
trunc_M(distrib, theta)
```

## Arguments

- distrib:

  A truncated distribution object.

- theta:

  A named list of parameters.

## Value

A named list, one component per Hessian entry.

## See also

[`trunc_score_mean`](https://statmodels7.github.io/distributions7/reference/trunc_score_mean.md)
