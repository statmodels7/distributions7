# Mean of the Parent's Score Under the Truncated Law

\\m_i = \mathbb{E}\_T\[s_i\]\\, the quantity that recentres the parent's
score: the truncated score is \\d_i \ell_T = s_i(y) - m_i\\.

## Usage

``` r
trunc_score_mean(distrib, theta)
```

## Arguments

- distrib:

  A truncated distribution object.

- theta:

  A named list of parameters.

## Value

A named list, one component per parameter.

## Details

Taken from the cdf derivatives where those are exact, and from
quadrature otherwise.

## See also

[`trunc_score_mean_quad`](https://statmodels7.github.io/distributions7/reference/trunc_score_mean_quad.md),
[`trunc_M`](https://statmodels7.github.io/distributions7/reference/trunc_M.md)
