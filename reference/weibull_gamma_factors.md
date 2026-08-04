# Gamma Factors of a Weibull's Moments

Returns \\g_k = \Gamma(1 + k/\sigma)\\ for \\k = 1, \ldots, 4\\, from
which every moment of a Weibull follows as \\E\[Y^k\] = \mu^k g_k\\.

## Usage

``` r
weibull_gamma_factors(sigma, k = 4L)
```

## Arguments

- sigma:

  The shape parameter.

- k:

  How many factors to return.

## Value

A list of numeric vectors, `g1` to `gk`.
