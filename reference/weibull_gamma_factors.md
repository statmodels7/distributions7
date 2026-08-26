# Gamma Factors of a Weibull's Moments

Returns \\g_k = \Gamma(1 + k/\sigma)\\ for \\k = 1, \ldots, K\\. Every
raw moment of a Weibull is one of these times a power of the scale,
\\E\[Y^k\] = \mu^k g_k\\, so the four moment methods of the family share
this one helper and differ only in which combination of the factors they
assemble.

## Usage

``` r
weibull_gamma_factors(sigma, k = 4L)
```

## Arguments

- sigma:

  The shape parameter, a positive numeric vector. Each factor diverges
  as the shape approaches zero, the corresponding moment not existing
  there.

- k:

  How many factors to return. A single whole number, 4 by default, the
  number the excess kurtosis needs; the variance needs 2 and the
  skewness 3, and each caller asks for only what it uses.

## Value

A named list of `k` numeric vectors, `g1` to `gk`, each the length of
`sigma`.

## See also

[`variance.Weibull1Distrib()`](https://statmodels7.github.io/distributions7/reference/variance.Weibull1Distrib.md),
[`skewness.Weibull1Distrib()`](https://statmodels7.github.io/distributions7/reference/skewness.Weibull1Distrib.md)
and
[`kurtosis.Weibull1Distrib()`](https://statmodels7.github.io/distributions7/reference/kurtosis.Weibull1Distrib.md)
for the three consumers.

## Examples

``` r
# At shape 1 the family is exponential and g_k is k factorial.
distributions7:::weibull_gamma_factors(1, 4)
#> $g1
#> [1] 1
#> 
#> $g2
#> [1] 2
#> 
#> $g3
#> [1] 6
#> 
#> $g4
#> [1] 24
#> 
```
