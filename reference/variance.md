# Variance of a Distribution or Sample

Computes the variance. For `distrib` objects the second central moment
is evaluated numerically (analytical methods may override this for
specific distributions); for numeric vectors the sample variance
[`var`](https://rdrr.io/r/stats/cor.html) is returned.

## Usage

``` r
variance(x, ...)
```

## Arguments

- x:

  An object inheriting from class `"distrib"`, or a numeric vector.

- ...:

  For `distrib` objects: `theta` (a named list of parameters) and
  further arguments passed to
  [`moment`](https://statmodels7.github.io/distributions7/reference/moment.md).
  For numeric vectors: `na.rm`.

## Value

A numeric vector.

## Details

\$\$\operatorname{Var}(Y) = \mathbb{E}\left\[(Y -
\mathbb{E}\[Y\])^{2}\right\].\$\$

A family with a closed form registers its own method; otherwise the
expectation is the quadrature or the exact sum of
[`expectation`](https://statmodels7.github.io/distributions7/reference/expectation.md).

## Examples

``` r
variance(gaussian1_distrib(), list(mu = 0, sigma = 2))
#> [1] 4
variance(poisson_distrib(), list(mu = 3))
#> [1] 3
```
