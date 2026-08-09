# Excess Kurtosis of a Distribution or Sample

Computes the excess kurtosis, the fourth standardized moment less the
three a Gaussian has:

\$\$\gamma_2 = \mathbb{E}\\\left\[\left(\frac{Y -
\mathbb{E}\[Y\]}{\operatorname{sd}(Y)}\right)^{4}\right\] - 3.\$\$

For `distrib` objects it is evaluated numerically via
[`moment`](https://statmodels7.github.io/distributions7/reference/moment.md);
for numeric vectors the sample excess kurtosis (population denominator)
is returned.

## Usage

``` r
kurtosis(x, ...)
```

## Arguments

- x:

  An object inheriting from class `"distrib"`, or a numeric vector.

- ...:

  For `distrib` objects: `theta` and further arguments passed to
  [`moment`](https://statmodels7.github.io/distributions7/reference/moment.md).
  For numeric vectors: `na.rm`.

## Value

A numeric vector.

## See also

[`expectation`](https://statmodels7.github.io/distributions7/reference/expectation.md),
[`moment`](https://statmodels7.github.io/distributions7/reference/moment.md),
[`variance`](https://statmodels7.github.io/distributions7/reference/variance.md),
[`std_dev`](https://statmodels7.github.io/distributions7/reference/std_dev.md),
[`skewness`](https://statmodels7.github.io/distributions7/reference/skewness.md)

## Examples

``` r
kurtosis(gaussian1_distrib(), list(mu = 0, sigma = 1))
#> [1] 0
kurtosis(gamma2_distrib(), list(mu = 2, sigma2 = 1))
#> [1] 1.5
```
