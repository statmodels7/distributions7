# Skewness of a Distribution or Sample

Computes the skewness (third standardized moment):

\$\$\gamma_1 = \mathbb{E}\\\left\[\left(\frac{Y -
\mathbb{E}\[Y\]}{\operatorname{sd}(Y)}\right)^{3}\right\].\$\$

For `distrib` objects it is evaluated numerically via
[`moment`](https://statmodels7.github.io/distributions7/reference/moment.md);
for numeric vectors the sample skewness (population denominator) is
returned.

## Usage

``` r
skewness(x, ...)
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
[`kurtosis`](https://statmodels7.github.io/distributions7/reference/kurtosis.md)

## Examples

``` r
skewness(gaussian1_distrib(), list(mu = 0, sigma = 1))
#> [1] 0
skewness(gamma2_distrib(), list(mu = 2, sigma2 = 1))
#> [1] 1
```
