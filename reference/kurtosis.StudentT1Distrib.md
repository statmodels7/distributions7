# Excess Kurtosis of the Student t Distribution

Closed form: \\\gamma_2 = 6/(\nu-4)\\ for \\\nu \> 4\\, `Inf` for \\2 \<
\nu \le 4\\, and `NaN` at or below two. The value is the toolkit's
standard example of tail weight running with a shape parameter: it is 6
at \\\nu = 5\\, 1 at \\\nu = 10\\ and tends to the Gaussian's zero as
the degrees of freedom grow.

## Arguments

- x:

  A `StudentT1Distrib`, from
  [`student_t1_distrib()`](https://statmodels7.github.io/distributions7/reference/student_t1_distrib.md).

- theta:

  A named list with components `mu`, `sigma` (positive) and `nu`
  (positive), each a numeric vector of length 1 or `n`. Only `nu` enters
  the value.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of excess kurtoses, of length equal to the longest of
the three components, `Inf` where \\2 \< \nu \le 4\\ and `NaN` where
\\\nu \le 2\\.

## Notation

\\\nu \> 0\\ is the degrees of freedom. Neither the location nor the
scale enters a standardized moment.

## See also

[`skewness.StudentT1Distrib()`](https://statmodels7.github.io/distributions7/reference/skewness.StudentT1Distrib.md),
whose threshold is 3;
[`kurtosis.CauchyDistrib()`](https://statmodels7.github.io/distributions7/reference/kurtosis.CauchyDistrib.md),
the \\\nu = 1\\ case;
[`kurtosis.SkewTDistrib()`](https://statmodels7.github.io/distributions7/reference/kurtosis.SkewTDistrib.md);
[`student_t1_distrib()`](https://statmodels7.github.io/distributions7/reference/student_t1_distrib.md).

## Examples

``` r
d <- student_t1_distrib()

# NaN, infinite and finite, in that order, as nu crosses its two thresholds.
kurtosis(d, list(mu = 0, sigma = 1, nu = c(2, 3, 4, 5)))
#> [1] NaN Inf Inf   6

# 6 / (nu - 4), tending to the Gaussian's zero.
kurtosis(d, list(mu = 0, sigma = 1, nu = c(10, 100, 1e4)))
#> [1] 1.0000000000 0.0625000000 0.0006002401
```
