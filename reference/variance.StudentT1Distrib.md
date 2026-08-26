# Variance of the Student t Distribution

Closed form: \\\operatorname{Var}(Y) = \sigma^2\nu/(\nu-2)\\ for \\\nu
\> 2\\, `Inf` for \\1 \< \nu \le 2\\, and `NaN` at or below one. The
three answers are three different statements: a finite value, a second
moment that diverges to infinity while the first exists, and a second
moment that is undefined because the first is.

## Arguments

- x:

  A `StudentT1Distrib`, from
  [`student_t1_distrib()`](https://statmodels7.github.io/distributions7/reference/student_t1_distrib.md).

- theta:

  A named list with components `mu`, `sigma` (positive) and `nu`
  (positive), each a numeric vector of length 1 or `n`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of variances, of length equal to the longest of the
three components, `Inf` where \\1 \< \nu \le 2\\ and `NaN` where \\\nu
\le 1\\.

## Details

The factor \\\nu/(\nu-2)\\ exceeds one at every finite \\\nu\\ and falls
onto it as the degrees of freedom grow, so the scale is a lower bound on
the standard deviation and is attained only in the Gaussian limit. It is
3 at \\\nu = 3\\ and 1.25 at \\\nu = 10\\, so the inflation is large
only close to the threshold.

## Notation

\\\sigma \> 0\\ is the scale and \\\nu \> 0\\ the degrees of freedom.

## See also

[`kurtosis.StudentT1Distrib()`](https://statmodels7.github.io/distributions7/reference/kurtosis.StudentT1Distrib.md),
whose threshold is 4;
[`variance.CauchyDistrib()`](https://statmodels7.github.io/distributions7/reference/variance.CauchyDistrib.md),
the \\\nu = 1\\ case;
[`student_t1_distrib()`](https://statmodels7.github.io/distributions7/reference/student_t1_distrib.md).

## Examples

``` r
d <- student_t1_distrib()

# NaN, infinite and finite, in that order, as nu crosses its two thresholds.
variance(d, list(mu = 0, sigma = 1, nu = c(1, 1.5, 2, 3)))
#> [1] NaN Inf Inf   3

# The scale is a lower bound, approached as the degrees of freedom grow.
variance(d, list(mu = 0, sigma = 2, nu = c(3, 10, 1e4)))
#> [1] 12.0000  5.0000  4.0008
```
