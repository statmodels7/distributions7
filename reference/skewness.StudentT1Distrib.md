# Skewness of the Student t Distribution

Zero for \\\nu \> 3\\, and `NaN` at or below three, where the third
moment does not exist. The density is symmetric about \\\mu\\, so the
value is zero wherever it is defined. The threshold records where the
third moment exists.

## Arguments

- x:

  A `StudentT1Distrib`, from
  [`student_t1_distrib()`](https://statmodels7.github.io/distributions7/reference/student_t1_distrib.md).

- theta:

  A named list with components `mu`, `sigma` (positive) and `nu`
  (positive), each a numeric vector of length 1 or `n`. Only `nu` enters
  the value; settings with \\\nu \le 3\\ give `NaN`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector, of length equal to the longest of the three
components, zero where \\\nu \> 3\\ and `NaN` elsewhere.

## Notation

\\\mu\\ is the location and \\\nu \> 0\\ the degrees of freedom.

## See also

[`kurtosis.StudentT1Distrib()`](https://statmodels7.github.io/distributions7/reference/kurtosis.StudentT1Distrib.md),
whose threshold is 4;
[`skewness.SkewTDistrib()`](https://statmodels7.github.io/distributions7/reference/skewness.SkewTDistrib.md),
the asymmetric extension;
[`student_t1_distrib()`](https://statmodels7.github.io/distributions7/reference/student_t1_distrib.md).

## Examples

``` r
d <- student_t1_distrib()

# Undefined at or below three degrees of freedom, then zero by symmetry.
skewness(d, list(mu = 0, sigma = 1, nu = c(2, 3, 4)))
#> [1] NaN NaN   0
```
