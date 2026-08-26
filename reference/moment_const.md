# Recycle a Constant Moment to the Length of the Parameters

Returns `value` repeated to the length the parameters imply. A moment
that does not depend on the parameters still has to come back one value
per setting, so that every moment method returns the same shape whatever
the family; this helper is what enforces that for the constants.

## Usage

``` r
moment_const(theta, k, value)
```

## Arguments

- theta:

  An aligned named list of parameters, as
  [`align_theta()`](https://statmodels7.github.io/distributions7/reference/align_theta.md)
  returns.

- k:

  How many of its components to read the length from. A single whole
  number, usually the number of parameters the family carries.

- value:

  The constant to recycle. A numeric vector of length 1, which is `NaN`
  for a moment that does not exist and a number for one that is fixed by
  the family.

## Value

A numeric vector of length `max(lengths(theta[seq_len(k)]))`, every
element `value`.

## See also

[`skewness.CauchyDistrib()`](https://statmodels7.github.io/distributions7/reference/skewness.CauchyDistrib.md),
which uses it for a moment that does not exist, and
[`kurtosis.GumbelDistrib()`](https://statmodels7.github.io/distributions7/reference/kurtosis.GumbelDistrib.md),
which uses it for one that is constant.

## Examples

``` r
# One value per parameter setting, even though the value is fixed.
skewness(gumbel_distrib(), list(mu = c(0, 1, 2), sigma = 1))
#> [1] 1.139547 1.139547 1.139547
```
