# The Expectation of a Smooth Function Under a Truncated Family

Computes \\E_T\[f(Y)\]\\ for a smooth `f` by
[`trunc_rule_expectation()`](https://statmodels7.github.io/distributions7/reference/trunc_rule_expectation.md),
with the masses of the parent's atoms inside the retained interval
added, and by
[`expectation()`](https://statmodels7.github.io/distributions7/reference/expectation.md)
where the rule is not taken. The derivatives of a truncated family's
expected information call it.

## Usage

``` r
trunc_expectation(distrib, f, theta, ...)
```

## Arguments

- distrib:

  A truncated continuous family.

- f:

  The function whose expectation is taken.

- theta:

  A named list of parameters.

- ...:

  Further arguments passed to `f`.

## Value

A numeric vector with one value per parameter combination.
