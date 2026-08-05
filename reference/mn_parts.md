# The Pieces a Multinomial Derivative Needs

The probability vector and the simplex's first two derivatives, computed
once per call.

## Usage

``` r
mn_parts(distrib, theta)
```

## Arguments

- distrib:

  A
  [`MultinomialDistrib`](https://statmodels7.github.io/distributions7/reference/MultinomialDistrib.md)
  object.

- theta:

  A named list of parameters.

## Value

A list with `prob`, `A`, `B`, `idx` and `k`.

## See also

[`multinomial_distrib`](https://statmodels7.github.io/distributions7/reference/multinomial_distrib.md)
