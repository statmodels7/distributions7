# Expected Derivatives of the Beta-Binomial by Exact Summation

Averages a derivative over the support \\\\0, \dots, n\\\\ weighted by
the mass function, which is exact because the support is finite.

## Usage

``` r
betabinom2_expected(distrib, y, theta, order)
```

## Arguments

- distrib:

  A
  [`BetaBinom2Distrib`](https://statmodels7.github.io/distributions7/reference/BetaBinom2Distrib.md)
  object.

- y:

  A numeric vector, used only for its length.

- theta:

  A list with `alpha` and `beta`.

- order:

  The derivative order.

## Value

A named list of component vectors.

## See also

[`betabinom2_distrib`](https://statmodels7.github.io/distributions7/reference/betabinom2_distrib.md)
