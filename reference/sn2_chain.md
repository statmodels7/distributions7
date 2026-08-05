# Derivatives of the Skew Normal in Its Centred Parametrisation

The parent's derivatives carried into the centred coordinates by the
partition sum of
[`chain_derivatives`](https://statmodels7.github.io/distributions7/reference/chain_derivatives.md).

## Usage

``` r
sn2_chain(distrib, y, theta, order, expected = FALSE)
```

## Arguments

- distrib:

  A
  [`SkewNormal2Distrib`](https://statmodels7.github.io/distributions7/reference/SkewNormal2Distrib.md)
  object.

- y:

  The response.

- theta:

  A list with `mu`, `sigma` and `gamma1`.

- order:

  The derivative order.

- expected:

  Logical; if `TRUE`, carries the expected derivatives.

## Value

A named list of component vectors.

## See also

[`skewnormal2_distrib`](https://statmodels7.github.io/distributions7/reference/skewnormal2_distrib.md)
