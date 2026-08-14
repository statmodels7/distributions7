# Derivatives of the Skew Normal in Its Centered Parametrization

The parent's derivatives carried into the centered coordinates by the
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

## Details

The map to the direct parametrization runs through \\r =
\sqrt\[3\]{2\gamma_1/(4-\pi)}\\, whose derivative grows like
\\\gamma_1^{-2/3}\\, so at zero skewness the map is not differentiable
and the chain rule is asked for a quantity that does not exist. The
first derivatives of the log-density survive the limit – the map's
factor cancels and they approach a finite value from both sides – but
the second ones diverge at that rate, which is a property of the
CENTERED parametrization and not of the family. It is rejected here,
where the map is used and the reason can be named, rather than left to
reach a comparison against `NA` several frames further on.

## See also

[`skewnormal2_distrib`](https://statmodels7.github.io/distributions7/reference/skewnormal2_distrib.md)
