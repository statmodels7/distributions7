# Parameter Derivatives of a Reparametrized Distribution

The body the four derivative methods share: aligns the new parameters,
maps them to the parent's, fetches the map's tables and hands all three
to
[`chain_derivatives()`](https://statmodels7.github.io/distributions7/reference/chain_derivatives.md).
The `order` and `expected` arguments are what separate the four
registrations.

## Usage

``` r
reparam_chain(
  distrib,
  y,
  theta,
  order,
  expected = FALSE,
  approx = "opg",
  nsim = 10000
)
```

## Arguments

- distrib:

  A reparametrized distribution.

- y:

  The response, a numeric vector.

- theta:

  A named list of the new parameters, on the new parameter scale.

- order:

  The derivative order, 1 to 4.

- expected:

  Should the expected derivatives be carried? A single logical, `FALSE`
  by default. `TRUE` is meaningful from order 2 up, the score having
  mean zero.

## Value

A named list of numeric vectors, keyed as
[`deriv_names(distrib@params, order)`](https://statmodels7.github.io/distributions7/reference/deriv_names.md).

## See also

[`chain_derivatives()`](https://statmodels7.github.io/distributions7/reference/chain_derivatives.md),
which does the work;
[`distrib_gradient.ReparamContinuousDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.ReparamContinuousDistrib.md)
and its three siblings, the registrations;
[`reparametrize()`](https://statmodels7.github.io/distributions7/reference/reparametrize.md).
