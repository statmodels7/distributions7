# Third-Order Derivatives of a Reparametrized Distribution

The partition sum at order three, observed or expected.

## Usage

``` r
reparam_deriv3(
  distrib,
  y,
  theta,
  expected = FALSE,
  scale = c("parameter", "link"),
  approx = c("integrate", "bartlett", "mc", "opg"),
  nsim = 10000,
  ...
)
```

## Arguments

- distrib:

  A reparametrized distribution.

- y:

  The response.

- theta:

  A named list of the new parameters.

- expected:

  Logical; if `TRUE`, carries the expected derivatives.

- scale:

  Either `"parameter"` or `"link"`; handled by the generic.

- approx:

  Passed to the parent.

- nsim:

  Passed to the parent.

- ...:

  Unused.

## Value

A named list of third-derivative components.

## See also

[`reparametrize`](https://statmodels7.github.io/distributions7/reference/reparametrize.md)
