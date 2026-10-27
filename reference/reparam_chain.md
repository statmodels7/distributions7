# Parameter Derivatives of a Reparametrized Distribution

The body that the derivative methods share: it aligns the new
parameters, maps them to the parent's, fetches the map's tables and
hands all three to
[`chain_derivatives()`](https://statmodels7.github.io/distributions7/reference/chain_derivatives.md).
The registrations differ in the `order` and `expected` arguments that
they pass.

## Usage

``` r
reparam_chain(distrib, y, theta, order, expected = FALSE, method = NULL)
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

  Should the expected derivatives be carried? It is a single logical,
  with default `FALSE`. `TRUE` is meaningful from order 2 up, because
  the score has mean zero.

- method:

  `NULL` (the default) or an object inheriting from
  [ExpectMethod](https://statmodels7.github.io/distributions7/reference/ExpectMethod.md).
  The generic applies an object before dispatching on the class, so this
  method receives `NULL` and forwards `NULL` to the parent, which
  returns its closed form where it has one. See
  [expected_derivative_methods](https://statmodels7.github.io/distributions7/reference/expected_derivative_methods.md).

## Value

A named list of numeric vectors, keyed as
[`deriv_names(distrib@params, order)`](https://statmodels7.github.io/distributions7/reference/deriv_names.md).

## See also

[`chain_derivatives()`](https://statmodels7.github.io/distributions7/reference/chain_derivatives.md),
which does the work;
[`distrib_gradient.ReparamContinuousDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.ReparamContinuousDistrib.md)
and its three siblings, the registrations;
[`reparametrize()`](https://statmodels7.github.io/distributions7/reference/reparametrize.md).
