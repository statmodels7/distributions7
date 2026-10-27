# Route an Explicit Expectation Method

Computes the expected derivative of an order by the method that the
object selects. A univariate family goes to
[`expected_derivative()`](https://statmodels7.github.io/distributions7/reference/expected_derivative.md).
A multivariate family has no quadrature over the response, so it goes to
[`mv_expected_sampled()`](https://statmodels7.github.io/distributions7/reference/mv_expected_sampled.md),
which takes the expectation by sampling.

## Usage

``` r
route_expected(distrib, y, theta, order, method)
```

## Arguments

- distrib:

  An object inheriting from `distrib`.

- y:

  A numeric vector or matrix of observations.

- theta:

  A named list of parameters.

- order:

  The derivative order, 2 to 4.

- method:

  An object inheriting from
  [ExpectMethod](https://statmodels7.github.io/distributions7/reference/ExpectMethod.md).

## Value

A named list of expected derivative component vectors.

## See also

[`explicit_method()`](https://statmodels7.github.io/distributions7/reference/explicit_method.md),
[expected_derivative_methods](https://statmodels7.github.io/distributions7/reference/expected_derivative_methods.md)
