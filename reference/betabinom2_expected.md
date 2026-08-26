# Expected Derivatives of the Beta-Binomial by Exact Summation

Averages every component of a derivative of the given order over the
support \\\\0, \dots, n\\\\, weighted by the mass function. The support
being finite, the average is the expectation exactly, with no quadrature
error and no sampling error.

## Usage

``` r
betabinom2_expected(distrib, y, theta, order)
```

## Arguments

- distrib:

  A
  [`BetaBinom2Distrib()`](https://statmodels7.github.io/distributions7/reference/BetaBinom2Distrib.md)
  object, read for its `size` and `params`.

- y:

  A numeric vector, used only for its length: the result is one constant
  per component, repeated to that length.

- theta:

  A named list with components `alpha` and `beta`, each a numeric vector
  of length 1, strictly positive.

- order:

  The derivative order, an integer from 2 to 4.

## Value

A named list of component vectors, each of length `length(y)` and
constant along it. At order 2 the names and their order are
[`hess_names()`](https://statmodels7.github.io/distributions7/reference/hess_names.md)'s;
above it they are
[`deriv_names()`](https://statmodels7.github.io/distributions7/reference/deriv_names.md)'s.

## See also

[`distrib_expected_hessian.BetaBinom2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.BetaBinom2Distrib.md),
[`distrib_deriv3.BetaBinom2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.BetaBinom2Distrib.md)
and
[`distrib_deriv4.BetaBinom2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.BetaBinom2Distrib.md),
which call this;
[`betabinom2_derivs()`](https://statmodels7.github.io/distributions7/reference/betabinom2_derivs.md)
for the components being averaged.
