# Expected Derivatives of the Beta-Binomial by Exact Summation

Averages every component of a derivative of the given order over the
support \\\\0, \dots, n\\\\, weighted by the mass function. The support
being finite, the average is the expectation exactly, with no quadrature
error and no sampling error.

## Usage

``` r
betabinom2_expected(
  distrib,
  y,
  theta,
  order,
  derivs = function(yy, aa, bb) {
     betabinom2_derivs(yy, aa, bb, distrib@size, order,
    distrib@params)
 }
)
```

## Arguments

- distrib:

  A
  [`BetaBinom2Distrib()`](https://statmodels7.github.io/distributions7/reference/BetaBinom2Distrib.md)
  or
  [`BetaBinom1Distrib()`](https://statmodels7.github.io/distributions7/reference/BetaBinom1Distrib.md)
  object, read for its `size` and `params`.

- y:

  A numeric vector, used only for its length.

- theta:

  A named list of the two parameters, each a numeric vector of length 1
  or of the length of `y`. Parameters of length 1 give one constant per
  component, repeated to the length of `y`; longer ones give one sum per
  observation.

- order:

  The derivative order, an integer from 2 to 4.

- derivs:

  A function of the counts and the two parameters, in that order, that
  returns the components being averaged. The default is
  [`betabinom2_derivs()`](https://statmodels7.github.io/distributions7/reference/betabinom2_derivs.md)
  at `order`;
  [`betabinom1_expected()`](https://statmodels7.github.io/distributions7/reference/betabinom1_expected.md)
  passes
  [`betabinom1_components()`](https://statmodels7.github.io/distributions7/reference/betabinom1_components.md).

## Value

A named list of component vectors, each of the recycled length of `y`
and the parameters. At order 2 the names and their order are
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
