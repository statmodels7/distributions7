# Observed Derivatives of a Given Order

Routes to
[`distrib_gradient()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.md),
[`distrib_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.md),
[`distrib_deriv3()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.md),
[`distrib_deriv4()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.md)
or
[`distrib_deriv5()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv5.md)
according to `order`, so that code working at an order fixed only at run
time does not have to branch. The methods of
[expected_derivative_methods](https://statmodels7.github.io/distributions7/reference/expected_derivative_methods.md)
read it:
[`expect_integrate()`](https://statmodels7.github.io/distributions7/reference/expect_integrate.md)
and
[`expect_mc()`](https://statmodels7.github.io/distributions7/reference/expect_mc.md)
for the quantity that they average, and the score identity for the lower
orders that it multiplies together.

## Usage

``` r
observed_deriv(distrib, y, theta, order)
```

## Arguments

- distrib:

  An object inheriting from `distrib`.

- y:

  A numeric vector of observations.

- theta:

  A named list of parameters, aligned to `distrib@params`.

- order:

  The derivative order, a single integer from 1 to 5. Anything else
  signals an error.

## Value

A named list of derivative component vectors; each has length
`length(y)`. The keys are `distrib@params` at order 1,
[`hess_names()`](https://statmodels7.github.io/distributions7/reference/hess_names.md)
at order 2 (diagonal first) and
[`deriv_names()`](https://statmodels7.github.io/distributions7/reference/deriv_names.md)
above it (lexicographic), so a caller pairing two orders must match by
name.

## See also

[expected_derivative_methods](https://statmodels7.github.io/distributions7/reference/expected_derivative_methods.md),
which reads it;
[`deriv_names()`](https://statmodels7.github.io/distributions7/reference/deriv_names.md)
and
[`hess_names()`](https://statmodels7.github.io/distributions7/reference/hess_names.md)
for the two keyings;
[`expected_derivative()`](https://statmodels7.github.io/distributions7/reference/expected_derivative.md),
the router above it.
