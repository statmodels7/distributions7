# Compute an Expected Derivative by a Method

Validates `method` and routes to the computation that the object
selects. It is the shared entry point for the expected derivatives that
a family does not write out, and for those that a caller requests with
an explicit method. A `NULL` method is resolved to the default of the
order:
[`expect_opg()`](https://statmodels7.github.io/distributions7/reference/expect_opg.md)
at order 2 and
[`expect_integrate()`](https://statmodels7.github.io/distributions7/reference/expect_integrate.md)
above it.

## Usage

``` r
expected_derivative(distrib, y, theta, order, method = NULL)
```

## Arguments

- distrib:

  An object inheriting from class `"distrib"`.

- y:

  A numeric vector of observations.

- theta:

  A named list of parameters.

- order:

  The derivative order, 2 to 4.

- method:

  An object inheriting from
  [ExpectMethod](https://statmodels7.github.io/distributions7/reference/ExpectMethod.md),
  or `NULL` (the default). Any other value, a character string included,
  signals an error.

## Value

A named list of expected derivative component vectors.

## Details

The computation of each class is the method of the internal generic
[`expected_by_method()`](https://statmodels7.github.io/distributions7/reference/expected_by_method.md);
see
[expected_derivative_methods](https://statmodels7.github.io/distributions7/reference/expected_derivative_methods.md)
for what the three classes compute, their costs and their defaults.

## See also

[expected_derivative_methods](https://statmodels7.github.io/distributions7/reference/expected_derivative_methods.md)
