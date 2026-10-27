# Expected Derivative by the Route of a Method

Internal generic that dispatches on the class of `method` and computes
the expected derivative of a given order. The method for
[ExpectOPG](https://statmodels7.github.io/distributions7/reference/ExpectMethod.md)
calls
[`expected_by_opg()`](https://statmodels7.github.io/distributions7/reference/expected_by_opg.md).
The method for
[ExpectIntegrate](https://statmodels7.github.io/distributions7/reference/ExpectMethod.md)
calls
[`expected_by_integrate()`](https://statmodels7.github.io/distributions7/reference/expected_by_integrate.md)
for `identity = "hessian"` and
[`expected_by_bartlett()`](https://statmodels7.github.io/distributions7/reference/expected_by_bartlett.md)
for `identity = "score"`. The method for
[ExpectMC](https://statmodels7.github.io/distributions7/reference/ExpectMethod.md)
calls
[`expected_by_mc()`](https://statmodels7.github.io/distributions7/reference/expected_by_mc.md)
with the identity of the object. A subclass of
[ExpectMethod](https://statmodels7.github.io/distributions7/reference/ExpectMethod.md)
defined elsewhere registers its own method.

## Usage

``` r
expected_by_method(method, distrib, y, theta, order)
```

## Arguments

- method:

  An object inheriting from
  [ExpectMethod](https://statmodels7.github.io/distributions7/reference/ExpectMethod.md).

- distrib:

  An object inheriting from `distrib`.

- y:

  A numeric vector of observations; only its length is used by the
  methods that take an expectation.

- theta:

  A named list of parameters.

- order:

  The derivative order, 2 to 4.

## Value

A named list of expected derivative component vectors.

## See also

[`expected_derivative()`](https://statmodels7.github.io/distributions7/reference/expected_derivative.md),
[expected_derivative_methods](https://statmodels7.github.io/distributions7/reference/expected_derivative_methods.md)
