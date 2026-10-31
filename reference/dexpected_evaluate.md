# Evaluate Polynomials of Observed Derivatives by a Method

Computes, for each polynomial, its value at the observations
([`expect_opg()`](https://statmodels7.github.io/distributions7/reference/expect_opg.md)),
its expectation by quadrature or summation
([`expect_integrate()`](https://statmodels7.github.io/distributions7/reference/expect_integrate.md))
or its Monte Carlo mean
([`expect_mc()`](https://statmodels7.github.io/distributions7/reference/expect_mc.md)).

## Usage

``` r
dexpected_evaluate(distrib, y, theta, polys, method, scale = "parameter")
```

## Arguments

- distrib:

  A univariate distribution.

- y:

  A numeric vector of observations; only its length is used by the
  methods that take an expectation.

- theta:

  A named list of parameters, aligned.

- polys:

  A list of polynomials.

- method:

  An object inheriting from
  [ExpectMethod](https://statmodels7.github.io/distributions7/reference/ExpectMethod.md).

- scale:

  `"parameter"` (the default) or `"link"`, the scale of the observed
  derivatives in the polynomials.

## Value

A list of numeric vectors, each of length `length(y)`.

## See also

[`dexpected_by_method()`](https://statmodels7.github.io/distributions7/reference/dexpected_by_method.md)
