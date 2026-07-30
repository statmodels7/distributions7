# Expected Derivatives by Numerical Integration

Integrates each observed derivative component directly against the
density, component by component, through
[`expectation`](https://statmodels7.github.io/distributions7/reference/expectation.md).

## Usage

``` r
expected_by_integrate(distrib, y, theta, order)
```

## Arguments

- distrib:

  An object inheriting from class `"distrib"`.

- y:

  A numeric vector of observations; only its length is used, to recycle
  the result.

- theta:

  A named list of parameters.

- order:

  The derivative order, 2 to 4.

## Value

A named list of expected derivative component vectors, each of length
`length(y)`.

## Details

This estimates \\\mathbb{E}\[\partial^k \ell\]\\ literally. For a
regular model that is the quantity wanted; for a non-regular one it is
not the information, and
[`expected_by_bartlett`](https://statmodels7.github.io/distributions7/reference/expected_by_bartlett.md)
is the route that stays valid (see
[`laplace_distrib`](https://statmodels7.github.io/distributions7/reference/laplace_distrib.md)).

Quadrature is unreliable when the observed derivative is itself a finite
difference, since it then integrates numerical noise, so the error is
caught and re-raised naming the component and pointing at the
alternatives rather than surfacing as an opaque failure from the
integrator.

## See also

[`expected_derivative_methods`](https://statmodels7.github.io/distributions7/reference/expected_derivative_methods.md)
