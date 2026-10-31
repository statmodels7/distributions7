# Carry an Analytic Fifth Derivative of the Information to the Link Scale

Composes the fifth derivative of the expected information of a family
that registers its own
[`distrib_d5expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_d5expected_hessian.md)
method with the parameter-scale information and its derivatives of
orders one to four, through
[`dexpected_link_higher()`](https://statmodels7.github.io/distributions7/reference/dexpected_link_higher.md),
which uses the derivatives of the inverse links up to the sixth.

## Usage

``` r
d5expected_scale(distrib, y, theta, res)
```

## Arguments

- distrib:

  A distribution object.

- y, theta:

  The observations and the parameters, as given to the generic.

- res:

  The fifth derivative on the parameter scale.

## Value

A named list keyed as
[`d5expected_names()`](https://statmodels7.github.io/distributions7/reference/d5expected_names.md).

## See also

[`has_exact_d5expected()`](https://statmodels7.github.io/distributions7/reference/has_exact_d5expected.md),
[`d34expected_analytic()`](https://statmodels7.github.io/distributions7/reference/d34expected_analytic.md).
