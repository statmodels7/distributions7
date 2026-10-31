# Carry an Analytic Sixth Order to the Link Scale

Composes the sixth-order parameter-scale derivatives of a family that
registers its own
[`distrib_deriv6()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv6.md)
method with the parameter-scale orders one to five, through
[`to_link_scale()`](https://statmodels7.github.io/distributions7/reference/to_link_scale.md),
which uses the derivatives of the inverse links up to the sixth.

## Usage

``` r
deriv6_scale(distrib, y, theta, res)
```

## Arguments

- distrib:

  A distribution object.

- y, theta:

  The observations and the parameters, as given to the generic.

- res:

  The sixth order on the parameter scale.

## Value

A named list of link-scale sixth-order components.

## See also

[`has_exact_deriv6()`](https://statmodels7.github.io/distributions7/reference/has_exact_deriv6.md),
[`deriv5_scale()`](https://statmodels7.github.io/distributions7/reference/deriv5_scale.md).
