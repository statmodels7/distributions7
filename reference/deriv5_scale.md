# Carry an Analytic Fifth Order to the Link Scale

The
[`distrib_deriv5()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv5.md)
generic, unlike those of orders one to four, does not apply the link
scale itself, because its numerical default differentiates on the link
scale directly. A family with its own fifth order returns the parameter
scale and passes it here, which composes it with the parameter-scale
orders one to four through
[`to_link_scale()`](https://statmodels7.github.io/distributions7/reference/to_link_scale.md).

## Usage

``` r
deriv5_scale(distrib, y, theta, res, scale)
```

## Arguments

- distrib:

  A distribution object.

- y, theta:

  The observations and the parameters, as given to the method.

- res:

  The fifth order on the parameter scale.

- scale:

  `"parameter"` or `"link"`.

## Value

`res` itself on the parameter scale, its link-scale composition
otherwise.

## See also

[`link_scale_lower_orders()`](https://statmodels7.github.io/distributions7/reference/link_scale_lower_orders.md),
[`to_link_scale()`](https://statmodels7.github.io/distributions7/reference/to_link_scale.md).
