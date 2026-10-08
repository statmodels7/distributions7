# Lognormal Derivatives in the Response, Mean and Variance of Y

The derivatives of the log-density in the response to order four, and
the mixed derivatives of orders one and two in the response and one and
two in \\(m, v)\\, each a closed form from its own compiled kernel.

## Arguments

- distrib:

  A `Lognormal2Distrib` object.

- y:

  A numeric vector of strictly positive observations.

- theta:

  A named list with components `mean` and `var`.

- scale:

  `"parameter"` or `"link"`, for the mixed derivatives.

- ...:

  Unused.

## Value

A numeric vector for the derivatives in the response; a named list,
keyed by parameter or by parameter pair, for the mixed derivatives.

## See also

[`distrib_gradient.Lognormal2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.Lognormal2Distrib.md).

## Examples

``` r
d <- lognormal2_distrib()
th <- list(mean = 3, var = 2)
distrib_grad_y(d, c(0.7, 2.5), th)
#> [1]  8.2173091 -0.2365756
```
