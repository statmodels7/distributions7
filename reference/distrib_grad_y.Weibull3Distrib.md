# Weibull Derivatives in the Response, Mean and Shape

The derivatives of the log-density in the response to order four, and
the mixed derivatives of orders one and two in the response and one and
two in \\(m, \sigma)\\, each a closed form from its own compiled kernel.

## Arguments

- distrib:

  A `Weibull3Distrib` object.

- y:

  A numeric vector of strictly positive observations.

- theta:

  A named list with components `mean` and `sigma`.

- scale:

  `"parameter"` or `"link"`, for the mixed derivatives.

- ...:

  Unused.

## Value

A numeric vector for the derivatives in the response; a named list,
keyed by parameter or by parameter pair, for the mixed derivatives.

## See also

[`distrib_gradient.Weibull3Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.Weibull3Distrib.md).

## Examples

``` r
d <- weibull3_distrib()
distrib_grad_y(d, c(0.7, 2.5), list(mean = 4, sigma = 1.7))
#> [1] 0.89664341 0.02804236
```
