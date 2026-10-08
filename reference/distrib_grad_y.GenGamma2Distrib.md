# Generalized Gamma Derivatives in the Response, Mean

The derivatives of the log-density in the response to order four, and
the mixed derivatives of orders one and two in the response and one and
two in \\(m, d, p)\\, each a closed form from its own compiled kernel.

## Arguments

- distrib:

  A `GenGamma2Distrib` object.

- y:

  A numeric vector of strictly positive observations.

- theta:

  A named list with components `mean`, `d` and `p`.

- scale:

  `"parameter"` or `"link"`, for the mixed derivatives.

- ...:

  Unused.

## Value

A numeric vector for the derivatives in the response; a named list,
keyed by parameter or by parameter pair, for the mixed derivatives.

## See also

[`distrib_gradient.GenGamma2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.GenGamma2Distrib.md).

## Examples

``` r
d <- gengamma2_distrib()
distrib_grad_y(d, c(1.2, 4), list(mean = 5, d = 3, p = 1.5))
#> [1] 1.395430343 0.004792491
```
