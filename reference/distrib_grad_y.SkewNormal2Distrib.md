# Skew Normal Derivatives in the Response, Centered Parametrization

The first and second derivatives of the log-density in the response, and
the mixed derivatives of orders one and two in the response and one and
two in \\(\mu, \sigma, \gamma_1)\\, each from its own compiled kernel,
written as the parameter derivatives are (see
[`distrib_gradient.SkewNormal2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.SkewNormal2Distrib.md)).

## Arguments

- distrib:

  A `SkewNormal2Distrib` object.

- y:

  A numeric vector of observations.

- theta:

  A named list with components `mu`, `sigma` and `gamma1`.

- scale:

  `"parameter"` or `"link"`, for the mixed derivatives.

- ...:

  Unused.

- threads:

  A single positive integer, how many threads the kernel may use.
  Defaults to `1L`.

## Value

A numeric vector for the derivatives in the response; a named list for
the mixed derivatives.

## See also

[`distrib_gradient.SkewNormal2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.SkewNormal2Distrib.md).

## Examples

``` r
d <- skewnormal2_distrib()
distrib_grad_y(d, c(-1.7, 0.3), list(mu = 0.2, sigma = 1.3, gamma1 = 0.4))
#> [1]  1.425500 -0.243771
```
