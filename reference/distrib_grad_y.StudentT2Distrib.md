# Student t Derivatives in the Response, Standard Deviation

The derivatives of the log-density in the response to order four, and
the mixed derivatives of orders one and two in the response and one and
two in \\(\mu, \sigma, \nu)\\, each a closed form from its own compiled
kernel.

## Arguments

- distrib:

  A `StudentT2Distrib` object.

- y:

  A numeric vector of observations.

- theta:

  A named list with components `mu`, `sigma` and `nu`.

- scale:

  `"parameter"` or `"link"`, for the mixed derivatives.

- ...:

  Unused.

- threads:

  A single positive integer, how many threads the kernel may use.
  Defaults to `1L`.

## Value

A numeric vector for the derivatives in the response; a named list,
keyed by parameter or by parameter pair, for the mixed derivatives.

## See also

[`distrib_gradient.StudentT2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.StudentT2Distrib.md).

## Examples

``` r
d <- student_t2_distrib()
distrib_grad_y(d, c(-2, 0.5), list(mu = 1, sigma = 2, nu = 6))
#> [1] 0.8400000 0.2153846
```
