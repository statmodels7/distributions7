# Gaussian Log-CDF Gradient in Mean and Variance

Closed form, by the chain rule on the scale parametrization's
derivatives through \\\sigma = \sqrt{\sigma^2}\\. The mean component is
unchanged, \\-f(q)\\, and the variance component is the scale one
divided by \\2\sigma\\, which is the map's only non-zero partial.

## Arguments

- distrib:

  A `Gaussian2Distrib` object, from
  [`gaussian2_distrib()`](https://statmodels7.github.io/distributions7/reference/gaussian2_distrib.md).

- q:

  A numeric vector of quantiles.

- theta:

  A named list with components `mu` and `sigma2` (positive), each a
  numeric vector of length 1 or `n`.

- lower.tail:

  Is the lower tail wanted? A single logical, `TRUE` by default.

- log:

  Are derivatives of the log probability wanted? A single logical,
  `TRUE` by default.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A named list of two numeric vectors, `mu` and `sigma2`, each the length
of `q` recycled against `theta`.

## Notation

\\\mu\\ is the mean, \\\sigma^2 \> 0\\ the variance, \\f\\ the density
and \\F\\ the distribution function.

## See also

[`distrib_hess_cdf.Gaussian2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_cdf.Gaussian2Distrib.md)
for the second order;
[`distrib_grad_cdf.Gaussian1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_cdf.Gaussian1Distrib.md),
the parent;
[`mapped_cdf_deriv()`](https://statmodels7.github.io/distributions7/reference/mapped_cdf_deriv.md);
[`gaussian2_distrib()`](https://statmodels7.github.io/distributions7/reference/gaussian2_distrib.md).

## Examples

``` r
q <- c(-1, 0.5, 2)
g2 <- distrib_grad_cdf(gaussian2_distrib(), q,
                       list(mu = 0.3, sigma2 = 1.44), log = FALSE)
g1 <- distrib_grad_cdf(gaussian1_distrib(), q,
                       list(mu = 0.3, sigma = 1.2), log = FALSE)

# The mean component is unchanged by the map.
all.equal(g2$mu, g1$mu)
#> [1] TRUE

# The variance component is the scale one over 2 sigma.
all.equal(g2$sigma2, g1$sigma / (2 * 1.2))
#> [1] TRUE
```
