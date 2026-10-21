# Gaussian Log-CDF Gradient in Mean and Variance

Computes the derivatives in closed form, by the chain rule on the scale
parametrization's derivatives through \\\sigma = \sqrt{\sigma^2}\\. The
mean component is unchanged, \\-f(q)\\, and the variance component is
the scale one divided by \\2\sigma\\. This factor is the derivative of
\\\sigma\\ in \\\sigma^2\\, the only non-zero partial of the map.

## Arguments

- distrib:

  A `Gaussian2Distrib` object, from
  [`gaussian2_distrib()`](https://statmodels7.github.io/distributions7/reference/gaussian2_distrib.md).

- q:

  A numeric vector of quantiles.

- theta:

  A named list with components `mu` and `sigma2` (positive); each is a
  numeric vector of length 1 or `n`.

- lower.tail:

  A single logical, `TRUE` by default; `TRUE` selects the lower tail and
  `FALSE` the upper tail.

- log:

  A single logical, `TRUE` by default; `TRUE` selects the derivatives of
  the log probability and `FALSE` those of the probability.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A named list of two numeric vectors, `mu` and `sigma2`; each has the
length of `q` recycled against `theta`.

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
