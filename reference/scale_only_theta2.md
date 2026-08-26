# The Scale Component of a Family With No Location

Computes \\\partial^2\ell^{(y)}/\partial\sigma^2\\ or
\\\partial^2\ell^{(yy)}/\partial\sigma^2\\ for a family whose response
enters only through \\z = y/\sigma\\. It is
[`loc_scale_theta2_block()`](https://statmodels7.github.io/distributions7/reference/loc_scale_theta2_block.md)'s
scale component, and the same formula: the derivation never used the
location, only that \\\sigma\\ is a scale.

## Usage

``` r
scale_only_theta2(distrib, y, theta, order = 1L, at = 1L)
```

## Arguments

- distrib:

  An object inheriting from class `distrib`, with a scale and no
  location.

- y:

  A numeric vector of observations.

- theta:

  A named list of parameters.

- order:

  `1` for the block of
  [`distrib_grad_y()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.md),
  `2` for that of
  [`distrib_hess_y()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y.md).
  Default `1`.

- at:

  The index of the scale parameter in `distrib@params`. Default `1`.

## Value

A numeric vector as long as `y`: the scale's own pair alone, not a list.

## Details

With \\A = \sigma\ell^{(y)}\\ and \\B = \sigma^2\ell^{(yy)}\\,
\$\$\frac{\partial^2\ell^{(y)}}{\partial\sigma^2} = \frac{z^2A'' +
4zA' + 2A}{\sigma^3}, \qquad
\frac{\partial^2\ell^{(yy)}}{\partial\sigma^2} = \frac{z^2B'' + 6zB' +
6B}{\sigma^4},\$\$ and the primed quantities are the family's own third
and fourth response derivatives times a power of \\\sigma\\, so nothing
is differentiated here. Any parameter beyond the scale is a SHAPE and is
not covered;
[`scale_only_theta2_methods()`](https://statmodels7.github.io/distributions7/reference/scale_only_theta2_methods.md)
differences those.

## Notation

\\\ell\\ is the log-density of one observation, \\\ell^{(y)}\\ and
\\\ell^{(yy)}\\ its first and second derivatives in the response,
\\\sigma\\ the scale and \\z\\ the standardized response, \\y/\sigma\\
where there is no location and \\(y-\mu)/\sigma\\ where there is.

## See also

[`loc_scale_theta2_block()`](https://statmodels7.github.io/distributions7/reference/loc_scale_theta2_block.md),
the location-scale form,
[`scale_only_theta2_methods()`](https://statmodels7.github.io/distributions7/reference/scale_only_theta2_methods.md),
which splices it into the fallback, and
[`scale_only_cross2_method()`](https://statmodels7.github.io/distributions7/reference/scale_only_cross2_method.md)
for the first order in \\\theta\\.

## Examples

``` r
# The exponential has a scale and nothing else.
d <- exponential_distrib()
y <- c(0.4, 1.1, 2.3)
theta <- list(mu = 1.5)

distributions7:::scale_only_theta2(d, y, theta, 1L, 1L)
#> [1] -0.5925926 -0.5925926 -0.5925926

# The identity written out.
s <- 1.5
z <- y / s
A <- s * distrib_grad_y(d, y, theta)
A1 <- s^2 * distrib_hess_y(d, y, theta)
A2 <- s^3 * distrib_deriv3_y(d, y, theta)
(z^2 * A2 + 4 * z * A1 + 2 * A) / s^3
#> [1] -0.5925926 -0.5925926 -0.5925926

# And what the family's own method reports.
distrib_grad_y_hess(d, y, theta)$mu_mu
#> [1] -0.5925926 -0.5925926 -0.5925926
```
