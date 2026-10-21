# Gaussian Log-CDF Hessian

Computes the Hessian in closed form, from the same location-scale
structure. With \\\ell_y = \partial\log f/\partial y\\, which for a
Gaussian is \\-z/\sigma\\, \$\$\frac{\partial^2 F}{\partial\mu^2} =
f\\\ell_y, \qquad \frac{\partial^2 F}{\partial\mu\\\partial\sigma} =
f\left(z\ell_y + \frac{1}{\sigma}\right), \qquad \frac{\partial^2
F}{\partial\sigma^2} = f\left(z^2\ell_y + \frac{2z}{\sigma}\right).\$\$
The method is
[`loc_scale_hess_cdf()`](https://statmodels7.github.io/distributions7/reference/loc_scale_hess_cdf.md)
itself, shared with the logistic, the Laplace and the Gumbel.

## Arguments

- distrib:

  A `Gaussian1Distrib` object, from
  [`gaussian1_distrib()`](https://statmodels7.github.io/distributions7/reference/gaussian1_distrib.md).

- q:

  A numeric vector of quantiles.

- theta:

  A named list with components `mu` and `sigma` (positive); each is a
  numeric vector of length 1 or `n`.

- lower.tail:

  A single logical, `TRUE` by default; `TRUE` selects the lower tail and
  `FALSE` the upper tail.

- log:

  A single logical, `TRUE` by default; `TRUE` selects the derivatives of
  the log probability and `FALSE` those of the probability.

## Value

A named list of three numeric vectors keyed as
[`hess_names()`](https://statmodels7.github.io/distributions7/reference/hess_names.md),
`mu_mu`, `sigma_sigma` and `mu_sigma`, in which each vector has the
length of `q` recycled against `theta`. The gradient is not returned
alongside.

## Notation

\\\mu\\ is the mean, \\\sigma \> 0\\ the standard deviation, \\z =
(q-\mu)/\sigma\\, \\f\\ the density and \\\ell_y = \partial\log
f/\partial y\\ its response derivative.

## See also

[`loc_scale_hess_cdf()`](https://statmodels7.github.io/distributions7/reference/loc_scale_hess_cdf.md)
for the shared body;
[`distrib_grad_cdf.Gaussian1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_cdf.Gaussian1Distrib.md)
for the first order;
[`gaussian1_distrib()`](https://statmodels7.github.io/distributions7/reference/gaussian1_distrib.md).

## Examples

``` r
d <- gaussian1_distrib()
th <- list(mu = 0.3, sigma = 1.2)
q <- c(-1, 0.5, 2)

distrib_hess_cdf(d, q, th)
#> $mu_mu
#> [1] -0.5627614 -0.4157630 -0.1735906
#> 
#> $sigma_sigma
#> [1] -3.05624806  0.14930692 -0.03617539
#> 
#> $mu_sigma
#> [1]  1.7154051  0.4132739 -0.1357278
#> 

# The result is compared with a central difference of the cdf, which shares
# no arithmetic with it.
exact <- distrib_hess_cdf(d, q, th, log = FALSE)
fd <- numerical_cdf_deriv(d, q, th, order = 2)
max(abs(unlist(exact[names(fd)]) / unlist(fd) - 1))
#> [1] 1.656387e-07
```
