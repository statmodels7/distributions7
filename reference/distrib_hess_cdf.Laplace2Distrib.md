# Laplace Log-CDF Hessian in Location and Rate

Closed form. Writing \\r = q - \mu\\ and \\a = \|r\|\\,
\$\$\frac{\partial^2 F}{\partial\mu^2} = -\lambda\\\mathrm{sign}(r)\\f,
\qquad \frac{\partial^2 F}{\partial\lambda^2} =
-\frac{\mathrm{sign}(r)\\a^2 f}{\lambda}, \qquad \frac{\partial^2
F}{\partial\mu\\\partial\lambda} = -\left(\frac{1}{\lambda} - a\right)
f.\$\$ The sign function is what carries the kink: all three jump at \\q
= \mu\\, as they do in the scale parametrization.

## Arguments

- distrib:

  A `Laplace2Distrib` object, from
  [`laplace2_distrib()`](https://statmodels7.github.io/distributions7/reference/laplace2_distrib.md).

- q:

  A numeric vector of quantiles. At \\q = \mu\\ exactly, `sign(0)` is 0
  and the components are the average of their two one-sided limits.

- theta:

  A named list with components `mu` and `lambda` (positive), each a
  numeric vector of length 1 or `n`.

- lower.tail:

  Is the lower tail wanted? A single logical, `TRUE` by default.

- log:

  Are derivatives of the log probability wanted? A single logical,
  `TRUE` by default.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A named list of three numeric vectors keyed as
[`hess_names()`](https://statmodels7.github.io/distributions7/reference/hess_names.md),
each the length of `q` recycled against `theta`.

## Notation

\\\mu\\ is the location, \\\lambda \> 0\\ the rate, \\r = q - \mu\\, \\a
= \|r\|\\, \\f\\ the density and \\F\\ the distribution function.

## See also

[`distrib_grad_cdf.Laplace2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_cdf.Laplace2Distrib.md)
for the first order;
[`distrib_hess_cdf.LaplaceDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_cdf.LaplaceDistrib.md)
for the scale parametrization;
[`laplace2_distrib()`](https://statmodels7.github.io/distributions7/reference/laplace2_distrib.md).

## Examples

``` r
d <- laplace2_distrib()
th <- list(mu = 0.3, lambda = 1 / 1.2)

# All three components jump across q = mu.
distrib_hess_cdf(d, 0.3 + c(-1e-6, 1e-6), th, log = FALSE)
#> $mu_mu
#> [1]  0.3472219 -0.3472219
#> 
#> $lambda_lambda
#> [1]  4.999996e-13 -4.999996e-13
#> 
#> $mu_lambda
#> [1] -0.4999992 -0.4999992
#> 

# Away from the kink they agree with a central difference of the cdf.
q <- c(-1, 2)
exact <- distrib_hess_cdf(d, q, th, log = FALSE)
fd <- numerical_cdf_deriv(d, q, th, order = 2)
max(abs(unlist(exact[names(fd)]) / unlist(fd) - 1))
#> [1] 2.17143e-07
```
