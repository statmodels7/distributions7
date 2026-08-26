# Weibull Log-CDF Derivatives

Closed form at every order from one to four, from the survival function
\\S = \exp\\-(q/\mu)^{\sigma}\\\\. Writing \\h = \sigma(\log q -
\log\mu)\\ the exponent is \\L = -e^{h}\\, so its partial derivatives
are \\-e^{h}\\ times the complete Bell polynomial in the partials of
\\h\\, and those are elementary.

## Arguments

- distrib:

  A `Weibull1Distrib` object, from
  [`weibull1_distrib()`](https://statmodels7.github.io/distributions7/reference/weibull1_distrib.md).

- q:

  A numeric vector of quantiles. Values at or below zero give
  derivatives of exactly zero.

- theta:

  A named list with components `mu` (the scale, positive) and `sigma`
  (the shape, positive), each a numeric vector of length 1 or `n`.

- lower.tail:

  Is the lower tail wanted? A single logical, `TRUE` by default.

- log:

  Are derivatives of the log probability wanted? A single logical,
  `TRUE` by default.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A named list of numeric vectors of the order the generic asked for,
keyed as
[`deriv_names(distrib@params, order)`](https://statmodels7.github.io/distributions7/reference/deriv_names.md):
two components for the gradient, three for the Hessian, four at order 3
and five at order 4.

## The partials of h

\\\partial^j h/\partial\mu^j = \sigma(-1)^j(j-1)!/\mu^j\\, the same
without the factor \\\sigma\\ when one index names the shape, and
exactly zero when two do: \\h\\ is linear in the shape. That is what
keeps the expansion short at the higher orders.

## An exact zero worth knowing about

At \\q = \mu\\ the exponent \\h\\ vanishes and so does \\\partial
h/\partial\sigma = \log q - \log\mu\\, so the shape component of the
gradient is exactly zero there. A relative comparison against a
numerical derivative at that point measures nothing; an absolute one is
what to use, and it puts the closed route within \\2.3\times10^{-11}\\
of a central difference of the cdf.

## Notation

\\\mu \> 0\\ is the scale, \\\sigma \> 0\\ the shape, \\h = \sigma(\log
q - \log\mu)\\, \\F\\ the distribution function and \\S = 1 - F\\ the
survival function. The mean is \\\mu\\\Gamma(1+1/\sigma)\\.

## See also

[`surv_cdf_deriv_k()`](https://statmodels7.github.io/distributions7/reference/surv_cdf_deriv_k.md)
for the identity;
[`distrib_grad_cdf.ExponentialDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_cdf.ExponentialDistrib.md),
the shape-1 case;
[`weibull1_distrib()`](https://statmodels7.github.io/distributions7/reference/weibull1_distrib.md).

## Examples

``` r
d <- weibull1_distrib()
q <- c(0.5, 2, 5)
th <- list(mu = 2, sigma = 3)

# Against a central difference of the cdf, on an absolute scale.
fd <- numerical_cdf_deriv(d, q, th, order = 1)
max(abs(unlist(distrib_grad_cdf(d, q, th, log = FALSE)) - unlist(fd)))
#> [1] 2.329359e-11

# The shape component is exactly zero at q = mu.
distrib_grad_cdf(d, 2, th, log = FALSE)$sigma
#> [1] 0
```
