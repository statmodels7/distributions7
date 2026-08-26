# Lognormal Log-CDF Hessian

Closed form. On the log scale the family is location-scale, so with \\z
= (\log q - \mu)/\sigma\\ and \\\varphi\\ the standard normal density,
\$\$\frac{\partial^2 F}{\partial\mu^2} = -\frac{z\varphi}{\sigma^2},
\qquad \frac{\partial^2 F}{\partial\mu\\\partial\sigma^2} =
\frac{\varphi(1-z^2)}{2\sigma^3}, \qquad \frac{\partial^2
F}{\partial(\sigma^2)^2} = \frac{\varphi z(3-z^2)}{4\sigma^4}.\$\$

## Arguments

- distrib:

  A `Lognormal1Distrib` object, from
  [`lognormal1_distrib()`](https://statmodels7.github.io/distributions7/reference/lognormal1_distrib.md).

- q:

  A numeric vector of quantiles. At or below zero every component is
  zero, the distribution function being flat there.

- theta:

  A named list with components `mu` (any real value) and `sigma2`
  (positive), each a numeric vector of length 1 or `n`.

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

## Details

The formulas are written in \\\varphi(z)\\, so the Jacobian factor \\q\\
of the gradient does not appear: at second order it would have to be
differentiated too, and expressing everything on the log scale avoids
that.

The three components vanish at the points where their polynomial in
\\z\\ does, so `mu_mu` is exactly zero at the median \\q = e^{\mu}\\. A
relative comparison against a numerical derivative is meaningless there
and an absolute one is what to use.

## Notation

\\\mu\\ is the mean of \\\log Y\\, \\\sigma^2 \> 0\\ its variance, \\z =
(\log q - \mu)/\sigma\\, \\\varphi\\ the standard normal density and
\\F\\ the distribution function of \\Y\\.

## See also

[`distrib_grad_cdf.Lognormal1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_cdf.Lognormal1Distrib.md)
for the first order;
[`distrib_hess_cdf.Gaussian1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_cdf.Gaussian1Distrib.md),
the family on the log scale;
[`lognormal1_distrib()`](https://statmodels7.github.io/distributions7/reference/lognormal1_distrib.md).

## Examples

``` r
d <- lognormal1_distrib()
th <- list(mu = 0, sigma2 = 1)
q <- c(0.5, 1, 3)

# Against a central difference of the cdf, on an absolute scale: mu_mu is
# exactly zero at the median q = exp(mu) = 1.
exact <- distrib_hess_cdf(d, q, th, log = FALSE)
fd <- numerical_cdf_deriv(d, q, th, order = 2)
max(abs(unlist(exact[names(fd)]) - unlist(fd)))
#> [1] 7.450581e-09
```
