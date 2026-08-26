# Lognormal Log-CDF Gradient

Closed form. On the log scale the lognormal is a location-scale family,
so \\F(q) = \Phi(z)\\ with \\z = (\log q - \mu)/\sigma\\, and the two
derivatives are \$\$\frac{\partial F}{\partial\mu} = -q\\f(q), \qquad
\frac{\partial F}{\partial\sigma^2} = -\frac{q\\f(q)\\z}{2\sigma}.\$\$

## Arguments

- distrib:

  A `Lognormal1Distrib` object, from
  [`lognormal1_distrib()`](https://statmodels7.github.io/distributions7/reference/lognormal1_distrib.md).

- q:

  A numeric vector of quantiles, positive. At or below zero the
  distribution function is zero and its derivatives are too.

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

A named list of two numeric vectors, `mu` and `sigma2`, each the length
of `q` recycled against `theta`.

## Details

The factor \\q\\ is the Jacobian of the change of variable: the standard
normal density at \\z\\ is \\\varphi(z) = q\sigma f(q)\\, so the
location-scale formulas are rewritten in the density of \\Y\\ itself and
nothing on the log scale has to be evaluated. The extra \\1/(2\sigma)\\
in the second component is the chain rule onto the variance, this
parametrization carrying \\\sigma^2\\ rather than \\\sigma\\.

Only the gradient is registered; the second derivatives fall to
[`distrib_hess_cdf.continuous_distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_cdf.continuous_distrib.md),
which differences the cdf and reuses this closed gradient for its
first-order part.

## Notation

\\\mu\\ is the mean of \\\log Y\\, \\\sigma^2 \> 0\\ its variance, \\z =
(\log q - \mu)/\sigma\\, \\f\\ the density of \\Y\\ and \\\varphi\\ the
standard normal density.

## See also

[`distrib_hess_cdf.continuous_distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_cdf.continuous_distrib.md),
the second order;
[`distrib_grad_cdf.Gaussian1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_cdf.Gaussian1Distrib.md),
the family on the log scale;
[`lognormal1_distrib()`](https://statmodels7.github.io/distributions7/reference/lognormal1_distrib.md).

## Examples

``` r
d <- lognormal1_distrib()
th <- list(mu = 0, sigma2 = 1)
q <- c(0.5, 1, 3)

# The mean component is -q f(q), the Jacobian factor included.
all.equal(distrib_grad_cdf(d, q, th, log = FALSE)$mu,
          -q * distrib_pdf(d, q, th))
#> [1] TRUE

# Against a central difference of the cdf, which shares no arithmetic.
fd <- numerical_cdf_deriv(d, q, th, order = 1)
max(abs(distrib_grad_cdf(d, q, th, log = FALSE)$sigma2 / fd$sigma2 - 1))
#> [1] NaN
```
