# Elastic Net Log-CDF Derivatives

Closed form at every order from one to four. Each half of the
distribution function is a truncated Gaussian: with \\z = q - \mu\\, \\s
= \sqrt c\\ and \\x = a/\sqrt c\\ it is \\e^{w}\Phi(X)\\ below the
location and \\1 - e^{w}\Phi(X)\\ above it, for \\X = \pm sz - x\\ and a
weight \\w = -\log M(x) + x^2/2 + \mathrm{const}\\ written through the
Mills ratio the family already carries. The four registrations are made
together by
[`register_phi_terms_cdf()`](https://statmodels7.github.io/distributions7/reference/register_phi_terms_cdf.md).

## Arguments

- distrib:

  An `EnetDistrib` object, from
  [`enet_distrib()`](https://statmodels7.github.io/distributions7/reference/enet_distrib.md).

- q:

  A numeric vector of quantiles.

- theta:

  A named list with components `mu` (any real value), `lambda`
  (positive) and `alpha` (strictly between 0 and 1), each a numeric
  vector of length 1 or `n`.

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
three components for the gradient, six for the Hessian, ten at order 3
and fifteen at order 4.

## Why four orders are cheap here

Both \\s\\ and \\x\\ are a function of \\\lambda\\ times a function of
\\\alpha\\, so their mixed partials are products of one-variable ones
through
[`separable_deriv()`](https://statmodels7.github.io/distributions7/reference/separable_deriv.md),
and the one-variable ones are powers of affine functions through
[`dpow_affine()`](https://statmodels7.github.io/distributions7/reference/dpow_affine.md).
No multivariate expansion is formed at any order.

## The kink at the location

The location is the non-regular direction, as it is in the Laplace this
family contains: the second derivative in \\\mu\\ carries a point mass
at \\q = \mu\\, and the formulas hold on either side of it. A numerical
check that straddles the location is checking the arithmetic against a
reference that is not valid there.

## What the closed route is worth

Against a product stencil on the same cdf: \\1.6\times10^{-9}\\ at order
1, \\3.5\times10^{-7}\\ at order 2, \\3.0\times10^{-5}\\ at order 3 and
\\3.9\times10^{-4}\\ at order 4.

## Notation

\\\mu\\ is the location, \\\lambda \> 0\\ and \\\alpha \in (0,1)\\ the
two hyperparameters, \\a = \lambda\alpha\\ and \\c = \lambda(1-\alpha)\\
the two rates, \\M\\ the Mills ratio, \\\Phi\\ the standard normal
distribution function and \\F\\ the elastic net's.

## See also

[`phi_terms_cdf_deriv_k()`](https://statmodels7.github.io/distributions7/reference/phi_terms_cdf_deriv_k.md)
for the construction;
[`separable_deriv()`](https://statmodels7.github.io/distributions7/reference/separable_deriv.md)
and
[`dpow_affine()`](https://statmodels7.github.io/distributions7/reference/dpow_affine.md)
for the partials;
[`distrib_grad_cdf.InvGauss1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_cdf.InvGauss1Distrib.md),
the other family of this shape;
[`enet_distrib()`](https://statmodels7.github.io/distributions7/reference/enet_distrib.md).

## Examples

``` r
d <- enet_distrib()
q <- c(-1, 0.5, 2)
th <- list(mu = 0, lambda = 1, alpha = 0.5)

# Against a central difference of the cdf, which shares no arithmetic.
fd <- numerical_cdf_deriv(d, q, th, order = 1)
max(abs(unlist(distrib_grad_cdf(d, q, th, log = FALSE)) / unlist(fd) - 1))
#> [1] 1.580186e-09

# Three, six, ten and fifteen components as the order rises.
lengths(list(distrib_grad_cdf(d, q, th),
             distrib_hess_cdf(d, q, th),
             distrib_deriv3_cdf(d, q, th),
             distrib_deriv4_cdf(d, q, th)))
#> [1]  3  6 10 15
```
